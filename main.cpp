#include <iostream>
#include <boost/mpi.hpp>
#include "timer.h"

namespace mpi = boost::mpi;

#include "densitymat.h"
#include "mps_gen.h"

using std::cout;
using std::endl;

int main(int argc, char* argv[]){
  mpi::environment env(argc, argv);  
  cout.setf(std::ios::fixed, std::ios::floatfield);
  cout.precision(10);
  if (argc <= 1) {
    cout << "No input file specified" << endl;
    abort();
  }
  
  mpi::communicator world;  
  if (world.rank() == 0) {
    banner();
    params.path = string(argv[1]);
    params.temp = mktmpdir(params.temp_prefix);
    // read configure file
    read_config(params.path + "/config.in", params);
    coefs = read_orbitals(params.path + "/orbitals.in");
  }

  broadcast(world, params, 0);
  broadcast(world, coefs, 0);

  if (world.rank() == 0) {
    cout << "Calculation parameters" << endl;
    cout << params << endl;
  }
  
  auto dm = params.bcs ? 
    boost::shared_ptr<DensityMatrix>(new BCSDM(coefs)) : 
    boost::shared_ptr<DensityMatrix>(new SlaterDM(coefs));
  int nsites = dm -> get_nsites();
  vector<boost::shared_ptr<SchmidtBasis>> basis_set(nsites+1, nullptr);

  world.barrier();  
  Timer t_basis("Build Basis");
  t_basis.start();
  for (int i = 0; i < nsites+1; ++i) {
    if (i % world.size() == world.rank()) {
      cout << "Cut = " << i << " On processor " << world.rank() << endl;
      basis_set[i] = dm -> basis(i);
    }
  }
  t_basis.pause();
  
  world.barrier();
  if (world.rank() == 0) {
    cout << endl << "Build Schmidt basis" << endl;
    cout << "Processor 0 : " << t_basis.time() << " s" << endl;
    for (int i = 1; i < world.size(); ++i) {
      double time;
      world.recv(i, i, time);
      cout << "Processor " << i <<" : " << t_basis.time() << " s" << endl;
    }
    cout << endl;
    cout << "---------------------------" << endl;
    cout << "   Schmidt Basis Summary" << endl;
    cout << "---------------------------" << endl;    
    cout << endl;
  } else {
    world.send(0, world.rank(), t_basis.time());
  }

  Timer t_send("Broadcast Basis");
  if (world.rank() == 0) {
    t_send.start();
  }

  for (int i = 0; i < nsites+1; ++i) {
    broadcast(world, basis_set[i], i % world.size());
    if (world.rank() == 0) {
      cout << *basis_set[i] << endl;
    }
  }
  if (world.rank() == 0) {
    t_send.pause();
    t_send.print();
    cout << endl;
  }

  MPS<Quantum> A(nsites);
  for (int i = 0; i < nsites; ++i) {
    if (i % world.size() == world.rank()) {
      A[i] = generate_mps(basis_set[i], basis_set[i+1], i == nsites/2);
    }
  }
  return 0;
}
