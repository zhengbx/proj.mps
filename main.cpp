#include <iostream>
#include <boost/mpi.hpp>
#include "timer.h"

namespace mpi = boost::mpi;

#include "densitymat.h"
#include "schedule.h"


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
  if (world.size() < 2) {
    cout << "This program runs with at least 2 cores" << endl;
    abort();
  }
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

  Timer t_send("Gathering Basis");
  if (world.rank() == 0) {
    t_send.start();
  }

  if (world.rank() == 0) {
    for (int i = 0; i < nsites+1; ++i) {
      if (i % world.size() != 0) {
        world.recv(i % world.size(), i, basis_set[i]);
      }
      cout << *basis_set[i] << endl;      
    }
  } else {
    for (int i = 0; i < nsites+1; ++i) {
      if (i % world.size() == world.rank()) {
        world.send(0, i, basis_set[i]);
        basis_set[i].reset();
      }
    }
  }

  if (world.rank() == 0) {
    t_send.pause();
    t_send.print();
    cout << endl;
  }

  MPS<Quantum> A(nsites);

  dynamic_build(basis_set);
  return 0;
}
