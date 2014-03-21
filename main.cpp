#include <iostream>
#include <boost/mpi.hpp>

namespace mpi = boost::mpi;

#include "densitymat.h"

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
  vector<boost::shared_ptr<SchmidtBasis>> basis_set(nsites, nullptr);

  for (int i = 0; i < nsites; ++i) {
    if (i % world.size() == world.rank()) {
      cout << "Cut = " << i << " On processor " << world.rank() << endl;
      basis_set[i] = dm -> basis(i);
    }
  }

  for (int i = 0; i < nsites; ++i) {
    broadcast(world, basis_set[i], i % world.size());
    if (world.rank() == 0) {
      cout << *basis_set[i] << endl;
    }   
  }

  return 0;
}
