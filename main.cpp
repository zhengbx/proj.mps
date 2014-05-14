#include <iostream>
#include <boost/mpi.hpp>
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#include "timer.h"
#include "mps_op.h"

namespace mpi = boost::mpi;

#include "densitymat.h"
#include "schedule.h"


using std::cout;
using std::endl;

int main(int argc, char* argv[]){
  mpi::environment env(argc, argv);
  cout.setf(std::ios::fixed, std::ios::floatfield);
  cout.precision(4);
  if (argc <= 1) {
    cout << "No input file specified" << endl;
    abort();
  }
  
  mpi::communicator world;
  // the program cannot run with fewer than 3 cores
  if (world.size() < 3) {
    cout << "This program runs with at least 3 cores" << endl;
    abort();
  }

  // read configure file and orbital file
  if (world.rank() == 0) {
    banner();
    params.path = string(argv[1]);
    params.temp = mktmpdir(params.temp_prefix);
    // read configure file
    read_config(params.path + "/config.in", params);
    if (params.kspace) {
      read_orbitals_kspace(params.path + "/orbitals.in");
    } else {
      read_orbitals(params.path + "/orbitals.in");
    }
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
    //if (world.rank() == 0) {
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

  dynamic_build(basis_set);

  if (world.rank() == 0) {
    boost::filesystem::path mps_tmp_store(params.temp);
    if (!params.mem_test) {
      MPS<dtype, Quantum> A(nsites);
      //compress_on_disk(A, MPS_DIRECTION::Right, params.M, params.temp.c_str(), true);

      if (params.calc_spectra) {
        cout << "\nnow calculate entanglement spectrum\n";
        auto raw_spectra = Schmidt_on_disk(nsites, -1, params.temp.c_str(), nsites/2-2, nsites/2-2);
        spectra(raw_spectra);
      }

      if (params.savemps) {
        boost::filesystem::path mps_store(params.path + "/mps.out");
        boost::filesystem::create_directory(mps_store);
        for (int i = 0; i < nsites; ++i) {
          string filename = std::to_string(i) + ".mps";
          boost::filesystem::path p1(params.temp + "/" + filename);
          boost::filesystem::path p2(params.path + "/mps.out/" + filename);
          boost::filesystem::copy_file(p1, p2, boost::filesystem::copy_option::overwrite_if_exists);
        }
      }
    }
    boost::filesystem::remove_all(mps_tmp_store);
  }

  return 0;
}
