#include "schedule.h"
#include "mps_op.h"
#include <boost/mpi.hpp>

namespace mpi = boost::mpi;

void dynamic_build(vector<boost::shared_ptr<SchmidtBasis>> basis) {
  int nsites = basis.size()-1;
  mpi::communicator world;  
  if (world.rank() == 0) {
    for (int i = 0; i <= nsites/2; ++i) {
      int rank_ready;
      world.recv(mpi::any_source, nsites+1, rank_ready);
      //printf("Master: Generating MPS of Site %3d On processor %2d\n", i, rank_ready);
      world.send(rank_ready, -1, i);
      world.send(rank_ready, i, basis[i]);
      world.send(rank_ready, i+1, basis[i+1]);
      basis[i].reset();
    }
    for (int i = nsites-1; i > nsites/2; --i) {
      int rank_ready;
      world.recv(mpi::any_source, nsites+1, rank_ready);
      //printf("Master: Generating MPS of Site %3d On processor %2d\n", i, rank_ready);
      world.send(rank_ready, -1, i);
      world.send(rank_ready, i, basis[i]);
      world.send(rank_ready, i+1, basis[i+1]);
      basis[i+1].reset();
    }
    basis[nsites/2+1].reset();
    for (int i = 2; i < world.size(); ++i) {
      world.send(i, -1, -1);
    }
  } else if (world.rank() == 1) {
    partial_compress(nsites, MPS_DIRECTION::Left, params.M, params.temp.c_str(), nsites/2);
    partial_compress(nsites, MPS_DIRECTION::Right, params.M, params.temp.c_str(), nsites/2-2);
    MPS<Quantum> A(nsites);
    cout << "Now normalize" << endl;
    normalize_on_disk(A, params.temp.c_str(), nsites/2-2);
    //compress_on_disk(nsites, MPS_DIRECTION::Right, params.M, params.temp.c_str(), true, true);
  } else {
    int do_site;
    world.send(0, nsites+1, world.rank());
    world.recv(0, -1, do_site);
    while (do_site != -1) {
      world.recv(0, do_site, basis[do_site]);
      world.recv(0, do_site+1, basis[do_site+1]);
      QSDArray<3, Quantum> A = generate_mps(basis[do_site], basis[do_site+1], do_site == nsites/2);
      save_site(A, do_site, params.temp.c_str());
      basis[do_site].reset();
      basis[do_site+1].reset();
      world.send(0, nsites+1, world.rank());
      world.recv(0, -1, do_site);
    }
  }
  world.barrier();
}

void static_build(vector<boost::shared_ptr<SchmidtBasis>> basis) {
  int nsites = basis.size()-1;
  mpi::communicator world;  
  for (int i = 0; i < nsites; ++i) {
    if (i % world.size() == world.rank()) {
      printf("Generating MPS of Site %3d On processor %2d", i, world.rank());
      QSDArray<3, Quantum> A = generate_mps(basis[i], basis[i+1], i == nsites/2);
      A.clear();
    }
  }
}
