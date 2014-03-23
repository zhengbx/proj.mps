#ifndef SCHEDULE
#define SCHEDULE

#include <boost/mpi.hpp>
#include "mps_gen.h"

namespace mpi = boost::mpi;

void dynamic_build(MPS<Quantum>& A, vector<boost::shared_ptr<SchmidtBasis>> basis);

void static_build(MPS<Quantum>& A, vector<boost::shared_ptr<SchmidtBasis>> basis);

#endif
