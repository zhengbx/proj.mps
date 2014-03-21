#ifndef MPSGEN
#define MPSGEN

#include "include.h"
#include "schmidt.h"

QSDArray<3, Quantum> generate_mps(boost::shared_ptr<SchmidtBasis> s1, boost::shared_ptr<SchmidtBasis> s2, bool additional = false);

#endif
