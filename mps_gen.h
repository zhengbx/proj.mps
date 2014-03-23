#ifndef MPSGEN
#define MPSGEN

#include "include.h"
#include "schmidt.h"

QSDArray<3, Quantum> generate_mps(boost::shared_ptr<SchmidtBasis> s1, boost::shared_ptr<SchmidtBasis> s2, bool additional = false);


void compute_dense(DArray<3> d, int ql, int idx_p, int qr, boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, bool use_left);
#endif
