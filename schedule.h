#ifndef SCHEDULE
#define SCHEDULE

#include "mps_gen.h"

void dynamic_build(vector<boost::shared_ptr<SchmidtBasis>> basis);

void static_build(vector<boost::shared_ptr<SchmidtBasis>> basis);

#endif
