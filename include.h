#ifndef INCLUDE_H_
#define INCLUDE_H_

#ifdef _COMPLEX
#include <complex>
typedef std::complex<double> dtype;
typedef double d_real;
#else
typedef double dtype;
typedef double d_real;
#endif

#include "SpinQuantum.h"
namespace btas { typedef SpinQuantum Quantum; };

//#include <boost/archive/text_oarchive.hpp>
//#include <boost/archive/text_iarchive.hpp>
//
//#include <btas/common/blas_cxx_interface.h>
//
//#include <btas/common/TVector.h>
//
//#include <btas/DENSE/TArray.h>
//#include <btas/SPARSE/STConj.h>
//#include <btas/QSPARSE/QSTArray.h>
//#include <btas/QSPARSE/QSTCONTRACT.h>

#include <MPSblas.h>

#include "SpinHamiltonian.h"

#endif
