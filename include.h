#ifndef INCLUDE_H_
#define INCLUDE_H_

#include <complex>

#ifdef _SINGLE
typedef std::complex<float> d_comp;
typedef float d_real;
#else
typedef std::complex<double> d_comp;
typedef double d_real;
#endif


#ifdef _COMPLEX
#include <complex>
typedef d_comp dtype;
#else
typedef d_real dtype;
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
