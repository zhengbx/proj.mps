#ifndef MPS_OP
#define MPS_OP

#include "SpinQuantum.h"
#include "include.h"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <tuple>

using std::map;
using std::tuple;

size_t mps_size(const QSTArray<dtype, 3, Quantum>& A);

void check_existence(int i, const char* filename);

void save_site(const QSTArray<dtype, 3, Quantum>&, int, const char *);
void load_site(QSTArray<dtype, 3, Quantum>&, int, const char *);

void save_site(const MPS<dtype, Quantum>&, int, const char *);

void load_site(MPS<dtype, Quantum>&, int,const char *);

double norm_on_disk(const char*, int size);

void normalize_on_disk(const char* filename, int size, int site = -1);

typename remove_complex<dtype>::type partial_compress(int L,const MPS_DIRECTION& dir, int D, const char* filename, int first, int last);

typename remove_complex<dtype>::type compress_on_disk(int L,const MPS_DIRECTION& dir, int D, const char* filename, bool store = false, bool incomp = false);

// process entaglement spectra
void spectra(const tuple<STArray<typename remove_complex<dtype>::type, 1>, Qshapes<Quantum>>& raw);

tuple<STArray<typename remove_complex<dtype>::type, 1>, Qshapes<Quantum>> Schmidt_on_disk(int L, int site, const char* filename, int lc, int rc);
#endif
