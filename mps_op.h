#ifndef MPS_OP
#define MPS_OP

#include "SpinQuantum.h"
#include "include.h"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <tuple>

using std::map;
using std::tuple;

size_t mps_size(const QSDArray<3, Quantum>& A);

void check_existence(int i, const char* filename);

void save_site(const QSDArray<3, Quantum>&, int, const char *);
void load_site(QSDArray<3, Quantum>&, int, const char *);

void save_site(const MPS<Quantum>&, int, const char *);

void load_site(MPS<Quantum>&, int,const char *);

double norm_on_disk(MPS<Quantum>&, const char*);

void normalize_on_disk(MPS<Quantum>&, const char*, int site = -1);

void partial_compress(int L,const MPS_DIRECTION& dir, int D, const char* filename, int last);


void compress_on_disk(int L,const MPS_DIRECTION& dir, int D, const char* filename, bool store = false, bool on_the_fly = false);

// process entaglement spectra
void spectra(const tuple<SDArray<1>, Qshapes<Quantum>>& raw);

tuple<SDArray<1>, Qshapes<Quantum>> Schmidt_on_disk(MPS<Quantum>& mps, int site, const char* filename, int lc, int rc);
#endif
