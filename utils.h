#ifndef UTILS
#define UTILS

#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <boost/shared_ptr.hpp>

#define BOOST_NO_CXX11_SCOPED_ENUMS
//#include "newmat10/newmatap.h"
//#include "newmat10/newmatio.h"
//#include "newmat10/newmatutils.h"

#ifndef EIGEN_CONFIG_H_
#define EIGEN_CONFIG_H_
#include <boost/serialization/array.hpp>
#define EIGEN_DENSEBASE_PLUGIN "../plugins/EigenDenseBaseAddons.h"

#include <Eigen/Dense>
#endif // EIGEN_CONFIG_H_

#include "include.h"
#define PI 3.14159265358979

using std::vector;
using std::string;

typedef Eigen::Matrix<dtype, Eigen::Dynamic, Eigen::Dynamic> Matrix;

class kpoint {
private:
  int N, n;
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version) {
    ar & N & n;
  }
  void init(const kpoint& other) {
    if (N == 0) {
      N = other.get_N();
      n = 0;
    } else if (N == other.get_N()) {
      return;
    } else {
      abort();
    }
  }
public:
  kpoint(int _N = 0, int _n = 0): N(_N), n(_n) {}
  kpoint operator+(const kpoint& other) const {
    
    assert(N == other.N);
    return kpoint(N, (n+other.n) % N);
  }
  void operator+=(const kpoint& other) {
    init(other);
    n = (n + other.n) % N;
  }
  bool operator==(const kpoint& other) const {
    return (N == other.N) && (n == other.n);
  }
  kpoint operator*(const int& s) const {
    return kpoint(N, (n*s) % N);
  }
  int get_N() const {  return N;}
  int get_n() const {  return n;}
  friend std::ostream& operator << (std::ostream& os, const kpoint& p) {
    os << "(" << p.N << ", " << p.n << ")";
    return os;
  }
  dtype eval() const {
    d_real theta = (double)n / N;
    d_real c = cos(2*PI*theta);
    d_real s = sin(2*PI*theta);
#ifdef _COMPLEX
    return dtype(c, s);
#else
    assert(fabs(s) < 1e-5);
    return c;
#endif
  }
};

struct Input {
  d_real thr1p, thrnp;
  int M;
  bool calc_spectra, savemps, bcs, mem_test, kspace;
  bool restart;
  string temp_prefix, temp, path; // temp is temporary dir, path is input file path
  vector<kpoint> kpoints;
  vector<int> use_k;

  Input():
    thr1p(1e-7),
    thrnp(1e-8),
    M(0),
    calc_spectra(true),
    savemps(false),
    mem_test(false),
    temp_prefix("/scratch/boxiao/MPSTemp"),
    restart(false) {}
  friend std::ostream& operator << (std::ostream& os, const Input& inp);
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version) {
    ar & thr1p & thrnp & M & bcs;
    ar & calc_spectra & savemps & mem_test;
    ar & temp_prefix & temp & path;
    ar & kpoints & use_k & kspace;
  }
};

// read configure file
void read_config(string file, Input& inp);

// read orbitals
void read_orbitals(string file);
void read_orbitals_kspace(string file);


// permute sites
void permute(Matrix& orbs, vector<int>& order);

// create temporary storage
string mktmpdir(const string& prefix);

// print banner
void banner();

extern Input params;
extern Matrix coefs; // orbital coef matrix
extern vector<Matrix> k_coefs; // orbital coef for each k point
extern vector<int> orb_k; // indicating which orbital belongs to which orbital
extern vector<int> norb_k; // number of orbitals for each k point
extern vector<int> order; // order of orbitals
#endif
