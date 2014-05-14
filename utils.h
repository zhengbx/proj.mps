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

using std::vector;
using std::string;

typedef Eigen::Matrix<dtype, Eigen::Dynamic, Eigen::Dynamic> Matrix;

struct Input {
  d_real thr1p, thrnp;
  int M;
  bool calc_spectra, savemps, bcs, mem_test, kspace;
  string temp_prefix, temp, path; // temp is temporary dir, path is input file path
  vector<d_real> kpoints;
  vector<int> use_k;

  Input():
    thr1p(1e-7),
    thrnp(1e-8),
    M(0),
    calc_spectra(true),
    savemps(false),
    mem_test(false),
    temp_prefix("/scratch/boxiao/MPSTemp") {}
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
dtype unity(d_real k);

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
#endif
