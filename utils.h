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
#include "newmat10/newmatap.h"
#include "newmat10/newmatio.h"
#include "newmat10/newmatutils.h"


using std::vector;
using std::string;

struct Input {
  double thr1p, thrnp;
  int M;
  bool calc_spectra, savemps, bcs;
  string temp_prefix, temp, path;

  Input() : thr1p(1e-7), thrnp(1e-8), M(0), calc_spectra(true), savemps(false) {
    temp_prefix = "/scratch/boxiao/MPSTemp";
  }
  friend std::ostream& operator << (std::ostream& os, const Input& inp);
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version) {
    ar & thr1p & thrnp & M & bcs;
    ar & calc_spectra & savemps;
    ar & temp_prefix & temp & path;
  }
};

// read configure file
void read_config(string file, Input& inp);

// read orbitals
Matrix read_orbitals(string file);

// permute sites
void permute(Matrix& orbs, vector<int>& order);

// create temporary storage
string mktmpdir(const string& prefix);

// print banner
void banner();

extern Input params;
extern Matrix coefs;

#endif
