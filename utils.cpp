#include "utils.h"
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <cmath>

#define PI 3.14159265358979

using boost::trim;
using boost::is_any_of;
using std::cout;
using std::endl;
using std::ifstream;

Input params;
Matrix coefs;
vector<Matrix> k_coefs;
vector<int> orb_k;
vector<int> norb_k;

void permute(Matrix& orbs, vector<int>& order) {
  Matrix orbs_i = orbs;
  for (int i = 0; i < order.size(); ++i) {
    orbs.row(i) = orbs_i.row(order[i]-1);
  }
  if (order.size() == orbs.rows()) {
    return;
  } else if (order.size() * 2 == orbs.rows()) {
    for (int i = 0; i < order.size(); ++i) {
      orbs.row(order.size()+i) = orbs_i.row(order.size() + order[i]-1);
    }
  } else {
    cout << "wrong reorder vector" << endl;
    abort();
  }
}

dtype unity(d_real k) {
  // compute exp(2*\pi*j*k)
  d_real c = cos(2*PI*k);
  d_real s = sin(2*PI*k);
#ifdef _COMPLEX
  return dtype(c, s);
#else
  assert(fabs(s) < 1e-5);
  return c;
#endif
}

void read_config(string file, Input& inp) {
  string line;
  ifstream in(file.c_str());

  vector<string> tokens;
  while(std::getline(in, line)) {
    if (line.find('!') != string::npos) {
      boost::erase_tail(line, line.size() - line.find('!'));
    }
    trim(line);
    boost::split(tokens, line, is_any_of("\t ="), boost::token_compress_on);
    if (tokens[0].size() == 0) {
      continue;
    }
    if (boost::iequals(tokens[0], "thr1p")) {
      inp.thr1p = atof(tokens[1].c_str());
    } else if (boost::iequals(tokens[0], "thrnp")) {
      inp.thrnp = atof(tokens[1].c_str());
    } else if (boost::iequals(tokens[0], "M")) {
      inp.M = atoi(tokens[1].c_str());
    } else if (boost::iequals(tokens[0], "nospectra")) {
      inp.calc_spectra = false;
    } else if (boost::iequals(tokens[0], "savemps")) {
      inp.savemps = true;
    } else if (boost::iequals(tokens[0], "memtest")) {
      inp.mem_test = true;
    } else if (boost::iequals(tokens[0], "prefix")) {
      inp.temp_prefix = tokens[1];
    } 
#ifdef _COMPLEX    
    else if (boost::iequals(tokens[0], "kpoints")) {
      if (boost::iequals(tokens[1], "pbc") || boost::iequals(tokens[1], "apbc")) {
        int nkpts = atoi(tokens[2].c_str());
        d_real base = 1./(d_real)nkpts;
        d_real shift = boost::iequals(tokens[1], "apbc") ? 0.5*base:0;
        for (int i = 0; i < nkpts; ++i) {
          inp.kpoints.push_back(base*i+shift);
        }
      } else {
        // set kpoints manually
        for (int i = 1; i < tokens.size(); ++i) {
          inp.kpoints.push_back((d_real)atof(tokens[i].c_str()));
        }
      }
    } else if (boost::iequals(tokens[0], "use_kpoints")) {
      for (int i = 1; i < tokens.size(); ++i) {
        inp.use_k.push_back(atoi(tokens[i].c_str()));
      }
    } 
#endif
    else {
      cout << "\nUnrecognized option in config file:" << endl;
      cout << "\t" << tokens[0] << endl;
      abort();
    }
  }
  // check sanity
  if ((inp.kpoints.size() == 0) != (inp.use_k.size() == 0)) {
    cout << "\nk-points and sites using k-points must be set at the same time" << endl;
    abort();
  }
  for (int i = 0; i < inp.use_k.size(); ++i) {
    if (inp.use_k[i] % inp.kpoints.size() != 0) {
      cout << "\ninvalid position with translational symmetry" << endl;
      abort();      
    }
  }
  inp.kspace = (inp.kpoints.size() != 0);
}

void read_orbitals(string file) {
  string line;
  ifstream in(file.c_str());
  std::getline(in, line);  // the first line is an explanation of the file.
    
  int nsites, norbs;
  vector<string> tokens;
  // nsites/norbs line
  std::getline(in, line);  
  trim(line);
  boost::split(tokens, line, is_any_of("\t "), boost::token_compress_on);
  if (tokens.size() == 1) {
    params.bcs = true;
    nsites = atoi(tokens[0].c_str());
  } else if (tokens.size() == 2) {
    params.bcs = false;
    nsites = atoi(tokens[0].c_str());
    norbs = atoi(tokens[1].c_str());
  } else {
    cout << "nsites and norbs line error" << endl;
    abort();
  }
  // reorder line
  vector<int> order;
  std::getline(in, line);
  trim(line);
  boost::split(tokens, line, is_any_of("\t "), boost::token_compress_on);
  if (tokens.size() == 1 && boost::iequals(tokens[0], "default")) {
    for (int i = 0; i < nsites; ++i) {
      order.push_back(i+1);
    }
  } else if (tokens.size() == nsites) {
    for (int i = 0; i < nsites; ++i) {
      order.push_back(atoi(tokens[i].c_str()));
    }
  } else {
    cout << "order line error" << endl;
    abort();
  }

  int ncol = params.bcs ? nsites : norbs;
  int nrow = params.bcs ? 2*nsites : nsites;
  coefs.resize(nrow, ncol);
  d_real temp;
  for (int i = 0; i < ncol; i++) {
    for (int j = 0; j < nrow; j++) {
      in >> temp;
      coefs(j, i) = temp;
    }
  }
  permute(coefs, order);
}

void read_orbitals_kspace(string file) {
  string line;
  ifstream in(file.c_str());
  std::getline(in, line);  // the first line is an explanation of the file.
  int nsites, norbs;
  vector<string> tokens;
  // nsites/norbs line
  std::getline(in, line);  
  trim(line);
  boost::split(tokens, line, is_any_of("\t "), boost::token_compress_on);
  if (tokens.size() == 1) {
    params.bcs = true;
    nsites = atoi(tokens[0].c_str());
  } else if (tokens.size() == 2) {
    params.bcs = false;
    nsites = atoi(tokens[0].c_str());
    norbs = atoi(tokens[1].c_str());
  } else {
    cout << "nsites and norbs line error" << endl;
    abort();
  }
  assert(nsites % params.kpoints.size() == 0); // make sure the k-space symmetry is consistent with number of sites
  int ncsites = nsites / params.kpoints.size(); // number of sites in a cell
  
  // reorder line
  vector<int> order;
  std::getline(in, line);
  trim(line);
  boost::split(tokens, line, is_any_of("\t "), boost::token_compress_on);
  if (tokens.size() == 1 && boost::iequals(tokens[0], "default")) {
    for (int i = 0; i < ncsites; ++i) {
      for (int t = 0; t < params.kpoints.size(); ++t) {
        order.push_back(i+ncsites*t+1); // generic order for translational invariant wavefunction
      }
    }
  } else if (tokens.size() == nsites) {
    for (int i = 0; i < nsites; ++i) {
      order.push_back(atoi(tokens[i].c_str()));
    }
  } else {
    cout << "order line error" << endl;
    abort();
  }
  // eligibility of reorder vector
  for (int i = 0; i < params.use_k.size(); ++i) {
    int sites_per_cell = params.use_k[i] / params.kpoints.size();
    for (int t = 0; t < params.kpoints.size(); ++t) {
      for (int j = 0; j < sites_per_cell; ++j) {
        assert(std::find(order.begin(), order.begin()+params.use_k[i], j+t*ncsites+1) != order.begin()+params.use_k[i]);
      }
    }
  }

  int nrow = params.bcs ? 2*ncsites : ncsites;
  int ncols = 0;
  k_coefs.resize(params.kpoints.size());
  d_real temp_r, temp_i;

  for (int t = 0; t < params.kpoints.size(); ++t) { // read orbitals for each k-point
    // format:
    // norbs-for-this-k-point
    // orb1 (real imag real imag ...)
    // orb2
    // ...
    int ncol;
    in >> ncol;
    ncols += ncol;
    norb_k.push_back(ncol);
    k_coefs[t].resize(nrow, ncol);
    for (int i = 0; i < ncol; i++) {
      for (int j = 0; j < nrow; j++) {
#ifdef _COMPLEX
        in >> temp_r >> temp_i;
        k_coefs[t](j, i) = dtype(temp_r, temp_i);
#else
        in >> temp_r;
        k_coefs[t](j, i) = dtype(temp_r);
#endif
      }
      orb_k.push_back(t);
    }
  }

  assert(ncols == params.bcs ? nsites:norbs);
  // now build real space coef matrix
  coefs.resize(nrow*params.kpoints.size(), ncols);
  int offset = 0;
  for (int t = 0; t < params.kpoints.size(); ++t) {
    for (int s = 0; s < params.kpoints.size(); ++s) {
      coefs.block(s*ncsites, offset, ncsites, norb_k[t]) = k_coefs[t].block(0, 0, ncsites, norb_k[t]) * unity(params.kpoints[t]*s) / sqrt(params.kpoints.size());
      if (params.bcs) {
        coefs.block(nsites + s*ncsites, offset, ncsites, norb_k[t]) = k_coefs[t].block(ncsites, 0, ncsites, norb_k[t]) * unity(params.kpoints[t]) / sqrt(params.kpoints.size());
      }
    }
    offset += norb_k[t];
  }
  permute(coefs, order);
}

string mktmpdir(const string& prefix) {
  char* temp = new char[prefix.size() + 11];
  std::strcpy(temp, (prefix + "/tmpXXXXXX").c_str());
  mkdtemp(temp);
  // create this folder
  cout << "MPS Temporary Directory " << temp << endl;
  return string(temp);
}


void banner() {
  cout << "-----------------------------------------------------------------------\n";
  cout << "                           G P S - M P S                               \n";
  cout << "   (Gutzwiller Projection of Single Slater Determinants through MPS)   \n\n";
  cout << "                           Bo-Xiao Zheng                               \n";
  cout << "-----------------------------------------------------------------------\n\n";
}

std::ostream& operator << (std::ostream& os, const Input& inp) {
  os << "Projected BCS wfn       = " << (inp.bcs ? "True":"False") << std::endl;
  os << "one particle threshold  = " << inp.thr1p << std::endl;
  os << "many-body wfn threshold = " << inp.thrnp << std::endl;
  os << "limit M                 = " << inp.M << std::endl;
  os << "entanglement spectra    = " << (inp.calc_spectra ? "True":"False") << std::endl;
  os << "save MPS                = " << (inp.savemps ? "True":"False") << std::endl;
  os << "translational symm      = " << (inp.kspace ? "True":"False") << std::endl;
}
