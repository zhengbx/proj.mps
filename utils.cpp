#include "utils.h"
#include <boost/algorithm/string.hpp>

using boost::trim;
using boost::is_any_of;
using std::cout;
using std::endl;

Input params;
Matrix coefs;

void permute(Matrix& orbs, vector<int>& order) {
  Matrix orbs_i = orbs;
  for (int i = 0; i < order.size(); ++i) {
    orbs.Row(i+1) = orbs_i.Row(order[i]);
  }
  if (order.size() == orbs.Nrows()) {
    return;
  } else if (order.size() * 2 == orbs.Nrows()) {
    for (int i = 0; i < order.size(); ++i) {
      orbs.Row(order.size()+i+1) = orbs_i.Row(order.size() + order[i]);
    }
  } else {
    cout << "wrong reorder vector" << endl;
    abort();
  }
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
    } else if (boost::iequals(tokens[0], "prefix")) {
      inp.temp_prefix = tokens[1];
    } else {
      cout << "\nUnrecognized option in config file:" << endl;
      cout << "\t" << tokens[0] << endl;
      abort();
    }
  }
}

Matrix read_orbitals(string file) {
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
  Matrix orbs(nrow, ncol);
  double temp;
  for (int i = 0; i < ncol; i++) {
    for (int j = 0; j < nrow; j++) {
      in >> temp;
      orbs(j+1, i+1) = temp;
    }
  }
  permute(orbs, order);
  return std::move(orbs);
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
}
