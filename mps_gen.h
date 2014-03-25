#ifndef MPSGEN
#define MPSGEN

#include "include.h"
#include "schmidt.h"

QSDArray<3, Quantum> generate_mps(boost::shared_ptr<SchmidtBasis> s1, boost::shared_ptr<SchmidtBasis> s2, bool additional = false);

void compute_dense(DArray<3>& d, int ql, int idx_p, int qr, boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, bool use_left);

class Overlap {
public:
  virtual double operator() (const vector<bool>& c1, const vector<bool>& c2) = 0;
};

class Overlap_Slater_Left: public Overlap {
protected:
  Matrix work_a, work_b;
  Matrix m_ca, m_ac, m_aa, m_sa;
  int total_a, total_b;
  int lcore, lactive_a, lactive_b, lactive_size;
  int rcore, ractive_a, ractive_b, ractive_size;
  int parity;
  bool this_site_up;
public:
  Overlap_Slater_Left(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, int ql, int qr);
  double operator() (const vector<bool>& c1, const vector<bool>& c2);
};

class Overlap_Slater_Right: public Overlap {
protected:
  Matrix work_a, work_b;
  Matrix m_ca, m_ac, m_aa, m_sa;
  int total_a, total_b;
  int lcore, lactive_a, lactive_b, lactive_size;
  int rcore, ractive_a, ractive_b, ractive_size;
  int parity;
  bool this_site_up;
public:
  Overlap_Slater_Right(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, int ql, int qr);  
  double operator() (const vector<bool>& c1, const vector<bool>& c2);  
};

//class Overlap_BCS: public Overlap {
//protected:
//  int ncore, ntotal;
//  Matrix worksheet;
//public:
//  virtual double operator() (vector<bool>& c1, vector<bool>& c2) = 0;
//};
//
//class Overlap_BCS_Left: public Overlap_BCS {
//protected:
//public:
//  double operator() (vector<bool>& c1, vector<bool>& c2);  
//};
//
//class Overlap_BCS_Right: public Overlap_BCS {
//protected:
//public:
//  double operator() (vector<bool>& c1, vector<bool>& c2);    
//};
#endif
