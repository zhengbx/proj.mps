#ifndef MPSGEN
#define MPSGEN

#include "include.h"
#include "schmidt.h"

QSTArray<dtype, 3, Quantum> generate_mps(boost::shared_ptr<SchmidtBasis> s1, boost::shared_ptr<SchmidtBasis> s2, bool additional = false);

void compute_dense(DArray<3>& d, int ql, int idx_p, int qr, boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, bool use_left, bool additional = false);


int common_parity(boost::shared_ptr<SchmidtBasis> s, boost::shared_ptr<ActiveSpaceIterator> it);

int individual_parity(const vector<bool>& bits);

class Overlap {
protected:
  virtual void build_cc_block(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr) = 0;
  virtual void build_ac_block(const vector<bool>& c1, const vector<bool>& c2) = 0;
  virtual void build_ca_block(const vector<bool>& c1, const vector<bool>& c2) = 0;
  virtual void build_aa_block(const vector<bool>& c1, const vector<bool>& c2) = 0;
public:
  virtual dtype operator() (const vector<bool>& c1, const vector<bool>& c2) = 0;
};

class Overlap_Slater: public Overlap {
protected:
  Matrix work_a, work_b;
  Matrix m_ca, m_ac, m_aa, m_sa;
  int total_a, total_b, nsites;
  int lcore, lactive_a, lactive_b, lactive_size;
  int rcore, ractive_a, ractive_b, ractive_size;
  int parity;
  bool this_site_up;
public:
  Overlap_Slater() {};
  dtype operator() (const vector<bool>& c1, const vector<bool>& c2);
};

class Overlap_Slater_Left: public Overlap_Slater {
protected:
  void build_cc_block(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr);
  void build_ac_block(const vector<bool>& c1, const vector<bool>& c2);
  void build_ca_block(const vector<bool>& c1, const vector<bool>& c2);
  void build_aa_block(const vector<bool>& c1, const vector<bool>& c2);
public:
  Overlap_Slater_Left(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, int ql, int qr);
};

class Overlap_Slater_Right: public Overlap_Slater {
protected:
  void build_cc_block(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr);
  void build_ac_block(const vector<bool>& c1, const vector<bool>& c2);
  void build_ca_block(const vector<bool>& c1, const vector<bool>& c2);
  void build_aa_block(const vector<bool>& c1, const vector<bool>& c2);
public:
  Overlap_Slater_Right(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, int ql, int qr);
};

class Overlap_BCS: public Overlap {
protected:
  Matrix work;
  Matrix m_ca, m_ac, m_aa, m_sa;
  int ntotal, nsites;
  int lcore, lactive, lactive_size;
  int rcore, ractive, ractive_size;
  int parity;
  bool this_site_up;
public:
  Overlap_BCS() {}
  dtype operator() (const vector<bool>& c1, const vector<bool>& c2);
};

class Overlap_BCS_Left: public Overlap_BCS {
protected:
  void build_cc_block(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr);
  void build_ac_block(const vector<bool>& c1, const vector<bool>& c2);
  void build_ca_block(const vector<bool>& c1, const vector<bool>& c2);
  void build_aa_block(const vector<bool>& c1, const vector<bool>& c2);
public:
  Overlap_BCS_Left(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, int ql, int qr);
};

class Overlap_BCS_Right: public Overlap_BCS {
protected:
  void build_cc_block(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr);
  void build_ac_block(const vector<bool>& c1, const vector<bool>& c2);
  void build_ca_block(const vector<bool>& c1, const vector<bool>& c2);
  void build_aa_block(const vector<bool>& c1, const vector<bool>& c2);
public:
  Overlap_BCS_Right(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, int ql, int qr);
};
#endif
