#include "mps_gen.h"

QSDArray<3, Quantum> generate_mps(boost::shared_ptr<SchmidtBasis> s1, boost::shared_ptr<SchmidtBasis> s2, bool additional) {
  // assign quantums
  QSDArray<3, Quantum> A;
  TVector<Qshapes<Quantum>, 3> qshapes;
  TVector<Dshapes, 3> dshapes;
  vector<int> ql(s1->get_q()), qr(s2->get_q()), dl(s1->get_d()), dr(s2->get_d()), qp, dp;
  for (int i = 0; i < ql.size(); ++i) {
    qshapes[0].push_back(Quantum(ql[i]));
    dshapes[0].push_back(dl[i]);
  }
  for (int i = 0; i < qr.size(); ++i) {
    qshapes[2].push_back(-Quantum(qr[i]));
    dshapes[2].push_back(dr[i]);
  }
  qp = {-1, 1};
  dp = {1, 1};
  for (int i = 0; i < qp.size(); ++i) {
    qshapes[1].push_back(Quantum(qp[i]));
    dshapes[1].push_back(dp[i]);
  }
  
  A.resize(Quantum::zero(), qshapes, dshapes, false); // do not allocate
  bool use_left = s1 -> nlcore() < s1 -> nrcore();
  //cout << (use_left ? " Using  Left Block" : " Using Right Block") << endl;
  // now fill in the blocks
  size_t nelements = 0;
  for (int i = 0; i < ql.size(); ++i) {
    for (int j = 0; j < qp.size(); ++j) {
      for (int k = 0; k < qr.size(); ++k) {
        if (ql[i] + qp[j] == qr[k]) {
          IVector<3> idx = {i,j,k};
          A.reserve(idx);
          DArray<3> dense;
          dense.reference(*(A.find(idx) -> second));
          nelements += dense.size();          
          //compute_dense(dense, ql[i], j, qr[k], s1, s2, use_left);
        }
      }
    }
  }

  cout << "Site " << s1 -> lsites() <<  " Number of Elements " << nelements << "  Total memory " << (double)nelements / 1024 / 1024 / 1024 * 8 << " GB" << endl;

  if (additional) {}
  return A;
}

void compute_dense(DArray<3> d, int ql, int idx_p, int qr, boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, bool use_left) {
  // build overlap matrix
  
  auto it_l = sl -> iterator(ql); // ptr to active space iterator
  auto it_r = sr -> iterator(qr);
  
  for (int i = 0; i < d.shape(0); ++i) {
    vector<bool> c1 = it_l -> get_config(i, use_left);
    for (int j = 0; j < d.shape(2); ++j) {
      vector<bool> c2 = it_r -> get_config(j, use_left);
      //d(i, idx_p, j) = overlap(c1, c2);
    }
  }
}
