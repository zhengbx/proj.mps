#include "mps_gen.h"

QSDArray<3, Quantum> generate_mps(boost::shared_ptr<SchmidtBasis> s1, boost::shared_ptr<SchmidtBasis> s2, bool additional) {
  QSDArray<3, Quantum> A;
  
  // assign quantums
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

  if (s1 -> nlcore() > s1 -> nrcore()) {
    cout << " Using Right Block" << endl;
  } else {
    cout << " Using  Left Block" << endl;
  }
  if (additional) {}
  return std::move(A);
}

