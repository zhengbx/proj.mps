#include "mps_gen.h"
#include <omp.h>

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
  vector<IVector<3>> blocks;
  for (int i = 0; i < ql.size(); ++i) {
    for (int j = 0; j < qp.size(); ++j) {
      for (int k = 0; k < qr.size(); ++k) {
        if (ql[i] + qp[j] == qr[k]) {
          blocks.push_back(make_array(i, j, k));
          IVector<3> idx = {i, j, k};
          A.reserve(idx);
          nelements += A.find(idx) -> second -> size();
        }
      }
    }
  }

  cout << "Site " << s1 -> lsites() <<  " Number of Elements " << nelements << "  Total memory " << (double)nelements / 1024 / 1024 / 1024 * 8 << " GB" << endl;

  #pragma omp parallel for schedule(dynamic, 1) default(shared)
  for (int i = 0; i < blocks.size(); ++i) {
    auto idx = blocks[i];
    DArray<3> dense;
    dense.reference(*(A.find(idx) -> second));
    compute_dense(dense, ql[idx[0]], idx[1], qr[idx[2]], s1, s2, use_left);
  }

  if (additional) {}
  return A;
}

void compute_dense(DArray<3>& d, int ql, int idx_p, int qr, boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, bool use_left) {
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

Overlap_Slater_Left::Overlap_Slater_Left(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, int ql, int qr) {
  int nsites = sr -> lsites();
  total_a = (nsites+qr)/2; // total a electron
  total_b = (nsites-qr)/2; // total b electron
  lactive_size = sl -> nactive();
  ractive_size = sr -> nactive();
  lcore = sl -> nlcore();  // electrons in left_basis core
  rcore = sr -> nlcore();  // electrons in right_basis core
  lactive_a = (nsites-1+ql)/2 - lcore; // a electrons in left_basis active space
  lactive_b = (nsites-1-ql)/2 - lcore; // b electrons in left_basis active space
  ractive_a = total_a - rcore; // a electrons in right_basis active space
  ractive_b = total_b - rcore; // a electrons in right_basis active space
  Matrix work_a(total_a, total_a);
  Matrix work_b(total_b, total_b);

  if (total_a - lactive_a > total_b - lactive_b) { // the state on this site is up
    this_site_up = true;
    parity = ((total_a + total_b) % 2 == 1) ? 1 : -1;
    work_a.SubMatrix(1, 1, 1, rcore) = sr -> lcore().Row(nsites);
    work_a.SubMatrix(2, lcore+1, 1, rcore) = sl -> lcore().t() * sr -> lcore().Rows(1, nsites-1);
    work_b.SubMatrix(1, lcore, 1, rcore) = sl -> lcore().t() * sr -> lcore().Rows(1, nsites-1);
  } else {
    this_site_up = false;
    parity = (total_b % 2 == 1) ? 1 : -1;
    work_b.SubMatrix(1, 1, 1, rcore) = sr -> lcore().Row(nsites);
    work_b.SubMatrix(2, lcore+1, 1, rcore) = sl -> lcore().t() * sr -> lcore().Rows(1, nsites-1);    
    work_a.SubMatrix(1, lcore, 1, rcore) = sl -> lcore().t() * sr -> lcore().Rows(1, nsites-1);
  }
  m_ca = sl -> lcore().t() * sr -> lactive().Rows(1, nsites-1);
  m_ac = sl -> lactive().t() * sr -> lcore().Rows(1, nsites-1);
  m_aa = sl -> lactive().t() * sr -> lactive().Rows(1, nsites-1);
  m_sa = sr -> lactive().Row(nsites);
}

double Overlap_Slater_Left::operator() (const vector<bool>& c1, const vector<bool>& c2) {
  if (lactive_a) { // ac block in alpha spin
    int shift = this_site_up ? lcore+2:lcore+1, count = 0;
    for (int i = 0; i < lactive_size; ++i) {
      if (c1[i]) {
        work_a.SubMatrix(shift+count, shift+count, 1, rcore) = m_ac.Row(i+1);
        ++count;        
      }
    }
  }

  if (lactive_b) { // ac block in beta spin
    int shift = this_site_up ? lcore+1:lcore+2, count = 0;
    for (int i = 0; i < lactive_size; ++i) { // for beta spin
      if (c1[i+lactive_size]) {
        work_b.SubMatrix(shift+count, shift+count, 1, rcore) = m_ac.Row(i+1);
        ++count;        
      }
    }
  }

  if (ractive_a) { // ca, sa blocks in alpha spin
    int count = 0;
    for (int i = 0; i < ractive_size; ++i) {
      if (c2[i]) {
        if (this_site_up) {
          work_a(1, rcore+count+1) = m_sa(1, i+1);
          work_a.SubMatrix(2, lcore+1, rcore+count+1, rcore+count+1) = m_ca.Column(i+1);          
        } else {
          work_a.SubMatrix(1, lcore, rcore+count+1, rcore+count+1) = m_ca.Column(i+1); 
        }
        ++count;
      }
    }
  }

  if (ractive_b) { // ca, sa blocks in beta spin
    int count = 0;
    for (int i = 0; i < ractive_size; ++i) {
      if (c2[i+ractive_size]) {
        if (this_site_up) {
          work_b.SubMatrix(1, lcore, rcore+count+1, rcore+count+1) = m_ca.Column(i+1);
        } else {
          work_b(1, rcore+count+1) = m_sa(1, i+1);
          work_b.SubMatrix(2, lcore+1, rcore+count+1, rcore+count+1) = m_ca.Column(i+1);
        }
        ++count;
      }
    }
  }

  if (lactive_a && ractive_a) { // aa block in alpha spin
    int shift = this_site_up ? lcore+2:lcore+1, count_l = 0;
    for (int i = 0; i < lactive_size; ++i) {
      if (c1[i]) {
        int count_r = 0;
        for (int j = 0; j < ractive_size; ++j) {
          if (c2[j]) {
            work_a(shift+count_l, rcore+count_r+1) = m_aa(i+1, j+1);
            ++count_r;
          }
        }
        ++count_l;
      }
    }
  }

  if (lactive_b && ractive_b) {
    int shift = this_site_up ? lcore+1:lcore+2, count_l = 0;
    for (int i = 0; i < lactive_size; ++i) {
      if (c1[i+lactive_size]) {
        int count_r = 0;
        for (int j = 0; i < ractive_size; ++j) {
          if (c2[j+ractive_size]) {
            work_b(shift+count_l, rcore+count_r+1) = m_aa(i+1, j+1);
            ++count_r;
          }
        }
        ++count_l;
      }
    }
  }

  double detA = (total_a == 0) ? 1. : work_a.Determinant();
  double detB = (total_b == 0) ? 1. : work_b.Determinant();
  return parity * detA * detB;
}

Overlap_Slater_Right::Overlap_Slater_Right(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, int ql, int qr) {
  int nsites = sl -> rsites();
  total_a = (nsites+ql)/2; // total a electron
  total_b = (nsites-ql)/2; // total b electron
  lactive_size = sl -> nactive();
  ractive_size = sr -> nactive();
  lcore = sl -> nrcore();  // electrons in left_basis core
  rcore = sr -> nrcore();  // electrons in right_basis core
  lactive_a = total_a - lcore; // a electrons in left_basis active space
  lactive_b = total_b - lcore; // b electrons in left_basis active space
  ractive_a = (nsites-1+qr)/2 - rcore; // a electrons in right_basis active space
  ractive_b = (nsites-1-qr)/2 - rcore; // a electrons in right_basis active space
  Matrix work_a(total_a, total_a);
  Matrix work_b(total_b, total_b);

  if (total_a - lactive_a > total_b - lactive_b) { // the state on this site is up
    this_site_up = true;
    parity = 1;
    work_a.SubMatrix(1, 1, 1, lcore) = sl -> rcore().Row(1);
    work_a.SubMatrix(2, rcore+1, 1, lcore) = sr -> rcore().t() * sl -> rcore().Rows(2, nsites);
    work_b.SubMatrix(1, rcore, 1, lcore) = sr -> rcore().t() * sl -> rcore().Rows(2, nsites);
  } else {
    this_site_up = false;
    parity = (total_a % 2 == 0) ? 1 : -1;
    work_b.SubMatrix(1, 1, 1, lcore) = sl -> rcore().Row(1);
    work_b.SubMatrix(2, rcore+1, 1, lcore) = sr -> rcore().t() * sl -> rcore().Rows(2, nsites);
    work_a.SubMatrix(1, rcore, 1, lcore) = sr -> rcore().t() * sl -> rcore().Rows(2, nsites);
  }
  m_ca = sr -> rcore().t() * sl -> ractive().Rows(2, nsites);
  m_ac = sr -> ractive().t() * sl -> rcore().Rows(2, nsites);
  m_aa = sr -> ractive().t() * sl -> ractive().Rows(2, nsites);
  m_sa = sl -> ractive().Row(1);
}

double Overlap_Slater_Right::operator() (const vector<bool>& c1, const vector<bool>& c2) {
  if (ractive_a) { // ac block of spin alpha
    int shift = this_site_up ? rcore+2:rcore+1, count = 0;
    for (int i = 0; i < ractive_size; ++i) {
      if (c2[i]) {
        work_a.SubMatrix(shift+count, shift+count, 1, lcore) = m_ac.Row(i+1);
        ++count;
      }
    }
  }

  if (ractive_b) { // ac block of spin beta
    int shift = this_site_up ? rcore+1:rcore+2, count = 0;
    for (int i = 0; i < ractive_size; ++i) {
      if (c2[i+ractive_size]) {
        work_b.SubMatrix(shift+count, shift+count, 1, lcore) = m_ac.Row(i+1);
        ++count;
      }
    }
  }

  if (lactive_a) { // ca, sa blocks of spin alpha
    int count = 0;
    for (int i = 0; i < lactive_size; ++i) {
      if (c1[i]) {
        if (this_site_up) {
          work_a(1, lcore+count+1) = m_sa(1, i+1);
          work_a.SubMatrix(2, rcore+1, lcore+count+1, lcore+count+1) = m_ca.Column(i+1);
        } else {
          work_a.SubMatrix(1, rcore, lcore+count+1, lcore+count+1) = m_ca.Column(i+1);
        }
        ++count;
      }
    }
  }

  if (lactive_b) {
    int count = 0;
    for (int i = 0; i < lactive_size; ++i) {
      if (c1[i+lactive_size]) {
        if (this_site_up) {
          work_b.SubMatrix(1, rcore, lcore+count+1, lcore+count+1) = m_ca.Column(i+1);
        } else {
          work_b(1, lcore+count+1) = m_sa(1, i+1);
          work_b.SubMatrix(2, rcore+1, lcore+count+1, lcore+count+1) = m_ca.Column(i+1);
        }
        ++count;
      }
    }
  }

  if (lactive_a && ractive_a) { // aa block in alpha spin
    int shift = this_site_up ? rcore+2 : rcore+1, count_r = 0;
    for (int i = 0; i < ractive_size; ++i) {
      if (c2[i]) {
        int count_l = 0;
        for (int j = 0; j < lactive_size; ++j) {
          if (c1[i]) {
            work_a(shift+count_r, lcore+count_l+1) = m_aa(i+1, j+1);
            ++count_l;
          }
        }
        ++count_r;
      }
    }
  }

  if (lactive_b && ractive_b) {
    int shift = this_site_up ? rcore+1 : rcore+2, count_r = 0;
    for (int i = 0; i < ractive_size; ++i) {
      if (c2[i+ractive_size]) {
        int count_l = 0;
        for (int j = 0; j < lactive_size; ++j) {
          if (c1[i+lactive_size]) {
            work_b(shift+count_r, lcore+count_l+1) = m_aa(i+1, j+1);
            ++count_l;
          }
        }
        ++count_r;
      }
    }
  }

  double detA = (total_a == 0) ? 1. : work_a.Determinant();
  double detB = (total_b == 0) ? 1. : work_b.Determinant();

  return parity * detA * detB;
}
