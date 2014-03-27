#include "mps_gen.h"
#include <omp.h>
#include <cassert>

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

  printf("Site %3d has %16d elements in %3d blocks   Total Memory = %12.6f GB\n", s1 -> lsites(), nelements, blocks.size(), double(nelements) / 1024 / 1024 / 1024 * 8);

  #pragma omp parallel for schedule(dynamic, 1) default(shared)
  for (int i = 0; i < blocks.size(); ++i) {
    auto idx = blocks[i];
    DArray<3> dense;
    dense.reference(*(A.find(idx) -> second));
    compute_dense(dense, ql[idx[0]], idx[1], qr[idx[2]], s1, s2, use_left, additional);
  }
  return A;
}

void compute_dense(DArray<3>& d, int ql, int idx_p, int qr, boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, bool use_left, bool additional) {
  // build overlap matrix
  boost::shared_ptr<Overlap> overlap_calculator;
  if (params.bcs && use_left) {
    overlap_calculator = boost::shared_ptr<Overlap>(new Overlap_BCS_Left(sl, sr, ql, qr));
  } else if (params.bcs && !use_left) {
    overlap_calculator = boost::shared_ptr<Overlap>(new Overlap_BCS_Right(sl, sr, ql, qr));
  } else if (!params.bcs && use_left) {
    overlap_calculator = boost::shared_ptr<Overlap>(new Overlap_Slater_Left(sl, sr, ql, qr));
  } else {
    overlap_calculator = boost::shared_ptr<Overlap>(new Overlap_Slater_Right(sl, sr, ql, qr));
  }
  
  auto it_l = sl -> iterator(ql); // ptr to active space iterator
  auto it_r = sr -> iterator(qr);
 
  if (additional) {
    int com_parity = common_parity(sl, it_l);
    for (int i = 0; i < d.shape(0); ++i) {
      vector<bool> c1 = it_l -> get_config(i, use_left);
      double factor = it_l -> get_schmidt_coef(i) * individual_parity(c1) * com_parity;
      for (int j = 0; j < d.shape(2); ++j) {
        vector<bool> c2 = it_r -> get_config(j, use_left);
        d(i, 0, j) = (*overlap_calculator)(c1, c2) * factor;
      }
    }
  } else {
    for (int i = 0; i < d.shape(0); ++i) {
      vector<bool> c1 = it_l -> get_config(i, use_left);
      for (int j = 0; j < d.shape(2); ++j) {
        vector<bool> c2 = it_r -> get_config(j, use_left);
        d(i, 0, j) = (*overlap_calculator)(c1, c2);
      }
    }
  }
}

/*
 * compute parity
 *  Generic order of BCS wavefunction is defined as
 *  |core_l>|core_r>|active>
 *  the order of schmidt rep is
 *  |core_l>|active_l>|core_r>|active_r>
 *  so the common factor is nla_occ * nrc
 *  individual factor is the number of inverse of {|al>|ar>}
 * 
 *  Generic order of Slater determinant is
 *  |core_la>|cora_ra>|active_a>|core_lb>|cora_rb>|active_b>
 *  schmidt rep order is
 *  |core_la>|active_la>|core_lb>|active_lb>|cora_ra>|active_ra>|cora_rb>|active_rb>
 *  the common factor is thus nrc*(nla_a+nlc+nla_b)+nra_a*(nlc+nla_b)+nrc*nla_b
 *  = nrc*(nla_a+nlc) + nra_a*(nlc+nla_b)
 *  individual factor is then the number od inverse {|ala>|ara>} * {|alb>|arb>}
 */

int common_parity(boost::shared_ptr<SchmidtBasis> s, boost::shared_ptr<ActiveSpaceIterator> it) {
  int noi;
  if (params.bcs) {
    noi = it->occs()[0] * s->nrcore();
  } else {
    noi = s->nrcore() * (it->occs()[0] + s->nlcore()) + (s->nactive()-it->occs()[0]) * (s->nlcore() + it->occs()[1]);
  }
  return (noi % 2 == 0) ? 1 : -1;  
}

int individual_parity(const vector<bool>& bits) {
  int noi;
  if (params.bcs) {
    vector<int> l, r;
    for (int i = 0; i < bits.size(); ++i) {
      if (bits[i]) {
        l.push_back(i);
      } else {
        r.push_back(i);
      }
    }
    int pos = 0;
    for (int i = 0; i < r.size(); ++i) {
      while (pos < l.size() && r[i] > l[pos]) ++pos;
      noi += l.size() - pos;
    }
  } else {
    vector<int> la, ra, lb, rb;
    for (int i = 0; i < bits.size()/2; ++i) {
      if (bits[i]) {
        la.push_back(i);
      } else {
        ra.push_back(i);
      }
      if (bits[i+bits.size()/2]) {
        lb.push_back(i);
      } else {
        rb.push_back(i);
      }
    }
    int pos = 0;
    for (int i = 0; i < ra.size(); ++i) {
      while (pos < la.size() && ra[i] > la[pos]) ++pos;
      noi += la.size() - pos;
    }
    pos = 0;
    for (int i = 0; i < rb.size(); ++i) {
      while (pos < lb.size() && rb[i] > lb[pos]) ++pos;
      noi += lb.size() - pos;
    }
  }
  return (noi % 2 == 0) ? 1 : -1;
}
/* 
 * Overlap_Slater
*/

Overlap_Slater_Left::Overlap_Slater_Left(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, int ql, int qr) {
  nsites = sr -> lsites();
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
  work_a.ReSize(total_a, total_a);
  work_b.ReSize(total_b, total_b);
  work_a = 0; work_b = 0;

  m_ca = sl -> lcore().t() * sr -> lactive().Rows(1, nsites-1);
  m_ac = sl -> lactive().t() * sr -> lcore().Rows(1, nsites-1);
  m_aa = sl -> lactive().t() * sr -> lactive().Rows(1, nsites-1);
  m_sa = sr -> lactive().Row(nsites);

  build_cc_block(sl, sr);
}

Overlap_Slater_Right::Overlap_Slater_Right(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, int ql, int qr) {
  nsites = sl -> rsites();
  total_a = (nsites-ql)/2; // total a electron
  total_b = (nsites+ql)/2; // total b electron
  lactive_size = sl -> nactive();
  ractive_size = sr -> nactive();
  lcore = sl -> nrcore();  // electrons in left_basis core
  rcore = sr -> nrcore();  // electrons in right_basis core
  lactive_a = total_a - lcore; // a electrons in left_basis active space
  lactive_b = total_b - lcore; // b electrons in left_basis active space
  ractive_a = (nsites-1-qr)/2 - rcore; // a electrons in right_basis active space
  ractive_b = (nsites-1+qr)/2 - rcore; // a electrons in right_basis active space
  work_a.ReSize(total_a, total_a);
  work_b.ReSize(total_b, total_b);
  work_a = 0; work_b = 0;

  m_ca = sr -> rcore().t() * sl -> ractive().Rows(2, nsites);
  m_ac = sr -> ractive().t() * sl -> rcore().Rows(2, nsites);
  m_aa = sr -> ractive().t() * sl -> ractive().Rows(2, nsites);
  m_sa = sl -> ractive().Row(1);

  build_cc_block(sl, sr);
}

void Overlap_Slater_Left::build_cc_block(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr) {
  assert(total_a == lactive_a + lcore + 1 || total_b == lactive_b + lcore + 1);
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
}

void Overlap_Slater_Right::build_cc_block(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr) {
  assert(total_a == ractive_a + rcore + 1 || total_b == ractive_b + rcore + 1);
  if (total_a - ractive_a > total_b - ractive_b) { // the state on this site is up
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
}

double Overlap_Slater::operator() (const vector<bool>& c1, const vector<bool>& c2) {
  build_ac_block(c1, c2);
  build_ca_block(c1, c2);
  build_aa_block(c1, c2);

  double detA = (total_a == 0) ? 1. : work_a.Determinant();
  double detB = (total_b == 0) ? 1. : work_b.Determinant();
 
  return parity * detA * detB;
}

void Overlap_Slater_Left::build_ac_block(const vector<bool>& c1, const vector<bool>& c2) {
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
}

void Overlap_Slater_Right::build_ac_block(const vector<bool>& c1, const vector<bool>& c2) {
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

}

void Overlap_Slater_Left::build_ca_block(const vector<bool>& c1, const vector<bool>& c2) {
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
}

void Overlap_Slater_Right::build_ca_block(const vector<bool>& c1, const vector<bool>& c2) {
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
}

void Overlap_Slater_Left::build_aa_block(const vector<bool>& c1, const vector<bool>& c2) {
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
        for (int j = 0; j < ractive_size; ++j) {
          if (c2[j+ractive_size]) {
            work_b(shift+count_l, rcore+count_r+1) = m_aa(i+1, j+1);
            ++count_r;
          }
        }
        ++count_l;
      }
    }
  }
}

void Overlap_Slater_Right::build_aa_block(const vector<bool>& c1, const vector<bool>& c2) {
  if (lactive_a && ractive_a) { // aa block in alpha spin
    int shift = this_site_up ? rcore+2 : rcore+1, count_r = 0;
    for (int i = 0; i < ractive_size; ++i) {
      if (c2[i]) {
        int count_l = 0;
        for (int j = 0; j < lactive_size; ++j) {
          if (c1[j]) {
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
          if (c1[j+lactive_size]) {
            work_b(shift+count_r, lcore+count_l+1) = m_aa(i+1, j+1);
            ++count_l;
          }
        }
        ++count_r;
      }
    }
  }
}

/* 
 * Overlap_BCS
*/

Overlap_BCS_Left::Overlap_BCS_Left(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, int ql, int qr) {
  nsites = sr -> lsites(); // number of sites in the left basis
  ntotal = qr + nsites;  // total number of orbitals: each orbital raises spin by 1, from initial s=-nsites
  lactive_size = sl -> nactive();
  ractive_size = sr -> nactive();
  lcore = sl -> nlcore();
  rcore = sr -> nlcore();
  lactive = ql + (nsites-1) - lcore;
  ractive = ntotal - rcore;

  work.ReSize(ntotal, ntotal);
  work = 0;

  m_ca = sl -> lcore().t() * (sr -> lactive().Rows(1, nsites-1) & sr -> lactive().Rows(nsites+1, 2*nsites-1));
  m_ac = sl -> lactive().t() * (sr -> lcore().Rows(1, nsites-1) & sr -> lcore().Rows(nsites+1, 2*nsites-1));
  m_aa = sl -> lactive().t() * (sr -> lactive().Rows(1, nsites-1) & sr -> lactive().Rows(nsites+1, 2*nsites-1));
  m_sa = sr -> lactive().Row(2*nsites) & sr -> lactive().Row(nsites);
  build_cc_block(sl, sr);
}

Overlap_BCS_Right::Overlap_BCS_Right(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, int ql, int qr) {
  nsites = sl -> rsites(); // number of sites in the left basis
  ntotal = -ql + nsites; // total number of orbitals: each orbital raises spin by 1, from initial s=-nsites
  lactive_size = sl -> nactive();
  ractive_size = sr -> nactive();
  lcore = sl -> nrcore();
  rcore = sr -> nrcore();
  lactive = ntotal - lcore;
  ractive = -qr + (nsites-1) - rcore;

  work.ReSize(ntotal, ntotal);
  work = 0;

  m_ca = sr -> rcore().t() * (sl -> ractive().Rows(2, nsites) & sl -> ractive().Rows(nsites+2, 2*nsites));
  m_ac = sr -> ractive().t() * (sl -> rcore().Rows(2, nsites) & sl -> rcore().Rows(nsites+2, 2*nsites));
  m_aa = sr -> ractive().t() * (sl -> ractive().Rows(2, nsites) & sl -> ractive().Rows(nsites+2, 2*nsites));
  m_sa = sl -> ractive().Row(nsites+1) & sl -> ractive().Row(1);
  
  build_cc_block(sl, sr);
}

void Overlap_BCS_Left::build_cc_block(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr) {
  assert(ntotal == lcore + lactive || ntotal == lcore + lactive + 2);
  if (ntotal > lcore + lactive) { // then this site is spin up
    this_site_up = true;
    parity = 1;
    work.SubMatrix(1, 2, 1, rcore) = sr -> lcore().Row(2*nsites) & sr -> lcore().Row(nsites);
    work.SubMatrix(3, lcore + 2, 1, rcore) = sl -> lcore().t() * (sr -> lcore().Rows(1, nsites-1) & sr -> lcore().Rows(nsites+1, 2*nsites-1));
  } else {
    this_site_up = false;
    parity = 1;
    work.SubMatrix(1, lcore, 1, rcore) = sl -> lcore().t() * (sr -> lcore().Rows(1, nsites-1) & sr -> lcore().Rows(nsites+1, 2*nsites-1));
  }
}

void Overlap_BCS_Right::build_cc_block(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr) {
  assert(ntotal == rcore + ractive || ntotal == rcore + ractive + 2);
  if (ntotal > rcore + ractive) { // then this site is spin up
    this_site_up = true;
    parity = 1;
    work.SubMatrix(1, 2, 1, lcore) = sl -> rcore().Row(nsites+1) & sl -> rcore().Row(1);
    work.SubMatrix(3, rcore + 2, 1, lcore) = sr -> rcore().t() * (sl -> rcore().Rows(2, nsites) & sl -> rcore().Rows(nsites+2, 2*nsites));
  } else {
    this_site_up = false;
    parity = 1;
    work.SubMatrix(1, rcore, 1, lcore) = sr -> rcore().t() * (sl -> rcore().Rows(2, nsites) & sl -> rcore().Rows(nsites+2, 2*nsites));
  }
}

double Overlap_BCS::operator() (const vector<bool>& c1, const vector<bool>& c2) {
  build_ac_block(c1, c2);
  build_ca_block(c1, c2);
  build_aa_block(c1, c2);
  
  double det = (ntotal == 0) ? 1. : work.Determinant();
  return parity * det;
}

void Overlap_BCS_Left::build_ac_block(const vector<bool>& c1, const vector<bool>& c2) {
  if (lactive) {
    int shift = this_site_up ? lcore+3 : lcore+1, count = 0;
    for (int i = 0; i < lactive_size; ++i) {
      if (c1[i]) {
        work.SubMatrix(shift+count, shift+count, 1, rcore) = m_ac.Row(i+1);
        ++count;
      }
    }
  }
}

void Overlap_BCS_Right::build_ac_block(const vector<bool>& c1, const vector<bool>& c2) {
  if (ractive) {
    int shift = this_site_up ? rcore+3 : rcore+1, count = 0;
    for (int i = 0; i < ractive_size; ++i) {
      if (c2[i]) {
        work.SubMatrix(shift+count, shift+count, 1, lcore) = m_ac.Row(i+1);
        ++count;
      }
    }
  }
}

void Overlap_BCS_Left::build_ca_block(const vector<bool>& c1, const vector<bool>& c2) {
  if (ractive) {
    int count = 0;
    for (int i = 0; i < ractive_size; ++i) {
      if (c2[i]) {
        if (this_site_up) {
          work.SubMatrix(1, lcore+2, rcore+count+1, rcore+count+1) = m_sa.Column(i+1) & m_ca.Column(i+1);
        } else {
          work.SubMatrix(1, lcore, rcore+count+1, rcore+count+1) = m_ca.Column(i+1);
        }
        ++count;
      }
    }
  }
}

void Overlap_BCS_Right::build_ca_block(const vector<bool>& c1, const vector<bool>& c2) {
  if (lactive) {
    int count = 0;
    for (int i = 0; i < lactive_size; ++i) {
      if (c1[i]) {
        if (this_site_up) {
          work.SubMatrix(1, rcore+2, lcore+count+1, lcore+count+1) = m_sa.Column(i+1) & m_ca.Column(i+1);
        } else {
          work.SubMatrix(1, rcore, lcore+count+1, lcore+count+1) = m_ca.Column(i+1);
        }
        ++count;
      }
    }
  }
}

void Overlap_BCS_Left::build_aa_block(const vector<bool>& c1, const vector<bool>& c2) {
  if (lactive && ractive) {
    int shift = this_site_up ? lcore+3 : lcore+1, count_l = 0;
    for (int i = 0; i < lactive_size; ++i) {
      if (c1[i]) {
        int count_r = 0;
        for (int j = 0; j < ractive_size; ++j) {
          if (c2[j]) {
            work(shift+count_l, rcore+count_r+1) = m_aa(i+1, j+1);
            ++count_r;
          }
        }
        ++count_l;
      }
    }
  }
}

void Overlap_BCS_Right::build_aa_block(const vector<bool>& c1, const vector<bool>& c2) {
  if (lactive && ractive) {
    int shift = this_site_up ? rcore+3 : rcore+1, count_r = 0;
    for (int i = 0; i < ractive_size; ++i) {
      if (c2[i]) {
        int count_l = 0;
        for (int j = 0; j < lactive_size; ++j) {
          if (c1[j]) {
            work(shift+count_r, lcore+count_l+1) = m_aa(i+1, j+1);
            ++count_l;
          }
        }
        ++count_r;
      }
    }
  }
}
