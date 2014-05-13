#include "mps_gen.h"
#include <omp.h>
#include <boost/mpi.hpp>
#include <cassert>
#include "timer.h"

namespace mpi = boost::mpi;

QSTArray<dtype, 3, Quantum> generate_mps(boost::shared_ptr<SchmidtBasis> s1, boost::shared_ptr<SchmidtBasis> s2, bool additional) {
  Timer t;
  t.start();
  // assign quantums
  QSTArray<dtype, 3, Quantum> A;
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

  mpi::communicator world;
  
  int size = 8;
#ifdef _SINGLE
  size /= 2;
#endif

#ifdef _COMPLEX
  size *= 2;
#endif
  
  printf("On Processor %3d: Site %3d has %12d elements in %3d blocks    Total Memory = %12.6f GB\n", world.rank(), s1 -> lsites(), nelements, blocks.size(), d_real(nelements) / 1024 / 1024 / 1024 * size);
  if (!params.mem_test) {
    //#pragma omp parallel for schedule(dynamic, 1) default(shared)
    for (int i = 0; i < blocks.size(); ++i) {
      auto idx = blocks[i];
      TArray<dtype, 3> dense;
      dense.reference(*(A.find(idx) -> second));
      compute_dense(dense, ql[idx[0]], idx[1], qr[idx[2]], s1, s2, use_left, additional);
    }
  }
  t.pause();
  printf("On Processor %3d: Site %3d Time = %6.2f\n", world.rank(), s1 -> lsites(), t.time());
  
  return A;
}

void compute_dense(TArray<dtype, 3>& d, int ql, int idx_p, int qr, boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr, bool use_left, bool additional) {
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
      d_real factor = it_l -> get_schmidt_coef(i) * individual_parity(c1) * com_parity;
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
  int noi = 0;
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
  work_a.resize(total_a, total_a);
  work_b.resize(total_b, total_b);
  work_a.setZero();
  work_b.setZero();

  m_ca = sl -> lcore().transpose() * sr -> lactive().topRows(nsites-1);
  m_ac = sl -> lactive().transpose() * sr -> lcore().topRows(nsites-1);
  m_aa = sl -> lactive().transpose() * sr -> lactive().topRows(nsites-1);
  m_sa = sr -> lactive().row(nsites-1);

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
  work_a.resize(total_a, total_a);
  work_b.resize(total_b, total_b);
  work_a.setZero();
  work_b.setZero();

  m_ca = sr -> rcore().transpose() * sl -> ractive().bottomRows(nsites-1);
  m_ac = sr -> ractive().transpose() * sl -> rcore().bottomRows(nsites-1);
  m_aa = sr -> ractive().transpose() * sl -> ractive().bottomRows(nsites-1);
  m_sa = sl -> ractive().row(0);

  build_cc_block(sl, sr);
}

void Overlap_Slater_Left::build_cc_block(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr) {
  assert(total_a == lactive_a + lcore + 1 || total_b == lactive_b + lcore + 1);
  if (total_a - lactive_a > total_b - lactive_b) { // the state on this site is up
    this_site_up = true;
    parity = ((total_a + total_b) % 2 == 1) ? 1 : -1;
    work_a.topLeftCorner(1,rcore) = sr->lcore().row(nsites-1);
    work_a.block(1, 0, lcore, rcore) = sl->lcore().transpose() * sr->lcore().topRows(nsites-1);
    work_b.topLeftCorner(lcore, rcore) = sl->lcore().transpose() * sr->lcore().topRows(nsites-1);
  } else {
    this_site_up = false;
    parity = (total_b % 2 == 1) ? 1 : -1;
    work_b.topLeftCorner(1, rcore) = sr->lcore().row(nsites-1);
    work_b.block(1, 0, lcore, rcore) = sl->lcore().transpose() * sr->lcore().topRows(nsites-1);
    work_a.topLeftCorner(lcore, rcore) = sl->lcore().transpose() * sr->lcore().topRows(nsites-1);
  }
}

void Overlap_Slater_Right::build_cc_block(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr) {
  assert(total_a == ractive_a + rcore + 1 || total_b == ractive_b + rcore + 1);
  if (total_a - ractive_a > total_b - ractive_b) { // the state on this site is up
    this_site_up = true;
    parity = 1;
    work_a.topLeftCorner(1, lcore) = sl->rcore().row(0);
    work_a.block(1, 0, rcore, lcore) = sr->rcore().transpose() * sl->rcore().bottomRows(nsites-1);
    work_b.topLeftCorner(rcore, lcore) = sr -> rcore().transpose() * sl -> rcore().bottomRows(nsites-1);
  } else {
    this_site_up = false;
    parity = (total_a % 2 == 0) ? 1 : -1;
    work_b.topLeftCorner(1, lcore) = sl->rcore().row(0);
    work_b.block(1, 0, rcore, lcore) = sr->rcore().transpose() * sl->rcore().bottomRows(nsites-1);
    work_a.topLeftCorner(rcore, lcore) = sr->rcore().transpose() * sl->rcore().bottomRows(nsites-1);
  }
}

dtype Overlap_Slater::operator() (const vector<bool>& c1, const vector<bool>& c2) {
  build_ac_block(c1, c2);
  build_ca_block(c1, c2);
  build_aa_block(c1, c2);

  dtype detA = (total_a == 0) ? 1. : work_a.determinant();
  dtype detB = (total_b == 0) ? 1. : work_b.determinant();
 
  return (d_real)parity * detA * detB;
}

void Overlap_Slater_Left::build_ac_block(const vector<bool>& c1, const vector<bool>& c2) {
  if (lactive_a) { // ac block in alpha spin
    int shift = this_site_up ? lcore+1:lcore, count = 0;
    for (int i = 0; i < lactive_size; ++i) {
      if (c1[i]) {
        work_a.block(shift+count, 0, 1, rcore) = m_ac.row(i);
        ++count;
      }
    }
  }
  if (lactive_b) { // ac block in beta spin
    int shift = this_site_up ? lcore:lcore+1, count = 0;
    for (int i = 0; i < lactive_size; ++i) { // for beta spin
      if (c1[i+lactive_size]) {
        work_b.block(shift+count, 0, 1, rcore) = m_ac.row(i);
        ++count;        
      }
    }
  }
}

void Overlap_Slater_Right::build_ac_block(const vector<bool>& c1, const vector<bool>& c2) {
  if (ractive_a) { // ac block of spin alpha
    int shift = this_site_up ? rcore+1:rcore, count = 0;
    for (int i = 0; i < ractive_size; ++i) {
      if (c2[i]) {
        work_a.block(shift+count, 0, 1, lcore) = m_ac.row(i);
        ++count;
      }
    }
  }
  if (ractive_b) { // ac block of spin beta
    int shift = this_site_up ? rcore:rcore+1, count = 0;
    for (int i = 0; i < ractive_size; ++i) {
      if (c2[i+ractive_size]) {
        work_b.block(shift+count, 0, 1, lcore) = m_ac.row(i);
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
          work_a(0, rcore+count) = m_sa(0, i);
          work_a.col(rcore+count).segment(1, lcore) = m_ca.col(i);
        } else {
          work_a.col(rcore+count).head(lcore) = m_ca.col(i);
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
          work_b.col(rcore+count).head(lcore) = m_ca.col(i);
        } else {
          work_b(0, rcore+count) = m_sa(0, i);
          work_b.col(rcore+count).segment(1, lcore) = m_ca.col(i);
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
          work_a(0, lcore+count) = m_sa(0, i);
          work_a.col(lcore+count).segment(1, rcore) = m_ca.col(i);
        } else {
          work_a.col(lcore+count).head(rcore) = m_ca.col(i);
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
          work_b.col(lcore+count).head(rcore) = m_ca.col(i);
        } else {
          work_b(0, lcore+count) = m_sa(0, i);
          work_b.col(lcore+count).segment(1, rcore) = m_ca.col(i);
        }
        ++count;
      }
    }
  }
}

void Overlap_Slater_Left::build_aa_block(const vector<bool>& c1, const vector<bool>& c2) {
  if (lactive_a && ractive_a) { // aa block in alpha spin
    int shift = this_site_up ? lcore+1:lcore, count_l = 0;
    for (int i = 0; i < lactive_size; ++i) {
      if (c1[i]) {
        int count_r = 0;
        for (int j = 0; j < ractive_size; ++j) {
          if (c2[j]) {
            work_a(shift+count_l, rcore+count_r) = m_aa(i, j);
            ++count_r;
          }
        }
        ++count_l;
      }
    }
  }

  if (lactive_b && ractive_b) {
    int shift = this_site_up ? lcore:lcore+1, count_l = 0;
    for (int i = 0; i < lactive_size; ++i) {
      if (c1[i+lactive_size]) {
        int count_r = 0;
        for (int j = 0; j < ractive_size; ++j) {
          if (c2[j+ractive_size]) {
            work_b(shift+count_l, rcore+count_r) = m_aa(i, j);
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
    int shift = this_site_up ? rcore+1 : rcore, count_r = 0;
    for (int i = 0; i < ractive_size; ++i) {
      if (c2[i]) {
        int count_l = 0;
        for (int j = 0; j < lactive_size; ++j) {
          if (c1[j]) {
            work_a(shift+count_r, lcore+count_l) = m_aa(i, j);
            ++count_l;
          }
        }
        ++count_r;
      }
    }
  }
  if (lactive_b && ractive_b) {
    int shift = this_site_up ? rcore : rcore+1, count_r = 0;
    for (int i = 0; i < ractive_size; ++i) {
      if (c2[i+ractive_size]) {
        int count_l = 0;
        for (int j = 0; j < lactive_size; ++j) {
          if (c1[j+lactive_size]) {
            work_b(shift+count_r, lcore+count_l) = m_aa(i, j);
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

  work.resize(ntotal, ntotal);
  work.setZero();

  m_ca = sl->lcore().topRows(nsites-1).transpose() * sr->lactive().topRows(nsites-1) +
      sl->lcore().bottomRows(nsites-1).transpose() * sr->lactive().bottomRows(nsites).topRows(nsites-1);
  m_ac = sl->lactive().topRows(nsites-1).transpose() * sr->lcore().topRows(nsites-1) +
      sl->lactive().bottomRows(nsites-1).transpose() * sr->lcore().bottomRows(nsites).topRows(nsites-1);
  m_aa = sl->lactive().topRows(nsites-1).transpose() * sr->lactive().topRows(nsites-1) +
      sl->lactive().bottomRows(nsites-1).transpose() * sr->lactive().bottomRows(nsites).topRows(nsites-1);

  m_sa.resize(2, ractive_size);
  m_sa.row(0) = sr->lactive().row(2*nsites-1);
  m_sa.row(1) = sr->lactive().row(nsites-1);

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

  work.resize(ntotal, ntotal);
  work.setZero();

  m_ca = sr->rcore().topRows(nsites-1).transpose() * sl->ractive().topRows(nsites).bottomRows(nsites-1) + 
     sr->rcore().bottomRows(nsites-1).transpose() * sl->ractive().bottomRows(nsites-1);
  m_ac = sr->ractive().topRows(nsites-1).transpose() * sl->rcore().topRows(nsites).bottomRows(nsites-1) + 
     sr->ractive().bottomRows(nsites-1).transpose() * sl->rcore().bottomRows(nsites-1);
  m_aa = sr->ractive().topRows(nsites-1).transpose() * sl->ractive().topRows(nsites).bottomRows(nsites-1) + 
     sr->ractive().bottomRows(nsites-1).transpose() * sl->ractive().bottomRows(nsites-1);
  
  m_sa.resize(2, lactive_size);
  m_sa.row(0) = sl->ractive().row(nsites);
  m_sa.row(1) = sl->ractive().row(0);
  
  build_cc_block(sl, sr);
}

void Overlap_BCS_Left::build_cc_block(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr) {
  assert(ntotal == lcore + lactive || ntotal == lcore + lactive + 2);
  if (ntotal > lcore + lactive) { // then this site is spin up
    this_site_up = true;
    parity = 1;
    work.row(0).head(rcore) = sr->lcore().row(2*nsites-1);
    work.row(1).head(rcore) = sr->lcore().row(nsites-1);
    work.block(2, 0, lcore, rcore) = sl->lcore().topRows(nsites-1).transpose() * sr->lcore().topRows(nsites-1) +
        sl->lcore().bottomRows(nsites-1).transpose() * sr->lcore().bottomRows(nsites).topRows(nsites-1);
  } else {
    this_site_up = false;
    parity = 1;
    work.topLeftCorner(lcore, rcore) = sl->lcore().topRows(nsites-1).transpose() * sr->lcore().topRows(nsites-1) +
        sl->lcore().bottomRows(nsites-1).transpose() * sr->lcore().bottomRows(nsites).topRows(nsites-1);
  }
}

void Overlap_BCS_Right::build_cc_block(boost::shared_ptr<SchmidtBasis> sl, boost::shared_ptr<SchmidtBasis> sr) {
  assert(ntotal == rcore + ractive || ntotal == rcore + ractive + 2);
  if (ntotal > rcore + ractive) { // then this site is spin up
    this_site_up = true;
    parity = 1;
    work.row(0).head(lcore) = sl->rcore().row(nsites);
    work.row(1).head(lcore) = sl->rcore().row(0);
    work.block(2, 0, rcore, lcore) = sr->rcore().topRows(nsites-1).transpose() * sl->rcore().topRows(nsites).bottomRows(nsites-1) +
        sr->rcore().bottomRows(nsites-1).transpose() * sl->rcore().bottomRows(nsites-1);
  } else {
    this_site_up = false;
    parity = 1;
    work.topLeftCorner(rcore, lcore) = sr->rcore().topRows(nsites-1).transpose() * sl->rcore().topRows(nsites).bottomRows(nsites-1) +
        sr->rcore().bottomRows(nsites-1).transpose() * sl->rcore().bottomRows(nsites-1);
  }
}

dtype Overlap_BCS::operator() (const vector<bool>& c1, const vector<bool>& c2) {
  build_ac_block(c1, c2);
  build_ca_block(c1, c2);
  build_aa_block(c1, c2);
  
  dtype det = (ntotal == 0) ? 1. : work.determinant();
  return (d_real)parity * det;
}

void Overlap_BCS_Left::build_ac_block(const vector<bool>& c1, const vector<bool>& c2) {
  if (lactive) {
    int shift = this_site_up ? lcore+2 : lcore, count = 0;
    for (int i = 0; i < lactive_size; ++i) {
      if (c1[i]) {
        work.row(shift+count).head(rcore) = m_ac.row(i);
        ++count;
      }
    }
  }
}

void Overlap_BCS_Right::build_ac_block(const vector<bool>& c1, const vector<bool>& c2) {
  if (ractive) {
    int shift = this_site_up ? rcore+2 : rcore, count = 0;
    for (int i = 0; i < ractive_size; ++i) {
      if (c2[i]) {
        work.row(shift+count).head(lcore) = m_ac.row(i);
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
          work.col(rcore+count).head(2) = m_sa.col(i);
          work.col(rcore+count).segment(2, lcore) = m_ca.col(i);
        } else {
          work.col(rcore+count).head(lcore) = m_ca.col(i);
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
          work.col(lcore+count).head(2) = m_sa.col(i);
          work.col(lcore+count).segment(2, rcore) = m_ca.col(i);
        } else {
          work.col(lcore+count).head(rcore) = m_ca.col(i);
        }
        ++count;
      }
    }
  }
}

void Overlap_BCS_Left::build_aa_block(const vector<bool>& c1, const vector<bool>& c2) {
  if (lactive && ractive) {
    int shift = this_site_up ? lcore+2 : lcore, count_l = 0;
    for (int i = 0; i < lactive_size; ++i) {
      if (c1[i]) {
        int count_r = 0;
        for (int j = 0; j < ractive_size; ++j) {
          if (c2[j]) {
            work(shift+count_l, rcore+count_r) = m_aa(i, j);
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
    int shift = this_site_up ? rcore+2 : rcore, count_r = 0;
    for (int i = 0; i < ractive_size; ++i) {
      if (c2[i]) {
        int count_l = 0;
        for (int j = 0; j < lactive_size; ++j) {
          if (c1[j]) {
            work(shift+count_r, lcore+count_l) = m_aa(i, j);
            ++count_l;
          }
        }
        ++count_r;
      }
    }
  }
}
