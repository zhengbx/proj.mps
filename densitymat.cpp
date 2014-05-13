#include "densitymat.h"

boost::shared_ptr<SchmidtBasis> SlaterDM::basis(int cut) const {
  Matrix la, lc, ra, rc; // active and core orbitals
  vector<d_real> lweight; // this is the square of natural orbital norm

  if (cut < 0 || cut > nsites) {
    cout << "Invalid cut position";
    abort();
  }

  if (cut <= nsites/2) {
    Matrix pRDM = m_coefs.topRows(cut).transpose() * m_coefs.topRows(cut);
    Eigen::SelfAdjointEigenSolver<Matrix> es(pRDM);
    auto D = es.eigenvalues();
    Matrix new_orbs = m_coefs * es.eigenvectors();

    vector<int> r, l, a;
    for (int i = 0; i < norbs; ++i) {
      if (D(i) < params.thr1p) {
        r.push_back(i);
      } else if (D(i) > 1. - params.thr1p) {
        l.push_back(i);
      } else {
        a.push_back(i);
        lweight.push_back(D(i));
      }
    }

    rc.resize(nsites-cut, r.size());
    lc.resize(cut, l.size());
    ra.resize(nsites-cut, a.size());
    la.resize(cut, a.size());

    for (int i = 0; i < r.size(); ++i) {
      rc.col(i) = new_orbs.col(r[i]).tail(nsites-cut) / sqrt(1.-D(r[i]));
    }
    for (int i = 0; i < l.size(); ++i) {
      lc.col(i) = new_orbs.col(l[i]).head(cut) / sqrt(D(l[i]));
    }
    for (int i = 0; i < a.size(); ++i) {
      la.col(i) = new_orbs.col(a[i]).head(cut) / sqrt(D(a[i]));
      ra.col(i) = new_orbs.col(a[i]).tail(nsites-cut) / sqrt(1.-D(a[i]));
    }

  } else {
    Matrix pRDM = m_coefs.bottomRows(nsites-cut).transpose() * m_coefs.bottomRows(nsites-cut);
    Eigen::SelfAdjointEigenSolver<Matrix> es(pRDM);    
    auto D = es.eigenvalues();
    Matrix new_orbs = m_coefs * es.eigenvectors();

    vector<int> r, l, a;
    for (int i = 0; i < norbs; ++i) {
      if (D(i) < params.thr1p) {
        l.push_back(i);
      } else if (D(i) > 1. - params.thr1p) {
        r.push_back(i);
      } else {
        a.push_back(i);
        lweight.push_back(1.-D(i));
      }
    }

    rc.resize(nsites-cut, r.size());
    lc.resize(cut, l.size());
    ra.resize(nsites-cut, a.size());
    la.resize(cut, a.size());

    for (int i = 0; i < l.size(); ++i) {
      lc.col(i) = new_orbs.col(l[i]).head(cut) / sqrt(1.-D(l[i]));
    }
    for (int i = 0; i < r.size(); ++i) {
      rc.col(i) = new_orbs.col(r[i]).tail(nsites-cut) / sqrt(D(r[i]));
    }
    for (int i = 0; i < a.size(); ++i) {
      la.col(i) = new_orbs.col(a[i]).head(cut) / sqrt(1.-D(a[i]));
      ra.col(i) = new_orbs.col(a[i]).tail(nsites-cut) / sqrt(D(a[i]));
    }

  }
  return boost::shared_ptr<SchmidtBasis>(new SchmidtBasis_Slater(lc, rc, la, ra, lweight));
}

boost::shared_ptr<SchmidtBasis> BCSDM::basis(int cut) const {
  Matrix la, lc, ra, rc; // left and right orbitals, core and active
  vector<d_real> lweight; // this is the square of natural orbital norm

  if (cut < 0 || cut > nsites) {
    cout << "Invalid cut position";
    abort();
  }

  if (cut <= nsites/2) {
    Matrix pRDM = m_coefs.topRows(cut).transpose() * m_coefs.topRows(cut) + 
        m_coefs.bottomRows(nsites).topRows(cut).transpose() * m_coefs.bottomRows(nsites).topRows(cut);
    Eigen::SelfAdjointEigenSolver<Matrix> es(pRDM);
    auto D = es.eigenvalues();
    Matrix new_orbs = m_coefs * es.eigenvectors();

    vector<int> r, l, a;
    for (int i = 0; i < nsites; ++i) {
      if (D(i) < params.thr1p) {
        r.push_back(i);
      } else if (D(i) > 1. - params.thr1p) {
        l.push_back(i);
      } else {
        a.push_back(i);
        lweight.push_back(D(i));
      }
    }

    rc.resize((nsites-cut)*2, r.size());
    lc.resize(cut*2, l.size());
    ra.resize((nsites-cut)*2, a.size());
    la.resize(cut*2, a.size());

    for (int i = 0; i < r.size(); ++i) {
      rc.col(i).head(nsites-cut) = new_orbs.col(r[i]).segment(cut, nsites-cut) / sqrt(1.-D(r[i]));
      rc.col(i).tail(nsites-cut) = new_orbs.col(r[i]).tail(nsites-cut) / sqrt(1.-D(r[i]));
    }
    for (int i = 0; i < l.size(); ++i) {
      lc.col(i).head(cut) = new_orbs.col(l[i]).head(cut) / sqrt(D(l[i]));
      lc.col(i).tail(cut) = new_orbs.col(l[i]).segment(nsites, cut) / sqrt(D(l[i]));
    }
    for (int i = 0; i < a.size(); ++i) {
      la.col(i).head(cut) = new_orbs.col(a[i]).head(cut) / sqrt(D(a[i]));
      la.col(i).tail(cut) = new_orbs.col(a[i]).segment(nsites, cut) / sqrt(D(a[i]));
      ra.col(i).head(nsites-cut) = new_orbs.col(a[i]).segment(cut, nsites-cut) / sqrt(1.-D(a[i]));
      ra.col(i).tail(nsites-cut) = new_orbs.col(a[i]).tail(nsites-cut) / sqrt(1.-D(a[i]));
    }
    
  } else {
    Matrix pRDM = m_coefs.topRows(nsites).bottomRows(nsites-cut).transpose() * m_coefs.topRows(nsites).bottomRows(nsites-cut) + 
        m_coefs.bottomRows(nsites-cut).transpose() * m_coefs.bottomRows(nsites-cut);
    Eigen::SelfAdjointEigenSolver<Matrix> es(pRDM);
    auto D = es.eigenvalues();
    Matrix new_orbs = m_coefs * es.eigenvectors();

    vector<int> r, l, a;
    for (int i = 0; i < nsites; ++i) {
      if (D(i) < params.thr1p) {
        l.push_back(i);
      } else if (D(i) > 1. - params.thr1p) {
        r.push_back(i);
      } else {
        a.push_back(i);
        lweight.push_back(1.-D(i));
      }
    }

    rc.resize((nsites-cut)*2, r.size());
    lc.resize(cut*2, l.size());
    ra.resize((nsites-cut)*2, a.size());
    la.resize(cut*2, a.size());

    for (int i = 0; i < l.size(); ++i) {
      lc.col(i).head(cut) = new_orbs.col(l[i]).head(cut) / sqrt(1.-D(l[i]));
      lc.col(i).tail(cut) = new_orbs.col(l[i]).segment(nsites, cut) / sqrt(1.-D(l[i]));
    }
    for (int i = 0; i < r.size(); ++i) {
      rc.col(i).head(nsites-cut) = new_orbs.col(r[i]).segment(cut, nsites-cut) / sqrt(D(r[i]));
      rc.col(i).tail(nsites-cut) = new_orbs.col(r[i]).tail(nsites-cut) / sqrt(D(r[i]));
    }
    for (int i = 0; i < a.size(); ++i) {
      la.col(i).head(cut) = new_orbs.col(a[i]).head(cut) / sqrt(1.-D(a[i]));
      la.col(i).tail(cut) = new_orbs.col(a[i]).segment(nsites, cut) / sqrt(1.-D(a[i]));
      ra.col(i).head(nsites-cut) = new_orbs.col(a[i]).segment(cut, nsites-cut) / sqrt(D(a[i]));
      ra.col(i).tail(nsites-cut) = new_orbs.col(a[i]).tail(nsites-cut) / sqrt(D(a[i]));
    }
  }

  return boost::shared_ptr<SchmidtBasis>(new SchmidtBasis_BCS(lc, rc, la, ra, lweight));
}

