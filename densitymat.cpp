#include "densitymat.h"
#include <utility>

boost::shared_ptr<SchmidtBasis> SlaterDM::basis(int cut) const {
  Matrix la, lc, ra, rc; // active and core orbitals
  vector<d_real> lweight; // this is the square of natural orbital norm

  if (cut < 0 || cut > nsites) {
    cout << "Invalid cut position";
    abort();
  }

  if (cut <= nsites/2) {
    Matrix pOverlap = m_coefs.topRows(cut).adjoint() * m_coefs.topRows(cut);
    Eigen::SelfAdjointEigenSolver<Matrix> es(pOverlap);
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
    Matrix pOverlap = m_coefs.bottomRows(nsites-cut).adjoint() * m_coefs.bottomRows(nsites-cut);
    Eigen::SelfAdjointEigenSolver<Matrix> es(pOverlap);    
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
    Matrix pOverlap = m_coefs.topRows(cut).adjoint() * m_coefs.topRows(cut) + 
        m_coefs.bottomRows(nsites).topRows(cut).adjoint() * m_coefs.bottomRows(nsites).topRows(cut);
    Eigen::SelfAdjointEigenSolver<Matrix> es(pOverlap);
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
    Matrix pOverlap = m_coefs.topRows(nsites).bottomRows(nsites-cut).adjoint() * m_coefs.topRows(nsites).bottomRows(nsites-cut) + 
        m_coefs.bottomRows(nsites-cut).adjoint() * m_coefs.bottomRows(nsites-cut);
    Eigen::SelfAdjointEigenSolver<Matrix> es(pOverlap);
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

boost::shared_ptr<SchmidtBasis> SlaterDM::basis_k(int cut) const {
  Matrix la, lc, ra, rc;
  vector<d_real> lweight;
  vector<Matrix> new_orbs;
  vector<Eigen::Matrix<d_real, Eigen::Dynamic, 1>> D;

  vector<std::pair<int, int>> r, l, a;
  int kmax = params.kpoints[0].get_N();
  kpoint kl(kmax, 0), kr(kmax, 0);
  vector<kpoint> ka;

  if (cut < 0 || cut > nsites) {
    cout << "Invalid cut position";
    abort();
  }

  int offset = 0;  

  if (cut <= nsites/2) {
    for (int i = 0; i < params.kpoints.size(); ++i) {
      if (norb_k[i] == 0) {
        Eigen::Matrix<d_real, Eigen::Dynamic,1> D1;
        Matrix new_orbs1;
        D.push_back(D1);
        new_orbs.push_back(new_orbs1);
        continue;
      }
      Matrix pOverlap = m_coefs.block(0, offset, cut, norb_k[i]).adjoint() * m_coefs.block(0, offset, cut, norb_k[i]);
      Eigen::SelfAdjointEigenSolver<Matrix> es(pOverlap);
      D.push_back(es.eigenvalues());
      new_orbs.push_back(m_coefs.block(0, offset, nsites, norb_k[i]) * es.eigenvectors());

      for (int j = 0; j < norb_k[i]; ++j) {
        if (D.back()(j) < params.thr1p) {
          r.push_back(std::pair<int, int>(i,j));
          kr += params.kpoints[i];
        } else if (D.back()(j) > 1.-params.thr1p) {
          l.push_back(std::pair<int, int>(i,j));
          kl += params.kpoints[i];
        } else {
          a.push_back(std::pair<int, int>(i,j));
          ka.push_back(params.kpoints[i]);
          lweight.push_back(D.back()(j));
        }
      }
      offset += norb_k[i];
    }
    
    rc.resize(nsites-cut, r.size());
    lc.resize(cut, l.size());
    ra.resize(nsites-cut, a.size());
    la.resize(cut, a.size());

    for (int i = 0; i < r.size(); ++i) {
      rc.col(i) = new_orbs[r[i].first].col(r[i].second).tail(nsites-cut) / sqrt(1.-D[r[i].first](r[i].second));
    }
    for (int i = 0; i < l.size(); ++i) {
      lc.col(i) = new_orbs[l[i].first].col(l[i].second).head(cut) / sqrt(D[l[i].first](l[i].second));
    }
    for (int i = 0; i < a.size(); ++i) {
      la.col(i) = new_orbs[a[i].first].col(a[i].second).head(cut) / sqrt(D[a[i].first](a[i].second));
      ra.col(i) = new_orbs[a[i].first].col(a[i].second).tail(nsites-cut) / sqrt(1.-D[a[i].first](a[i].second));
    }

  } else {
    for (int i = 0; i < params.kpoints.size(); ++i) {
      if (norb_k[i] == 0) {
        Eigen::Matrix<d_real, Eigen::Dynamic,1> D1;
        Matrix new_orbs1;
        D.push_back(D1);
        new_orbs.push_back(new_orbs1);
        continue;
      }
      Matrix pOverlap = m_coefs.block(cut, offset, nsites-cut, norb_k[i]).adjoint() * m_coefs.block(cut, offset, nsites-cut, norb_k[i]);
      Eigen::SelfAdjointEigenSolver<Matrix> es(pOverlap);
      D.push_back(es.eigenvalues());
      new_orbs.push_back(m_coefs.block(0, offset, nsites, norb_k[i]) * es.eigenvectors());

      for (int j = 0; j < norb_k[i]; ++j) {
        if (D.back()(j) < params.thr1p) {
          l.push_back(std::pair<int, int>(i,j));
          kl += params.kpoints[i];
        } else if (D.back()(j) > 1.-params.thr1p) {
          r.push_back(std::pair<int, int>(i,j));
          kr += params.kpoints[i];          
        } else {
          a.push_back(std::pair<int, int>(i,j));
          ka.push_back(params.kpoints[i]);   
          lweight.push_back(1.-D.back()(j));
        }
      }
      offset += norb_k[i];
    }
    
    rc.resize(nsites-cut, r.size());
    lc.resize(cut, l.size());
    ra.resize(nsites-cut, a.size());
    la.resize(cut, a.size());

    for (int i = 0; i < l.size(); ++i) {
      lc.col(i) = new_orbs[l[i].first].col(l[i].second).head(cut) / sqrt(1.-D[l[i].first](l[i].second));
    }
    for (int i = 0; i < r.size(); ++i) {
      rc.col(i) = new_orbs[r[i].first].col(r[i].second).tail(nsites-cut) / sqrt(D[r[i].first](r[i].second));
    }
    for (int i = 0; i < a.size(); ++i) {
      la.col(i) = new_orbs[a[i].first].col(a[i].second).head(cut) / sqrt(1.-D[a[i].first](a[i].second));
      ra.col(i) = new_orbs[a[i].first].col(a[i].second).tail(nsites-cut) / sqrt(D[a[i].first](a[i].second));
    }
  }
  return boost::shared_ptr<SchmidtBasis>(new SchmidtBasis_Slater(lc, rc, la, ra, lweight, kl, kr, ka)); // temp
}

boost::shared_ptr<SchmidtBasis> BCSDM::basis_k(int cut) const {
  Matrix la, lc, ra, rc;
  vector<d_real> lweight;
  vector<Matrix> new_orbs;
  vector<Eigen::Matrix<d_real, Eigen::Dynamic, 1>> D;

  vector<std::pair<int, int>> r, l, a;
  int kmax = params.kpoints[0].get_N();
  kpoint kl(kmax, 0), kr(kmax, 0);
  vector<kpoint> ka;

  if (cut < 0 || cut > nsites) {
    cout << "Invalid cut position";
    abort();
  }

  int offset = 0;
  
  if (cut <= nsites/2) {
    for (int i = 0; i < params.kpoints.size(); ++i) {
      if (norb_k[i] == 0) {
        Eigen::Matrix<d_real, Eigen::Dynamic,1> D1;
        Matrix new_orbs1;
        D.push_back(D1);
        new_orbs.push_back(new_orbs1);
        continue;
      }
      Matrix pOverlap = m_coefs.block(0, offset, cut, norb_k[i]).adjoint() * m_coefs.block(0, offset, cut, norb_k[i]) +  m_coefs.block(nsites, offset, cut, norb_k[i]).adjoint() * m_coefs.block(nsites, offset, cut, norb_k[i]);
      Eigen::SelfAdjointEigenSolver<Matrix> es(pOverlap);
      D.push_back(es.eigenvalues());
      new_orbs.push_back(m_coefs.block(0, offset, nsites*2, norb_k[i]) * es.eigenvectors());

      for (int j = 0; j < norb_k[i]; ++j) {
        if (D.back()(j) < params.thr1p) {
          r.push_back(std::pair<int, int>(i,j));
          kr += params.kpoints[i];
        } else if (D.back()(j) > 1.-params.thr1p) {
          l.push_back(std::pair<int, int>(i,j));
          kl += params.kpoints[i];
        } else {
          a.push_back(std::pair<int, int>(i,j));
          ka.push_back(params.kpoints[i]);   
          lweight.push_back(D.back()(j));
        }
      }
      offset += norb_k[i];
    }

    rc.resize((nsites-cut)*2, r.size());
    lc.resize(cut*2, l.size());
    ra.resize((nsites-cut)*2, a.size());
    la.resize(cut*2, a.size());

    for (int i = 0; i < r.size(); ++i) {
      rc.col(i).head(nsites-cut) = new_orbs[r[i].first].col(r[i].second).segment(cut, nsites-cut) / sqrt(1.-D[r[i].first](r[i].second));
      rc.col(i).tail(nsites-cut) = new_orbs[r[i].first].col(r[i].second).tail(nsites-cut) / sqrt(1.-D[r[i].first](r[i].second));
    }
    for (int i = 0; i < l.size(); ++i) {
      lc.col(i).head(cut) = new_orbs[l[i].first].col(l[i].second).head(cut) / sqrt(D[l[i].first](l[i].second));
      lc.col(i).tail(cut) = new_orbs[l[i].first].col(l[i].second).segment(nsites, cut) / sqrt(D[l[i].first](l[i].second));
    }
    for (int i = 0; i < a.size(); ++i) {
      la.col(i).head(cut) = new_orbs[a[i].first].col(a[i].second).head(cut) / sqrt(D[a[i].first](a[i].second));
      la.col(i).tail(cut) = new_orbs[a[i].first].col(a[i].second).segment(nsites, cut) / sqrt(D[a[i].first](a[i].second));
      ra.col(i).head(nsites-cut) = new_orbs[a[i].first].col(a[i].second).segment(cut, nsites-cut) / sqrt(1.-D[a[i].first](a[i].second));
      ra.col(i).tail(nsites-cut) = new_orbs[a[i].first].col(a[i].second).tail(nsites-cut) / sqrt(1.-D[a[i].first](a[i].second));
    }

  } else {
    for (int i = 0; i < params.kpoints.size(); ++i) {
      if (norb_k[i] == 0) {
        Eigen::Matrix<d_real, Eigen::Dynamic,1> D1;
        Matrix new_orbs1;
        D.push_back(D1);
        new_orbs.push_back(new_orbs1);
        continue;
      }
      Matrix pOverlap = m_coefs.block(cut, offset, nsites-cut, norb_k[i]).adjoint() * m_coefs.block(cut, offset, nsites-cut, norb_k[i]) +  m_coefs.block(nsites+cut, offset, nsites-cut, norb_k[i]).adjoint() * m_coefs.block(nsites+cut, offset, nsites-cut, norb_k[i]);
      Eigen::SelfAdjointEigenSolver<Matrix> es(pOverlap);
      D.push_back(es.eigenvalues());
      new_orbs.push_back(m_coefs.block(0, offset, nsites*2, norb_k[i]) * es.eigenvectors());

      for (int j = 0; j < norb_k[i]; ++j) {
        if (D.back()(j) < params.thr1p) {
          l.push_back(std::pair<int, int>(i,j));
          kl += params.kpoints[i];
        } else if (D.back()(j) > 1.-params.thr1p) {
          r.push_back(std::pair<int, int>(i,j));
          kr += params.kpoints[i];
        } else {
          a.push_back(std::pair<int, int>(i,j));
          ka.push_back(params.kpoints[i]);   
          lweight.push_back(1.-D.back()(j));
        }
      }
      offset += norb_k[i];
    }

    rc.resize((nsites-cut)*2, r.size());
    lc.resize(cut*2, l.size());
    ra.resize((nsites-cut)*2, a.size());
    la.resize(cut*2, a.size());

    for (int i = 0; i < l.size(); ++i) {
      lc.col(i).head(cut) = new_orbs[l[i].first].col(l[i].second).head(cut) / sqrt(1.-D[l[i].first](l[i].second));
      lc.col(i).tail(cut) = new_orbs[l[i].first].col(l[i].second).segment(nsites, cut) / sqrt(1.-D[l[i].first](l[i].second));
    }
    for (int i = 0; i < r.size(); ++i) {
      rc.col(i).head(nsites-cut) = new_orbs[r[i].first].col(r[i].second).segment(cut, nsites-cut) / sqrt(D[r[i].first](r[i].second));
      rc.col(i).tail(nsites-cut) = new_orbs[r[i].first].col(r[i].second).tail(nsites-cut) / sqrt(D[r[i].first](r[i].second));
    }
    for (int i = 0; i < a.size(); ++i) {
      la.col(i).head(cut) = new_orbs[a[i].first].col(a[i].second).head(cut) / sqrt(1.-D[a[i].first](a[i].second));
      la.col(i).tail(cut) = new_orbs[a[i].first].col(a[i].second).segment(nsites, cut) / sqrt(1.-D[a[i].first](a[i].second));
      ra.col(i).head(nsites-cut) = new_orbs[a[i].first].col(a[i].second).segment(cut, nsites-cut) / sqrt(D[a[i].first](a[i].second));
      ra.col(i).tail(nsites-cut) = new_orbs[a[i].first].col(a[i].second).tail(nsites-cut) / sqrt(D[a[i].first](a[i].second));
    }
  }
  return boost::shared_ptr<SchmidtBasis>(new SchmidtBasis_BCS(lc, rc, la, ra, lweight, kl, kr, ka));
}
