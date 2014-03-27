#include "densitymat.h"

boost::shared_ptr<SchmidtBasis> SlaterDM::basis(int cut) const {
  Matrix la(cut, 0), lc(cut, 0); // active orbitals
  Matrix ra(nsites-cut, 0), rc(nsites-cut, 0); // core orbitals
  vector<double> lweight; // this is the square of natural orbital norm
  SymmetricMatrix S;
  DiagonalMatrix D;
  Matrix rotmat, new_orbs;

  if (cut < 0 || cut > nsites) {
    cout << "Invalid cut position";
    abort();
  }

  if (cut <= nsites/2) {
    S << m_coefs.Rows(1, cut).t() * m_coefs.Rows(1, cut);
    Jacobi(S, D, rotmat);
    new_orbs = m_coefs * rotmat;
    for (int i = 0; i < norbs; ++i) {
      if (D(i+1, i+1) < params.thr1p) { // rc
        rc |= new_orbs.SubMatrix(cut+1, nsites, i+1, i+1) / sqrt(1.-D(i+1, i+1));
      } else if (D(i+1, i+1) > 1. - params.thr1p) { // lc
        lc |= new_orbs.SubMatrix(1, cut, i+1, i+1) / sqrt(D(i+1, i+1));
      } else { // la and ra
        lweight.push_back(D(i+1, i+1));
        la |= new_orbs.SubMatrix(1, cut, i+1, i+1) / sqrt(D(i+1, i+1));
        ra |= new_orbs.SubMatrix(cut+1, nsites, i+1, i+1) / sqrt(1.-D(i+1, i+1));
      }
    }
  } else {
    S << m_coefs.Rows(cut+1, nsites).t() * m_coefs.Rows(cut+1, nsites);
    Jacobi(S, D, rotmat);
    new_orbs = m_coefs * rotmat;
    for (int i = 0; i < norbs; ++i) {
      if (D(i+1, i+1) < params.thr1p) { // lc
        lc |= new_orbs.SubMatrix(1, cut, i+1, i+1) / sqrt(1.-D(i+1, i+1));
      } else if (D(i+1, i+1) > 1. - params.thr1p) { // rc
        rc |= new_orbs.SubMatrix(cut+1, nsites, i+1, i+1) / sqrt(D(i+1, i+1));
      } else { // la and ra
        lweight.push_back(1.-D(i+1, i+1));
        la |= new_orbs.SubMatrix(1, cut, i+1, i+1) / sqrt(1.-D(i+1, i+1));
        ra |= new_orbs.SubMatrix(cut+1, nsites, i+1, i+1) / sqrt(D(i+1, i+1));
      }
    }
  }
  return boost::shared_ptr<SchmidtBasis>(new SchmidtBasis_Slater(lc, rc, la, ra, lweight));
}

boost::shared_ptr<SchmidtBasis> BCSDM::basis(int cut) const {
  Matrix la(cut*2, 0), lc(cut*2, 0); // left orbitals
  Matrix ra((nsites-cut)*2, 0), rc((nsites-cut)*2, 0); // right orbitals
  vector<double> lweight; // this is the square of natural orbital norm
  SymmetricMatrix S;
  DiagonalMatrix D;
  Matrix rotmat, new_orbs;

  if (cut < 0 || cut > nsites) {
    cout << "Invalid cut position";
    abort();
  }

  if (cut <= nsites/2) {
    S << m_coefs.Rows(1, cut).t() * m_coefs.Rows(1, cut) + m_coefs.Rows(nsites+1, nsites+cut).t() * m_coefs.Rows(nsites+1, nsites+cut);
    Jacobi(S, D, rotmat);
    new_orbs = m_coefs * rotmat;
    for (int i = 0; i < nsites; ++i) {
      if (D(i+1, i+1) < params.thr1p) { // rc
        rc |= (new_orbs.SubMatrix(cut+1, nsites, i+1, i+1) & new_orbs.SubMatrix(nsites+cut+1, nsites*2, i+1, i+1)) / sqrt(1.-D(i+1, i+1));
      } else if (D(i+1, i+1) > 1. - params.thr1p) { // lc
        lc |= (new_orbs.SubMatrix(1, cut, i+1, i+1) & new_orbs.SubMatrix(nsites+1, nsites+cut, i+1, i+1)) / sqrt(D(i+1, i+1));
      } else { // la and ra
        lweight.push_back(D(i+1, i+1));
        la |= (new_orbs.SubMatrix(1, cut, i+1, i+1) & new_orbs.SubMatrix(nsites+1, nsites+cut, i+1, i+1)) / sqrt(D(i+1, i+1));
        ra |= (new_orbs.SubMatrix(cut+1, nsites, i+1, i+1) & new_orbs.SubMatrix(nsites+cut+1, nsites*2, i+1, i+1)) / sqrt(1.-D(i+1, i+1));
      }
    }
  } else {
    S << m_coefs.Rows(cut+1, nsites).t() * m_coefs.Rows(cut+1, nsites) + m_coefs.Rows(nsites+cut+1, nsites*2).t() * m_coefs.Rows(nsites+cut+1, nsites*2);
    Jacobi(S, D, rotmat);
    new_orbs = m_coefs * rotmat;
    for (int i = 0; i < nsites; ++i) {
      if (D(i+1, i+1) < params.thr1p) { // lc
        lc |= (new_orbs.SubMatrix(1, cut, i+1, i+1) & new_orbs.SubMatrix(nsites+1, nsites+cut, i+1, i+1)) / sqrt(1.-D(i+1, i+1));
      } else if (D(i+1, i+1) > 1. - params.thr1p) { // rc
        rc |= (new_orbs.SubMatrix(cut+1, nsites, i+1, i+1) & new_orbs.SubMatrix(nsites+cut+1, nsites*2, i+1, i+1)) / sqrt(D(i+1, i+1));
      } else { // la and ra
        lweight.push_back(1.-D(i+1, i+1));
        la |= (new_orbs.SubMatrix(1, cut, i+1, i+1) & new_orbs.SubMatrix(nsites+1, nsites+cut, i+1, i+1)) / sqrt(1.-D(i+1, i+1));
        ra |= (new_orbs.SubMatrix(cut+1, nsites, i+1, i+1) & new_orbs.SubMatrix(nsites+cut+1, nsites*2, i+1, i+1)) / sqrt(D(i+1, i+1));
      }
    }
  }

  return boost::shared_ptr<SchmidtBasis>(new SchmidtBasis_BCS(lc, rc, la, ra, lweight));
}

