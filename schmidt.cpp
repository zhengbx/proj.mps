#include "schmidt.h"
#include <algorithm>

std::ostream& operator <<(std::ostream& os, const SchmidtBasis& basis) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(10);
  os << "Left Sites " << basis.lsites() << "  Right Sites " << basis.rsites() << endl;
  os << "Core Orbitals  Left (" << basis.nlcore() << ")  Right (" << basis.nrcore() << ")" << endl;
  os << "Active Space (" << basis.nactive() << ") with Left Weights:\n";
  for (int i = 0; i < basis.nactive(); ++i) {
    os << basis.m_lweight[i] << "  ";
  }
  os << endl;
  os << "Quantum Numbers\t" << basis.quantums << endl;
  return os;
}

double SchmidtBasis::lweight_bound(int n) const {
  vector<double> temp(m_lweight);
  std::sort(temp.begin(), temp.end());
  double max = 1.;
  for (int i = 0; i < m_na; ++i) {
    max *= (i < m_na-n) ? (1.-temp[i]) : temp[i];
  }
  return max;
}

void SchmidtBasis_Slater::calc_q() {
  int nelec = m_lsites - 2 * m_nlc; // in active space
  for (int nelec_a = 0; nelec_a <= nelec; ++nelec_a) {
    if (lweight_bound(nelec_a) * lweight_bound(nelec-nelec_a) > params.thrnp) {
      quantums.push_back(nelec_a*2 - nelec);
    }
  }
}

void SchmidtBasis_BCS::calc_q() {
  for (int n = m_nlc; n <= m_nlc + m_na; ++n) {
    cout << n << " " << m_lsites - n << " " << lweight_bound(n - m_nlc) << endl;
    if (n % 2 == 0 && lweight_bound(n - m_nlc) > params.thrnp) {
      quantums.push_back(m_lsites - n);
    }
  }
}

void SchmidtBasis::calc_dim() {
  dims.resize(quantums.size(), 0);
  for (int i = 0; i < quantums.size(); ++i) {
    dims[i] = iterator(quantums[i]) -> size();
  }
}

std::shared_ptr<ActiveSpaceIterator> SchmidtBasis::iterator(int q) {
  auto it = m_it.find(q);
  if (it == m_it.end()) {
    auto ptr_asi = params.bcs ? 
      std::shared_ptr<ActiveSpaceIterator>(new ActiveSpaceIterator_Slater(q, this)):
      std::shared_ptr<ActiveSpaceIterator>(new ActiveSpaceIterator_BCS(q, this));
    m_it.insert(std::pair<int, std::shared_ptr<ActiveSpaceIterator>>(q, ptr_asi));
  }
  return m_it[q];
}
