#include "schmidt.h"
#include <algorithm>

vector<bool> complimentary(const vector<bool>& bits) {
  vector<bool> c(bits.size());
  for (int i = 0; i < bits.size(); ++i) {
    c[i] = (!bits[i]);
  }
  return c;
}

size_t choose(int iN, int iR){
  if (iR < 0 || iR > iN) {
      return 0;
  }
  if (iR > iN/2) {
    return choose(iN, iN-iR);
  }
  size_t iComb = 1;
  int i = 0;
  while (i < iR) {
      ++i;
      iComb *= iN - i + 1;
      iComb /= i;
  }
  return iComb;
}

vector<bool> idx2bit(size_t idx, int nsize, int nocc) {
  vector<bool> bit(nsize, false);
  for (int i = 0; i < nsize; ++i) {
    if (idx >= choose(nsize-1-i, nocc)) {
      bit[i] = true;
      idx -= choose(nsize-1-i, nocc--);
    }
  }
  return bit;
}

size_t bit2idx(const vector<bool>& bit) {
  int occ = 0, nsize = bit.size();
  size_t idx = 0;
  for (int i = 0; i < nsize; ++i) {
    if (bit[nsize-1-i]) {
      idx += choose(i, ++occ);
    }
  }
  return idx;
}

d_real weight_bound(const vector<d_real>& w, int n) {
  vector<d_real> temp(w);
  std::sort(temp.begin(), temp.end());
  d_real max = 1.;
  for (int i = 0; i < temp.size(); ++i) {
    max *= (i < temp.size()-n) ? (1.-temp[i]) : temp[i];
  }
  return max;
}

map<vector<bool>, d_real> simple_combinations(const vector<d_real>& w, int n, d_real thr) {
  // simple function to calculate combination weight, keep only the ones larger than threshold
  int nsize = w.size();
  size_t max = choose(nsize, n);
  map<vector<bool>, d_real> wtable;

  for (size_t i = 0; i < max; ++i) {
    auto bits = idx2bit(i, nsize, n);
    d_real weight = 1.;
    for (int j = 0; j < nsize; ++j) {
      weight *= bits[j] ? w[j] : (1-w[j]);
    }
    if (weight > thr) {
      wtable.insert(std::pair<vector<bool>, d_real>(bits, weight));
    }
  }
  return wtable;
}

map<vector<bool>, d_real> combinations(const vector<d_real>& w, int n, d_real thr) {
  if (weight_bound(w, n) < thr) {
    map<vector<bool>, d_real> wtable;
    return wtable;
  }
  
  if (choose(w.size(), n) < 10000) {
    return simple_combinations(w, n, thr);
  } else { // if too big, use recursion method to compute
    vector<d_real> w1, w2;
    map<vector<bool>, d_real> wtable;    
    int nsize1 = w.size()/2, nsize2 = w.size()-nsize1;

    copy(w.begin(), w.begin() + nsize1, std::back_inserter(w1));
    copy(w.begin() + nsize1, w.end(), std::back_inserter(w2));
    
    for (int k = 0; k <= n; ++k) {
      if (k <= nsize1 && n-k >= 0 && n-k <= nsize2) {
        auto comb1 = combinations(w1, k, thr / weight_bound(w2, n-k));
        auto comb2 = combinations(w2, n-k, thr / weight_bound(w1, k));
        for (auto it1 = comb1.cbegin(); it1 != comb1.cend(); ++it1) {
          for (auto it2 = comb2.cbegin(); it2 != comb2.cend(); ++it2) {
            if (it1->second * it2->second > thr) {
              vector<bool> merge = it1->first;
              merge.insert(merge.end(), it2->first.begin(), it2->first.end());
              wtable.insert(std::pair<vector<bool>, d_real>(merge, it1->second * it2->second));
            }
          }
        }
      }
    }
    return wtable;
  }
}

ActiveSpaceIterator::ActiveSpaceIterator(int q, const SchmidtBasis* _basis): quantum(q), basis(_basis), weight(_basis -> lweight()), nsize(_basis -> nactive()) {}

ActiveSpaceIterator_Slater::ActiveSpaceIterator_Slater(int q, const SchmidtBasis* _basis): ActiveSpaceIterator(q, _basis) {
  int ntotal = basis -> lsites() - basis -> nlcore() * 2;
  na = (ntotal + q) / 2;
  nb = (ntotal - q) / 2;
  build_iterator();
}

ActiveSpaceIterator_BCS::ActiveSpaceIterator_BCS(int q, const SchmidtBasis* _basis): ActiveSpaceIterator(q, _basis) {
  nqp = q - basis -> nlcore() + basis -> lsites(); // nqp + ncore - nsites = q
  build_iterator();  
}

void ActiveSpaceIterator_Slater::build_iterator() {
  auto wtable_a = combinations(weight, na, params.thrnp / weight_bound(weight, nb));
  auto wtable_b = combinations(weight, nb, params.thrnp / weight_bound(weight, na));

  for (auto it1 = wtable_a.begin(); it1 != wtable_a.end(); ++it1) {
    for (auto it2 = wtable_b.begin(); it2 != wtable_b.end(); ++it2) {
      if (it1->second * it2->second > params.thrnp) {
        vector<bool> merge = it1->first;
        merge.insert(merge.end(), it2->first.begin(), it2->first.end());
        l_list.push_back(merge);
        r_list.push_back(complimentary(merge));
        npweight.push_back(it1->second * it2->second);
      }
    }
  }
}

void ActiveSpaceIterator_BCS::build_iterator() {
  auto wtable = combinations(weight, nqp, params.thrnp);
  l_list.resize(wtable.size());
  r_list.resize(wtable.size());
  npweight.resize(wtable.size());
  int count = 0;
  for (auto it = wtable.begin(); it != wtable.end(); ++it) {
    l_list[count] = it->first;
    r_list[count] = complimentary(it->first);
    npweight[count] = it->second;
    ++count;
  }
}

std::ostream& operator <<(std::ostream& os, const SchmidtBasis& basis) {
  os.setf(std::ios::fixed, std::ios::floatfield);
  os.precision(10);
  os << "Left Sites " << basis.lsites() << "  Right Sites " << basis.rsites() << endl;
  os << "Core Orbitals  Left (" << basis.nlcore() << ")  Right (" << basis.nrcore() << ")" << endl;
  os << "Active Space (" << basis.nactive() << ") with Left Weights:\n";
  for (int i = 0; i < basis.nactive(); ++i) {
    if (i > 0 && i % 6 == 0) {
      os << endl;
    }
    os << basis.m_lweight[i] << "  ";
  }
  os << endl;
  os << "Quantum Numbers  " << basis.quantums << endl;
  os << "Dimensions       " << basis.dims << endl;
  return os;
}

void SchmidtBasis_Slater::calc_q() {
  int nelec = m_lsites - 2 * m_nlc; // in active space
  for (int nelec_a = 0; nelec_a <= nelec; ++nelec_a) {
    if (weight_bound(m_lweight, nelec_a) * weight_bound(m_lweight, nelec-nelec_a) > params.thrnp) {
      quantums.push_back(nelec_a*2 - nelec);
    }
  }
}

void SchmidtBasis_BCS::calc_q() {
  for (int n = m_nlc; n <= m_nlc + m_na; ++n) {
    if (n % 2 == 0 && weight_bound(m_lweight, n - m_nlc) > params.thrnp) {
      quantums.push_back(n - m_lsites);  // we always consider left
    }
  }
}

void SchmidtBasis::calc_dim() {
  dims.resize(quantums.size(), 0);
  for (int i = 0; i < quantums.size(); ++i) {
    dims[i] = iterator(quantums[i]) -> size();
  }
}

boost::shared_ptr<ActiveSpaceIterator> SchmidtBasis::iterator(int q) {
  auto it = m_it.find(q);
  if (it == m_it.end()) {
    auto ptr_asi = params.bcs ? 
      boost::shared_ptr<ActiveSpaceIterator>(new ActiveSpaceIterator_BCS(q, this)):
      boost::shared_ptr<ActiveSpaceIterator>(new ActiveSpaceIterator_Slater(q, this));
    m_it.insert(std::pair<int, boost::shared_ptr<ActiveSpaceIterator>>(q, ptr_asi));
  }
  return m_it[q];
}

std::ostream& operator << (std::ostream& os, const vector<bool>& bits) {
  for (int i = 0; i < bits.size(); ++i) {
    os << (bits[i] ? 1 : 0);
  }
  return os;
}
