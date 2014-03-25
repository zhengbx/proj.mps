#ifndef SCHMIDT
#define SCHMIDT

#include "include.h"
#include "utils.h"
#include <memory>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/map.hpp>
#include <map>

using std::map;

class SchmidtBasis;

class ActiveSpaceIterator {
protected:
  int quantum, nsize;
  const SchmidtBasis* basis;
  vector<vector<bool>> l_list, r_list;
  vector<double> npweight;
  vector<double> weight;

  friend class boost::serialization::access;
  ActiveSpaceIterator() {}  
  template<class Archive> void serialize(Archive & ar, const unsigned int version) {
    ar & quantum & nsize & l_list & r_list & weight & npweight;
  }
public:
  ActiveSpaceIterator(int q, const SchmidtBasis* _basis);
  void set_basis(const SchmidtBasis* _basis) {
    basis = _basis;
  }
  ~ActiveSpaceIterator() {
    basis = nullptr;
  }
  int size() const {  return l_list.size();}
  vector<bool> get_config(const int& i, bool left = true) const {  
    return std::move(left ? l_list[i]: r_list[i]);
  }
};

class ActiveSpaceIterator_Slater: public ActiveSpaceIterator {
protected:
  int na, nb; // number of alpha and beta electrons
  void build_iterator();

  friend class boost::serialization::access;
  ActiveSpaceIterator_Slater() {}  
  template<class Archive> void serialize(Archive & ar, const unsigned int version) {
    ar & boost::serialization::base_object<ActiveSpaceIterator>(*this);
    ar & na & nb;
  }
public:
  ActiveSpaceIterator_Slater(int q, const SchmidtBasis* _basis);
};


class ActiveSpaceIterator_BCS: public ActiveSpaceIterator {
protected:
  int nqp; // number of quasiparticle excitations
  void build_iterator();

  friend class boost::serialization::access;
  ActiveSpaceIterator_BCS() {}  
  template<class Archive> void serialize(Archive & ar, const unsigned int version) {
    ar & boost::serialization::base_object<ActiveSpaceIterator>(*this);
    ar & nqp;
  }
public:
  ActiveSpaceIterator_BCS(int q, const SchmidtBasis* _basis);
};


class SchmidtBasis {
protected:
  int m_lsites, m_rsites;
  int m_nlc, m_nrc, m_na;
  Matrix m_lc, m_rc, m_la, m_ra;
  vector<double> m_lweight, m_rweight;
  vector<int> quantums, dims;
  map<int, boost::shared_ptr<ActiveSpaceIterator>> m_it;

  
  // compute possible quantum numbers, for left basis,
  // if use for right basis, just add a minus sign
  virtual void calc_q() = 0;
  void calc_dim();

  friend class boost::serialization::access;
  SchmidtBasis() {}  
  template<class Archive> void serialize(Archive & ar, const unsigned int version) {
    ar & m_nlc & m_nrc & m_na & m_lsites & m_rsites;
    ar & m_lc & m_rc & m_la & m_ra & quantums & dims;
    ar & m_lweight & m_rweight;

    ar.register_type(static_cast<ActiveSpaceIterator_Slater*>(nullptr));
    ar.register_type(static_cast<ActiveSpaceIterator_BCS*>(nullptr));
    ar & m_it;
    for (auto iter = m_it.begin(); iter != m_it.end(); ++iter) {
      iter -> second -> set_basis(this);
    }
  }
public:
  SchmidtBasis(const Matrix& lcore, const Matrix& rcore, const Matrix& lactive, 
      const Matrix& ractive, const vector<double>& lweight): m_lc(lcore), m_rc(rcore), 
      m_la(lactive), m_ra(ractive), m_lweight(lweight), m_nlc(lcore.Ncols()),
      m_nrc(rcore.Ncols()), m_na(lactive.Ncols()) {
    m_rweight.resize(m_na);
    for (int i = 0; i < m_na; ++i) {
      m_rweight[i] = 1. - m_lweight[i];
    }
  }

  int lsites() const {  return m_lsites;};
  int rsites() const {  return m_rsites;};

  int nlcore() const {  return m_nlc;};
  int nrcore() const {  return m_nrc;};

  int nactive() const {  return m_na;};

  Matrix lcore() const {  return std::move(m_lc);};
  Matrix rcore() const {  return std::move(m_rc);};
  Matrix lactive() const {  return std::move(m_la);};
  Matrix ractive() const {  return std::move(m_ra);};

  const vector<int> get_q() const {  return std::move(quantums);}
  const vector<int> get_d() const {  return std::move(dims);}

  vector<double> lweight() const {
    return std::move(m_lweight);
  }

  friend std::ostream& operator << (std::ostream& os, const SchmidtBasis& basis);
  boost::shared_ptr<ActiveSpaceIterator> iterator(int q); // spin of left block, use -q if with right block
};

class SchmidtBasis_Slater: public SchmidtBasis {
protected:
  friend class boost::serialization::access;
  SchmidtBasis_Slater() {}    
  template<class Archive> void serialize(Archive & ar, const unsigned int version) {
    ar & boost::serialization::base_object<SchmidtBasis>(*this);
  }
public:
  SchmidtBasis_Slater(const Matrix& lcore, const Matrix& rcore, const Matrix& lactive, 
      const Matrix& ractive, const vector<double>& lweight): SchmidtBasis(lcore, rcore, 
        lactive, ractive, lweight) {
    m_lsites = lcore.Nrows();
    m_rsites = rcore.Nrows();
    calc_q();
    calc_dim();
  }
  void calc_q();
};

class SchmidtBasis_BCS: public SchmidtBasis {
protected:
  friend class boost::serialization::access;
  SchmidtBasis_BCS() {}    
  template<class Archive> void serialize(Archive & ar, const unsigned int version) {
    ar & boost::serialization::base_object<SchmidtBasis>(*this);
  }
public:
  SchmidtBasis_BCS(const Matrix& lcore, const Matrix& rcore, const Matrix& lactive, 
      const Matrix& ractive, const vector<double>& lweight): SchmidtBasis(lcore, rcore, 
        lactive, ractive, lweight) {
    m_lsites = lcore.Nrows()/2;
    m_rsites = rcore.Nrows()/2;
    calc_q();
    calc_dim();
  }
  void calc_q();
};

namespace boost { 
namespace serialization {
template <class Archive> void serialize(Archive & ar, boost::shared_ptr<SchmidtBasis>& t, const unsigned int file_version) {
  ar.register_type(static_cast<SchmidtBasis_Slater*>(nullptr));
  ar.register_type(static_cast<SchmidtBasis_BCS*>(nullptr));
  BOOST_STATIC_ASSERT(boost::serialization::tracking_level<SchmidtBasis>::value != boost::serialization::track_never);
  boost::serialization::split_free(ar, t, file_version);
}
}
}

std::ostream& operator << (std::ostream& os, const vector<bool>& bits);
#endif
