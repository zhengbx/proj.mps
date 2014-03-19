#ifndef DM
#define DM

#include "utils.h"
#include "schmidt.h"
#include <memory>

class DensityMatrix {
protected:
  Matrix m_coefs;

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version) {
    ar & m_coefs;
  }
public:
  DensityMatrix(Matrix& coefs) {
    m_coefs = std::move(coefs);
  }
  virtual std::shared_ptr<SchmidtBasis> basis(int cut) const = 0;
  virtual int get_nsites() const = 0;
  virtual int get_norbs() const = 0;
};

class SlaterDM: public DensityMatrix {
protected:
  int nsites, norbs;

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version) {
    // serialize base class information
    ar & boost::serialization::base_object<DensityMatrix>(*this);
    ar & nsites & norbs;
  }
public:
  SlaterDM(Matrix& coefs): DensityMatrix(coefs), 
      nsites(m_coefs.Nrows()), norbs(m_coefs.Ncols()) {};
  std::shared_ptr<SchmidtBasis> basis(int cut) const;
  int get_nsites() const {  return nsites;};
  int get_norbs() const {  return norbs;};
};

class BCSDM: public DensityMatrix {
protected:
  int nsites;

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version) {
    // serialize base class information
    ar & boost::serialization::base_object<DensityMatrix>(*this);
    ar & nsites;
  }
public:
  BCSDM(Matrix& coefs): DensityMatrix(coefs), nsites(coefs.Ncols()) {}
  std::shared_ptr<SchmidtBasis> basis(int cut) const;
  int get_nsites() const {  return nsites;};
  int get_norbs() const {  return nsites;};

};

#endif
