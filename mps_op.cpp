#include "mps_op.h"
#include <fstream>
#include <unistd.h>
#include <boost/filesystem.hpp>
#include "utils.h"

using std::ofstream;
using std::ifstream;

size_t mps_size(const QSTArray<dtype, 3, Quantum>& A) {
  size_t size = 0;
  for (auto it = A.begin(); it != A.end(); ++it) {
    size += it -> second -> size();
  }
  return size;
}

void save_site(const QSTArray<dtype, 3, Quantum>& A, int site, const char *filename){
  char name[50];
  sprintf(name,"%s/%d.mps",filename,site);
  ofstream fout(name);
  boost::archive::binary_oarchive oar(fout);
  oar << A;
}

void load_site(QSTArray<dtype, 3, Quantum>& A, int site ,const char *filename){
  char name[50];
  sprintf(name,"%s/%d.mps",filename,site);
  ifstream fin(name);
  boost::archive::binary_iarchive iar(fin);
  iar >> A;
}

void save_site(const MPS<dtype, Quantum>& mps, int site, const char *filename){
  save_site(mps[site], site, filename);
}

void load_site(MPS<dtype, Quantum>& mps, int site ,const char *filename){
  load_site(mps[site], site, filename);
}

d_real norm_on_disk(const char* filename, int size) {
  QSTArray<dtype, 2, Quantum> E;
  MPS<dtype, Quantum> mps(size);
  load_site(mps, 0, filename);
  Contract((dtype)1.0,mps[0],shape(0,1),mps[0].conjugate(),shape(0,1),(dtype)0.0,E);
  mps[0].clear();
  // intermediate
  QSTArray<dtype, 3, Quantum> I;
  for(int i = 1; i < size; ++i){
    load_site(mps, i, filename);
    //construct intermediate, i.e. past X to E
    Contract((dtype)1.0,E,shape(0),mps[i],shape(0),(dtype)0.0,I);
    //clear structure of E
    E.clear();
    //construct E for site i by contracting I with Y
    Contract((dtype)1.0,I,shape(0,1),mps[i].conjugate(),shape(0,1),(dtype)0.0,E);
    mps[i].clear();
    I.clear();
    //bad style: if no blocks remain, return zero
    if(E.begin() == E.end()) {
      return (d_real)0.0;
    }
  }
#ifdef _COMPLEX
  return abs(sqrt((*(E.find(shape(0,0))->second))(0,0)));
#else
  return sqrt((*(E.find(shape(0,0))->second))(0,0));
#endif
}


void normalize_on_disk(const char* filename, int size, int site) {
  d_real norm = norm_on_disk(filename, size);

  MPS<dtype, Quantum> mps(size);  
  if (site < 0) {
    d_real alpha = pow(1./norm, 1./(d_real)size);
    for (int i = 0; i < size; ++i) {
      load_site(mps, i, filename);
      Scal(alpha, mps[i]);
      save_site(mps, i, filename);
      mps[i].clear();
    }
  } else {
    load_site(mps, site, filename);
    Scal((d_real)1./norm, mps[site]);
    save_site(mps, site, filename);
    mps[site].clear();
  }
}

void check_existence(int i, const char* filename) {
  char name[50];
  sprintf(name,"%s/%d.mps",filename,i);
  while (!boost::filesystem::exists(string(name))) {
    sleep(1);
  }
  sleep(1);
}

typename remove_complex<dtype>::type partial_compress(int L,const MPS_DIRECTION& dir, int D, const char* filename, int first, int last) {
  d_real dweight = 0.0;

  MPS<dtype, Quantum> mps(L);
  dtype acc_norm = 1.;

  if(dir == MPS_DIRECTION::Left) {
    STArray<d_real, 1> S;//singular values
    QSTArray<dtype, 2, Quantum> V;//V^T
    QSTArray<dtype, 3, Quantum> U;//U --> unitary left normalized matrix
    check_existence(first, filename);

    load_site(mps, first, filename);  // load first site
    printf("site %3d before compress %12d\n", first, mps_size(mps[first]));
    //redistribute the norm over the chain: for stability reasons
    // note this is not complete, the sites on the left are not affected
    // need final renormalization
    dtype nrm = sqrt(Dotc(mps[first],mps[first]));
    acc_norm *= pow(nrm, 1./(d_real)L);
    Scal(acc_norm / nrm, mps[first]);
    
    for (int i = first; i < last; ++i) {
      //svd
      dweight += Gesvd<dtype,3,3,Quantum,btas::RightArrow>(mps[i],S,U,V,D);
      //copy unitary to mpx
      Copy(U,mps[i]);
      printf("site %3d  after compress %12d\n", i, mps_size(mps[i]));      
      save_site(mps, i, filename);
      mps[i].clear();
      //paste S and V together
      Dimm(S,V);
      // now read next site
      check_existence(i+1, filename);
      load_site(mps, i+1, filename);
      printf("site %3d before compress %12d\n", i+1, mps_size(mps[i+1]));
      //and multiply with mpx on the next site
      U = mps[i + 1];
      //when compressing dimensions will change, so reset:
      mps[i + 1].clear();
      Contract((dtype)1.0,V,shape(1),U,shape(0),(dtype)0.0,mps[i + 1]);
      dtype nrm = sqrt(Dotc(mps[i+1],mps[i+1]));
      acc_norm *= pow(nrm, 1./(d_real)L);
      Scal(acc_norm / nrm, mps[i+1]);
    }
    printf("site %3d not compressed\n", last);
    save_site(mps, last, filename);
    mps[last].clear();
  } else {
    STArray<dtype, 1> S;//singular values
    QSTArray<dtype, 3, Quantum> V;//V^T --> unitary right normalized matrix
    QSTArray<dtype, 2, Quantum> U;//U
    check_existence(first, filename);
    load_site(mps, first, filename);  // load first site
    printf("site %3d before compress %12d\n", first, mps_size(mps[first]));
    //redistribute the norm over the chain: for stability reasons
    dtype nrm = sqrt(Dotc(mps[first],mps[first]));
    acc_norm *= pow(nrm, 1./(d_real)L);
    Scal(acc_norm/nrm,mps[first]);

    for(int i = first;i > last;--i){
      //then SVD: 
      dweight += Gesvd<dtype, 3, 2, Quantum, btas::RightArrow>(mps[i],S,U,V,D);
      //copy unitary to mpx
      Copy(V,mps[i]);
      printf("site %3d  after compress %12d\n", i, mps_size(mps[i]));
      save_site(mps, i, filename);
      mps[i].clear();
      //paste U and S together
      Dimm(U,S);
      check_existence(i-1, filename);
      load_site(mps, i-1, filename);
      printf("site %3d before compress %12d\n", i-1, mps_size(mps[i-1]));
      //and multiply with mpx on the next site
      V = mps[i - 1];
      //when compressing dimensions will change, so reset:
      mps[i - 1].clear();
      Contract((dtype)1.0,V,shape(2),U,shape(0),(dtype)0.0,mps[i - 1]);

      dtype nrm = sqrt(Dotc(mps[i-1],mps[i-1]));
      acc_norm *= pow(nrm, 1./(d_real)L);
      Scal(acc_norm/nrm,mps[i-1]);
    }
    printf("site %3d not compressed\n", last);
    save_site(mps, last, filename);
    mps[last].clear();
  }
  // attention normalization is needed later!
}

typename remove_complex<dtype>::type compress_on_disk(int L,const MPS_DIRECTION &dir,int D, const char *filename, bool store, bool incomp){

  d_real dweight = 0.0;

  MPS<dtype, Quantum> mps(L);
  dtype acc_norm = 1.;

  if(dir == MPS_DIRECTION::Left) {
    STArray<d_real, 1> S;//singular values
    QSTArray<dtype, 2, Quantum> V;//V^T
    QSTArray<dtype, 3, Quantum> U;//U --> unitary left normalized matrix
    if (incomp) check_existence(0, filename);
    load_site(mps, 0, filename);  // load first site
    cout << "site 0 before compress " << mps_size(mps[0]) << endl;
    for(int i = 0;i < L - 1;++i){
      //redistribute the norm over the chain: for stability reasons
      // note this is not complete, the sites on the left are not affected
      // need final renormalization
      dtype nrm = sqrt(Dotc(mps[i],mps[i]));
      acc_norm *= pow(nrm, 1./(d_real)L);
      Scal(acc_norm / nrm, mps[i]);
      
      //then svd
      dweight += Gesvd<dtype,3,3,Quantum,btas::RightArrow>(mps[i],S,U,V,D);
      //copy unitary to mpx
      Copy(U,mps[i]);
      cout << "site " << i << " after compress  " << mps_size(mps[i]) << endl;      
      if (store) {
        save_site(mps, i, filename);
        mps[i].clear();
      }

      //paste S and V together
      Dimm(S,V);
      // now read next site
      if (incomp) check_existence(i+1, filename);
      load_site(mps, i+1, filename);
      cout << "site " << i+1 << " before compress " << mps_size(mps[i+1]) << endl;      
      //and multiply with mpx on the next site
      U = mps[i + 1];
      //when compressing dimensions will change, so reset:
      mps[i + 1].clear();
      Contract((dtype)1.0,V,shape(1),U,shape(0),(dtype)0.0,mps[i + 1]);
    }
    dtype nrm = sqrt(Dotc(mps[L-1],mps[L-1]));
    acc_norm *= pow(nrm, 1./(d_real)L);
    Scal(acc_norm/nrm,mps[L-1]);
    cout << "site " << L-1 << " after compress  " << mps_size(mps[L-1]) << endl;      
    if (store) {
      save_site(mps, L-1, filename);
      mps[L-1].clear();
    }
  } else {//right
    STArray<d_real, 1> S;//singular values
    QSTArray<dtype, 3, Quantum> V;//V^T --> unitary right normalized matrix
    QSTArray<dtype, 2, Quantum> U;//U
    if (incomp) check_existence(L-1, filename);
    load_site(mps, L-1, filename);  // load first site
    cout << "site " << L-1 << " before compress " << mps_size(mps[L-1]) << endl;
    for(int i = L - 1;i > 0;--i){
      //redistribute the norm over the chain: for stability reasons
      dtype nrm = sqrt(Dotc(mps[i],mps[i]));
      acc_norm *= pow(nrm, 1./(d_real)L);      
      Scal(acc_norm/nrm,mps[i]);
      //then SVD: 
      dweight += Gesvd<dtype, 3, 2, Quantum, btas::RightArrow>(mps[i],S,U,V,D);

      //copy unitary to mpx
      Copy(V,mps[i]);
      cout << "site" << i << " after compress  " << mps_size(mps[i]) << endl;
      if (store) {
        save_site(mps, i, filename);
        mps[i].clear();
      }
      //paste U and S together
      Dimm(U,S);

      // now read next site
      if (incomp) check_existence(i-1, filename);
      load_site(mps, i-1, filename);
      cout << "site" << i-1 << " before compress " << mps_size(mps[i-1]) << endl;
      //and multiply with mpx on the next site
      V = mps[i - 1];
      //when compressing dimensions will change, so reset:
      mps[i - 1].clear();
      Contract((dtype)1.0,V,shape(2),U,shape(0),(dtype)0.0,mps[i - 1]);
    }
    dtype nrm = sqrt(Dotc(mps[0],mps[0]));
    acc_norm *= pow(nrm, 1./(d_real)L);    
    Scal(acc_norm/nrm,mps[0]);
    cout << "site " << 0 << " after compress  " << mps_size(mps[0]) << endl;      
    if (store) {
      save_site(mps, 0, filename);
      mps[0].clear();
    }
  }
  // now normalize all the sites
  if (store) {
    normalize_on_disk(filename, L);
  } else {
    normalize(mps);
  }
}

// process entaglement spectra
void spectra(const tuple<STArray<typename remove_complex<dtype>::type, 1> , Qshapes<Quantum>>& raw) {
  STArray<d_real, 1> sc = std::get<0>(raw);
  Qshapes<Quantum> sq = std::get<1>(raw);

  auto iter = sc.begin();  
  map<int, vector<d_real>> coef;
  for (int i = 0; i < sq.size(); i++) {
    int sp = sq[i].gSz();
    coef[sp];
    for (auto it_d = iter -> second -> begin(); it_d != iter -> second -> end(); ++it_d) {
      coef[sp].push_back(-2.*log(*it_d));
    }
    std::sort(coef[sp].begin(), coef[sp].end());
    ++iter;    
  }
  ofstream fspec((params.path + "/spectra.out").c_str());
  fspec.setf(std::ios::fixed, std::ios::floatfield);
  fspec.precision(10);
  for (auto it = coef.begin(); it != coef.end(); it++) {
    fspec << "Section S=" << it -> first << endl;
    int count = 0;
    for (auto it_d = it -> second.begin(); it_d < it -> second.end(); ++it_d) {
      if (count > 0 && count % 6 == 0) {
        fspec << endl;
      }
      fspec << *it_d << "\t";
      ++count;
    }
    fspec << endl;
  }
  fspec.close();
}

tuple<STArray<typename remove_complex<dtype>::type, 1>, Qshapes<Quantum>> Schmidt_on_disk(int L, int site, const char* filename, int lc, int rc) {
  MPS<dtype, Quantum> mps(L);
  vector<d_real> coef;
  if (lc < 0) {
    lc = 0;
  }
  if (rc < 0) {
    rc = L-1;
  }
  if (site < 0) {
    site = L/2;
  }
  
  if (lc < site) {
    load_site(mps, lc, filename);
    for(int i = lc;i < site;++i){
      STArray<d_real, 1> S;//singular values
      QSTArray<dtype, 2, Quantum> V;//V^T
      QSTArray<dtype, 3, Quantum> U;//U --> unitary left normalized matrix

      //then svd
      Gesvd<dtype, 3, 3, Quantum, btas::RightArrow>(mps[i],S,U,V,0);
      //copy unitary to mps
      Copy(U,mps[i]);
      save_site(mps, i, filename);
      mps[i].clear();

      //paste S and V together
      Dimm(S,V);
      // now read next site
      load_site(mps, i+1, filename);
      //and multiply with mpx on the next site
      U = mps[i + 1];
      //when compressing dimensions will change, so reset:
      mps[i + 1].clear();
      Contract((dtype)1.0,V,shape(1),U,shape(0),(dtype)0.0,mps[i + 1]);
    }
    // save matrix[site]
    save_site(mps, site, filename);
    mps[site].clear();
  }

  if (rc > site) {
    load_site(mps, rc, filename);  // load first site
    for(int i = rc;i > site;--i){
      STArray<d_real, 1> S;//singular values
      QSTArray<dtype, 3, Quantum> V;//V^T --> unitary right normalized matrix
      QSTArray<dtype, 2, Quantum> U;//U
      //then SVD: 
      Gesvd<dtype, 3, 2, Quantum, btas::RightArrow>(mps[i],S,U,V,0);
      //copy unitary to mpx
      Copy(V,mps[i]);
      save_site(mps, i, filename);
      mps[i].clear();

      //paste U and S together
      Dimm(U,S);
      // now read next site
      load_site(mps, i-1, filename);
      //and multiply with mpx on the next site
      V = mps[i - 1];
      //when compressing dimensions will change, so reset:
      mps[i - 1].clear();
      Contract((dtype)1.0,V,shape(2),U,shape(0),(dtype)0.0,mps[i - 1]);
    }
    // save matrix[site]
    save_site(mps, site, filename);
    mps[site].clear();
  }
  
  load_site(mps, site, filename);  
  STArray<d_real, 1> S;//singular values
  QSTArray<dtype, 3, Quantum> V;//V^T
  QSTArray<dtype, 2, Quantum> U;//U --> unitary left normalized matrix
  Gesvd<dtype, 3, 2, Quantum, btas::RightArrow>(mps[site],S,U,V,0);
  // save matrix[site]
  save_site(mps, site, filename);
  mps[site].clear();
  return std::make_tuple(S, U.qshape()[0]);
}




