#include "mps_op.h"
#include <fstream>
#include <unistd.h>
#include <boost/filesystem.hpp>
#include "utils.h"

using std::ofstream;
using std::ifstream;

size_t mps_size(const QSDArray<3, Quantum>& A) {
  size_t size = 0;
  for (auto it = A.begin(); it != A.end(); ++it) {
    size += it -> second -> size();
  }
  return size;
}

void save_site(const QSDArray<3, Quantum>& A, int site, const char *filename){
  char name[50];
  sprintf(name,"%s/%d.mps",filename,site);
  ofstream fout(name);
  boost::archive::binary_oarchive oar(fout);
  oar << A;
}

void load_site(QSDArray<3, Quantum>& A, int site ,const char *filename){
  char name[50];
  sprintf(name,"%s/%d.mps",filename,site);
  ifstream fin(name);
  boost::archive::binary_iarchive iar(fin);
  iar >> A;
}

void save_site(const MPS<Quantum>& mps, int site, const char *filename){
  save_site(mps[site], site, filename);
}

void load_site(MPS<Quantum>& mps, int site ,const char *filename){
  load_site(mps[site], site, filename);
}

double norm_on_disk(MPS<Quantum> &mps, const char* filename) {
  QSDArray<2> E;
  int L = mps.size();
  load_site(mps, 0, filename);
  QSDcontract(1.0,mps[0],shape(0,1),mps[0].conjugate(),shape(0,1),0.0,E);
  mps[0].clear();
  // intermediate
  QSDArray<3> I;
  for(int i = 1; i < L; ++i){
    load_site(mps, i, filename);
    //construct intermediate, i.e. past X to E
    QSDcontract(1.0,E,shape(0),mps[i],shape(0),0.0,I);
    //clear structure of E
    E.clear();
    //construct E for site i by contracting I with Y
    QSDcontract(1.0,I,shape(0,1),mps[i].conjugate(),shape(0,1),0.0,E);
    mps[i].clear();
    I.clear();
    //bad style: if no blocks remain, return zero
    if(E.begin() == E.end()) {
      return 0.0;
    }
  }
  return sqrt((*(E.find(shape(0,0))->second))(0,0));
}


void normalize_on_disk(MPS<Quantum>& mps, const char* filename, int site) {
  double norm = norm_on_disk(mps, filename);
  if (site < 0) {
    int L = mps.size();
    double alpha = pow(1./norm, 1./(double)L);
    for (int i = 0; i < L; ++i) {
      load_site(mps, i, filename);
      QSDscal(alpha, mps[i]);
      save_site(mps, i, filename);
      mps[i].clear();
    }
  } else {
    load_site(mps, site, filename);
    QSDscal(1./norm, mps[site]);
    save_site(mps, site, filename);
    mps[site].clear();
  }
}

void check_existence(int i, const char* filename) {
  char name[50];
  sprintf(name,"%s/%d.mps",filename,i);
  while (!boost::filesystem::exists(string(name))) {
    sleep(5);
  }
  sleep(5);
}

void partial_compress(int L,const MPS_DIRECTION& dir, int D, const char* filename, int last) {
  MPS<Quantum> mps(L);
  double acc_norm = 1.;

  if(dir == MPS_DIRECTION::Left) {
    SDArray<1> S;//singular values
    QSDArray<2> V;//V^T
    QSDArray<3> U;//U --> unitary left normalized matrix
    check_existence(0, filename);

    load_site(mps, 0, filename);  // load first site
    printf("site %3d before compress %12d\n", 0, mps_size(mps[0]));
    //redistribute the norm over the chain: for stability reasons
    // note this is not complete, the sites on the left are not affected
    // need final renormalization
    double nrm = sqrt(QSDdotc(mps[0],mps[0]));
    acc_norm *= pow(nrm, 1./(double)L);
    QSDscal(acc_norm / nrm, mps[0]);
    
    for (int i = 0; i < last; ++i) {
      //svd
      QSDgesvd(RightArrow,mps[i],S,U,V,D);
      //copy unitary to mpx
      QSDcopy(U,mps[i]);
      printf("site %3d  after compress %12d\n", i, mps_size(mps[i]));      
      save_site(mps, i, filename);
      mps[i].clear();
      //paste S and V together
      SDdidm(S,V);
      // now read next site
      check_existence(i+1, filename);
      load_site(mps, i+1, filename);
      printf("site %3d before compress %12d\n", i+1, mps_size(mps[i+1]));
      //and multiply with mpx on the next site
      U = mps[i + 1];
      //when compressing dimensions will change, so reset:
      mps[i + 1].clear();
      QSDcontract(1.0,V,shape(1),U,shape(0),0.0,mps[i + 1]);
      double nrm = sqrt(QSDdotc(mps[i+1],mps[i+1]));
      acc_norm *= pow(nrm, 1./(double)L);
      QSDscal(acc_norm / nrm, mps[i+1]);
    }
    printf("site %3d not compressed\n", last);
    save_site(mps, last, filename);
    mps[last].clear();
  } else {
    SDArray<1> S;//singular values
    QSDArray<3> V;//V^T --> unitary right normalized matrix
    QSDArray<2> U;//U
    check_existence(L-1, filename);
    load_site(mps, L-1, filename);  // load first site
    printf("site %3d before compress %12d\n", L-1, mps_size(mps[L-1]));
    //redistribute the norm over the chain: for stability reasons
    double nrm = sqrt(QSDdotc(mps[L-1],mps[L-1]));
    acc_norm *= pow(nrm, 1./(double)L);
    QSDscal(acc_norm/nrm,mps[L-1]);

    for(int i = L - 1;i > last;--i){
      //then SVD: 
      QSDgesvd(RightArrow,mps[i],S,U,V,D);
      //copy unitary to mpx
      QSDcopy(V,mps[i]);
      printf("site %3d  after compress %12d\n", i, mps_size(mps[i]));
      save_site(mps, i, filename);
      mps[i].clear();
      //paste U and S together
      SDdimd(U,S);
      check_existence(i-1, filename);
      load_site(mps, i-1, filename);
      printf("site %3d before compress %12d\n", i-1, mps_size(mps[i-1]));
      //and multiply with mpx on the next site
      V = mps[i - 1];
      //when compressing dimensions will change, so reset:
      mps[i - 1].clear();
      QSDcontract(1.0,V,shape(2),U,shape(0),0.0,mps[i - 1]);

      double nrm = sqrt(QSDdotc(mps[i-1],mps[i-1]));
      acc_norm *= pow(nrm, 1./(double)L);
      QSDscal(acc_norm/nrm,mps[i-1]);
    }
    printf("site %3d not compressed\n", last);
    save_site(mps, last, filename);
    mps[last].clear();
  }
  // attention normalization is needed later!
}



void compress_on_disk(int L,const MPS_DIRECTION &dir,int D, const char *filename, bool store, bool on_the_fly){
  MPS<Quantum> mps(L);
  double acc_norm = 1.;

  if(dir == MPS_DIRECTION::Left) {
    SDArray<1> S;//singular values
    QSDArray<2> V;//V^T
    QSDArray<3> U;//U --> unitary left normalized matrix
    if (on_the_fly) check_existence(0, filename);
    load_site(mps, 0, filename);  // load first site
    cout << "site 0 before compress " << mps_size(mps[0]) << endl;
    for(int i = 0;i < L - 1;++i){
      //redistribute the norm over the chain: for stability reasons
      // note this is not complete, the sites on the left are not affected
      // need final renormalization
      double nrm = sqrt(QSDdotc(mps[i],mps[i]));
      acc_norm *= pow(nrm, 1./(double)L);
      QSDscal(acc_norm / nrm, mps[i]);
      
      //then svd
      QSDgesvd(RightArrow,mps[i],S,U,V,D);
      //copy unitary to mpx
      QSDcopy(U,mps[i]);
      cout << "site " << i << " after compress  " << mps_size(mps[i]) << endl;      
      if (store) {
        save_site(mps, i, filename);
        mps[i].clear();
      }

      //paste S and V together
      SDdidm(S,V);
      // now read next site
      if (on_the_fly) check_existence(i+1, filename);
      load_site(mps, i+1, filename);
      cout << "site " << i+1 << " before compress " << mps_size(mps[i+1]) << endl;      
      //and multiply with mpx on the next site
      U = mps[i + 1];
      //when compressing dimensions will change, so reset:
      mps[i + 1].clear();
      QSDcontract(1.0,V,shape(1),U,shape(0),0.0,mps[i + 1]);
    }
    double nrm = sqrt(QSDdotc(mps[L-1],mps[L-1]));
    acc_norm *= pow(nrm, 1./(double)L);
    QSDscal(acc_norm/nrm,mps[L-1]);
    cout << "site " << L-1 << " after compress  " << mps_size(mps[L-1]) << endl;      
    if (store) {
      save_site(mps, L-1, filename);
      mps[L-1].clear();
    }
  } else {//right
    SDArray<1> S;//singular values
    QSDArray<3> V;//V^T --> unitary right normalized matrix
    QSDArray<2> U;//U
    if (on_the_fly) check_existence(L-1, filename);
    load_site(mps, L-1, filename);  // load first site
    cout << "site " << L-1 << " before compress " << mps_size(mps[L-1]) << endl;
    for(int i = L - 1;i > 0;--i){
      //redistribute the norm over the chain: for stability reasons
      double nrm = sqrt(QSDdotc(mps[i],mps[i]));
      acc_norm *= pow(nrm, 1./(double)L);      
      QSDscal(acc_norm/nrm,mps[i]);
      //then SVD: 
      QSDgesvd(RightArrow,mps[i],S,U,V,D);
      //copy unitary to mpx
      QSDcopy(V,mps[i]);
      cout << "site" << i << " after compress  " << mps_size(mps[i]) << endl;
      if (store) {
        save_site(mps, i, filename);
        mps[i].clear();
      }

      //paste U and S together
      SDdimd(U,S);
      // now read next site
      if (on_the_fly) check_existence(i-1, filename);
      load_site(mps, i-1, filename);
      cout << "site" << i-1 << " before compress " << mps_size(mps[i-1]) << endl;
      //and multiply with mpx on the next site
      V = mps[i - 1];
      //when compressing dimensions will change, so reset:
      mps[i - 1].clear();
      QSDcontract(1.0,V,shape(2),U,shape(0),0.0,mps[i - 1]);

    }
    double nrm = sqrt(QSDdotc(mps[0],mps[0]));
    acc_norm *= pow(nrm, 1./(double)L);    
    QSDscal(acc_norm/nrm,mps[0]);
    cout << "site " << 0 << " after compress  " << mps_size(mps[0]) << endl;      
    if (store) {
      save_site(mps, 0, filename);
      mps[0].clear();
    }
  }
  // now normalize all the sites
  if (store) {
    normalize_on_disk(mps, filename);
  } else {
    normalize(mps);
  }
}

// process entaglement spectra
void spectra(const tuple<SDArray<1>, Qshapes<Quantum>>& raw) {
  SDArray<1> sc = std::get<0>(raw);
  Qshapes<Quantum> sq = std::get<1>(raw);

  auto iter = sc.begin();  
  map<int, vector<double>> coef;
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

tuple<SDArray<1>, Qshapes<Quantum>> Schmidt_on_disk(MPS<Quantum>& mps, int site, const char* filename, int lc, int rc) {
  vector<double> coef;
  int L = mps.size();
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
      SDArray<1> S;//singular values
      QSDArray<2> V;//V^T
      QSDArray<3> U;//U --> unitary left normalized matrix
      
      //then svd
      QSDgesvd(RightArrow,mps[i],S,U,V,0);
      //copy unitary to mps
      QSDcopy(U,mps[i]);
      save_site(mps, i, filename);
      mps[i].clear();

      //paste S and V together
      SDdidm(S,V);
      // now read next site
      load_site(mps, i+1, filename);
      //and multiply with mpx on the next site
      U = mps[i + 1];
      //when compressing dimensions will change, so reset:
      mps[i + 1].clear();
      QSDcontract(1.0,V,shape(1),U,shape(0),0.0,mps[i + 1]);
    }
    // save matrix[site]
    save_site(mps, site, filename);
    mps[site].clear();
  }

  if (rc > site) {
    load_site(mps, rc, filename);  // load first site
    for(int i = rc;i > site;--i){
      SDArray<1> S;//singular values
      QSDArray<3> V;//V^T --> unitary right normalized matrix
      QSDArray<2> U;//U
      //then SVD: 
      QSDgesvd(RightArrow,mps[i],S,U,V,0);
      //copy unitary to mpx
      QSDcopy(V,mps[i]);
      save_site(mps, i, filename);
      mps[i].clear();

      //paste U and S together
      SDdimd(U,S);
      // now read next site
      load_site(mps, i-1, filename);
      //and multiply with mpx on the next site
      V = mps[i - 1];
      //when compressing dimensions will change, so reset:
      mps[i - 1].clear();
      QSDcontract(1.0,V,shape(2),U,shape(0),0.0,mps[i - 1]);
    }
    // save matrix[site]
    save_site(mps, site, filename);
    mps[site].clear();
  }
  
  load_site(mps, site, filename);  
  SDArray<1> S;//singular values
  QSDArray<3> V;//V^T
  QSDArray<2> U;//U --> unitary left normalized matrix
  QSDgesvd(RightArrow,mps[site],S,U,V,0);

  // save matrix[site]
  save_site(mps, site, filename);
  mps[site].clear();
  cout << S << endl;
  return std::make_tuple(S, U.qshape()[0]);
}




