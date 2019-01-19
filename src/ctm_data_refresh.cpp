//#include<Matrix.h>
//using namespace std;
#include <iostream>
#include<transition_matrix.h>
#include<random>
#include <chrono>
#define KAPOPT
auto start3 = std::chrono::system_clock::now();
auto end3 = std::chrono::system_clock::now();
auto dur3 = end3 - start3;
auto msec3 = std::chrono::duration_cast<std::chrono::milliseconds>(dur3).count();
template <class T>
void print_vec(vector<T> vec){
  for( auto itr = vec.begin(); itr != vec.end(); itr++){
    cout<<"\t"<<*itr;
  }
  cout<<endl;
}
double no_minus(double val){
  if(val > 0){
    return(val);
  }else{
    return(0);
  }
}

template <class T1,class T2, class T3> 
T3 multi_each(T1 mat1,T2 mat2){
  T3 ans = T3::Zero(mat1.rows(),mat1.cols());
  for(int i = 0;i < mat1.rows();i++){
    for(int j = 0;j < mat1.cols();j++){
      ans(i,j) = mat1(i,j)*mat2(i,j);
    }
  }
  return(ans);
}

Mat_stc no_minus(Mat_stc mat){
  Mat_stc ans = Mat_stc::Zero(mat.rows(),mat.cols());
  for(int i = 0;i < mat.rows();i++){
    for(int j = 0;j < mat.cols();j++){
      if(mat(i,j).real() > 0){
	ans(i,j) = mat(i,j);
      }else{
	ans(i,j) = 0;
      }
    }
  }
  return(ans);
}
bool minus_check(Mat_st mat){
  bool ans = false;
  for(int i = 0;i < mat.rows();i++){
    for(int j = 0;j < mat.cols();j++){
      if(mat(i,j) < 0){
	ans = true;
      }
    }
  }
  return(ans);
}
template <class T>
Mat_gen<T> element_each(T el,int N,int M){
  Mat_gen<T> ans = Mat_gen<T>::Zero(N,M);
  for(int i = 0;i < N;i++){
    for(int j = 0;j < M;j++){
      ans(i,j) = el;
    }
  }
  return(ans);
}
double binom_ctm(double p,int a,int b){
  int n = a+b;
  //combination and p
  double ans = 1; 
  for(int i = 0;i < a;i++){
    ans *= ((double)(n-i)/(a-i))*p;
  }
  for(int i = 0;i < b;i++){
    ans *= 1-p;
  }  
  if(p == 1){
    if(b==0){
      ans = 1;
    }else{
      ans = 0;
    }
  }
  if(p == 0){
    if(a == 0){
      ans = 1;
    }else{
      ans = 0;
    }
  }
  return(ans);
}
template<class T>
Mat_gen<T> CTM<T>::diag(Mat_gen<T> vec){
  Mat_gen<T> ans = Mat_gen<T>::Zero(vec.rows(),vec.rows());
  for(int i = 0;i < vec.rows();i++){
    ans(i,i) = vec(i,0);
  }
  return(ans);
}
template <class T>
double CTM<T>::to_xi(int i){
  double h = (double)1/MSTNUM;
  return(h*i);
}
template <class T>
Mat_gen<T> CTM<T>::exp_eg(double t){
  Mat_gen<T> exp_eigens(STNUM,1);
  //time setting
  for(int i = 0;i < STNUM;i++){
    exp_eigens(i) = exp(t*ld(i));
  }
  return(exp_eigens);
}

//calucu etr
template<class T>
Mat_gen<T> CTM<T>::etr_make(double t){
  Mat_gen<T> exp_eigens = exp_eg(t);
  //Mat_st er_mat = element_each<double>(exp(t*ldmx)*epmx,STNUM,STNUM);
  Mat_gen<T> etq = diag(exp_eigens);
  //error correct from biased zooming dg
  //etq = no_minus(etq - er_mat);
  Mat_gen<T> etr = etq;
  return(etr);
}
// pre pab data

//caluculate data wise 
template<class T>
bool CTM<T>::lpe_refresh(vector<Vec_st> _pabl,vector<double> _tpl,Vec_st pz){
  tpl = _tpl;
  Vec_st tone_vec = element_each((double)1,1,STNUM);
  //Mat_st er_mat =  element_each<double>(epmx,STNUM,STNUM);
  //case
  dsize = tpl.size();
  pabx = _pabl;
  //calucurate p(alpha,beta|endpoint)
  //caluculate  etr
  vector<Mat_gen<T>> etr(dsize);
  for(int d = 0;d < dsize-1;d++){
    double t = tpl[d+1]-tpl[d];
    etr[d] = etr_make(t);

  }
  bool calc_fail = false;
  //forward
  fwd.resize(dsize);
  fwd[0] = multi_each<Mat_st,Mat_st,Mat_gen<T>>(pabx[0],pz);
  //cout<<fwd[0].str()<<endl;
  for(int d = 0;d < dsize-1;d++){
    Mat_gen<T> nfwd = u_mat*(etr[d]*(ui_mat*fwd[d]));
    //cout<<etr[d].str()<<endl;
    //cout<<pabx[d+1].str()<<endl;
    fwd[d+1] = multi_each<Mat_st,Mat_gen<T>,Mat_gen<T>>(pabx[d+1],nfwd);
    //if(minus_check(fwd[d+1].real())){
    /*
    if(fwd[d+1].real().minCoeff() < 0){
      //cout<<"negtive"<<endl;
      double prop = -fwd[d+1].real().minCoeff()/fwd[d+1].real().maxCoeff();
      if(prop > 1.0e-4){
	cout<<"negtive"<<endl;
	calc_fail = true;
      }
    }
    if(fwd[d+1].real().maxCoeff() > 1){
      //cout<<"over 1"<<endl;
      calc_fail = true;
    }
    */
  }
  //cout<<fwd[dsize-1].str()<<endl;
  //forward
  bwd.resize(dsize);
  bwd[dsize-1] = pabx[dsize-1];
  for(int d = 1;d < dsize;d++){
    int rd = dsize - d;
    Mat_gen<T> pbwdT = (((bwd[rd].transpose())*u_mat)*etr[rd-1])*ui_mat;
    bwd[rd-1] = multi_each<Mat_st,Mat_gen<T>,Mat_gen<T>>(pabx[rd-1],pbwdT.transpose());
    /*
    if(bwd[rd-1].real().minCoeff() < 0){
      //cout<<"back negtive"<<endl;
      double prop = -bwd[rd-1].real().minCoeff()/bwd[rd-1].real().maxCoeff();
      if(prop > 1.0e-4){
	cout<<"back negtive"<<endl;
	calc_fail = true;
      }
    }
    if(bwd[rd-1].real().maxCoeff() > 1){
      //cout<<"back over 1"<<endl;
      calc_fail = true;
    }
    */
    //cout<<bwd[rd-1].str()<<endl;
  }
  if(dsize == 1){
    pe = 1;
  }else{
    //get pe form onevec*fwd(D) into log
    pe = (element_each((double)1,1,STNUM)*fwd[dsize-1])(0,0);
  }
  //cout<<pe<<endl;
  return(calc_fail);

}

// to
template<class T>
void CTM<T>::kapd_refresh(){
  kapd_mat = Mat_gen<T>::Zero(STNUM,STNUM);
  for(int l = 0; l < dsize -1; l++){
    double t = tpl[l+1] - tpl[l];
    //cout<<t<<endl;
    Mat_gen<T> b_vec = (bwd[l+1].transpose()*u_mat).transpose();
    //cout<<b_vec.str()<<endl;
    Mat_gen<T> f_vec = ui_mat * fwd[l];
    Mat_gen<T> fb_mat = f_vec*b_vec.transpose();
    Mat_gen<T> tld = t*ld;
    Mat_gen<T> etld = tld;
    for(int u = 0; u < STNUM; u++){
      etld(u) = exp(tld(u));
    }
    T* tld_ptr = tld.data();
    T* etld_ptr = etld.data();
    T* kap_ptr = kapd_mat.data();
    T* fb_ptr = fb_mat.data();
    for(int v = 0; v < STNUM; v++){
      for(int u = 0; u < STNUM; u++){
	T tlu = tld_ptr[u];
	T tlv = tld_ptr[v];
	T val = (etld_ptr[u] - etld_ptr[v])/(tlu - tlv);
	if(abs(tlu - tlv) < 1.0e-50){
	  val = etld_ptr[u];
	}
	kap_ptr[u+v*STNUM] += t*fb_ptr[u+v*STNUM]*val/pe;
      }
    }
  }
  //cout<<kapd_mat.str()<<endl;
  //kapd_mat = kapd_mat.transpose();
  kap_mat = kap_mat + kapd_mat;
}
template<class T>
double CTM<T>::pe_get(){
  return(pe);
}
template<>
double CTM<std::complex<double>>::pe_get(){
  return(pe.real());
}

template<class T>
void CTM<T>::tkap_mat_refresh(){
  tkap_mat = kap_mat*ui_mat;
}
template<class T>
Vec_st CTM<T>::make_fs_vec(){
  Vec_st fs_vec = Vec_st::Zero(STNUM,1);
  for(int i =  1; i < STNUM-1; i++){
    fs_vec(i) = (u_mat.row(i)*tkap_mat.col(i))(0,0);
  }
  return(fs_vec);
}
    
template<class T>
Vec_st CTM<T>::make_np_vec(){
  Vec_st np_vec = Vec_st::Zero(STNUM,1);
  for(int i =  1; i < STNUM-1; i++){
    np_vec(i) = rp_vec(i)*(u_mat.row(i)*tkap_mat.col(i+1))(0,0);
  }
  return(np_vec);
}
template<class T>
Vec_st CTM<T>::make_nm_vec(){
  Vec_st nm_vec = Vec_st::Zero(STNUM,1);
  for(int i =  1; i < STNUM-1; i++){
    nm_vec(i) = (rm_vec(i)*(u_mat.row(i)*tkap_mat.col(i-1)))(0,0);
  }
  return(nm_vec);
}

template<>
Vec_st CTM<std::complex<double>>::make_fs_vec(){
  Vec_st fs_vec = Vec_st::Zero(STNUM,1);
  for(int i =  1; i < STNUM-1; i++){
    fs_vec(i) = (u_mat.row(i)*tkap_mat.col(i))(0,0).real();
  }
  return(fs_vec);
}
    
template<>
Vec_st CTM<std::complex<double>>::make_np_vec(){
  Vec_st np_vec = Vec_st::Zero(STNUM,1);
  for(int i =  1; i < STNUM-1; i++){
    np_vec(i) = rp_vec(i)*(u_mat.row(i)*tkap_mat.col(i+1))(0,0).real();
  }
  return(np_vec);
}
template<>
Vec_st CTM<std::complex<double>>::make_nm_vec(){
  Vec_st nm_vec = Vec_st::Zero(STNUM,1);
  for(int i =  1; i < STNUM-1; i++){
    nm_vec(i) = (rm_vec(i)*(u_mat.row(i)*tkap_mat.col(i-1)))(0,0).real();
  }
  return(nm_vec);
}

template<class T>
Vec_st CTM<T>::make_gam_vec(Vec_st &pz){
  Vec_st gam_vec = Vec_st::Zero(STNUM,1);
  for(int i =  0; i < STNUM; i++){
    gam_vec(i) = pz(i)*bwd[0](i)/pe;
    //gam_vec(i) = bwd[0](i)/bwd[0].sum();
    //gam_vec(i) = pabx[0](i);
  }
  return(gam_vec);
}
template<>
Vec_st CTM<std::complex<double>>::make_gam_vec(Vec_st &pz){
  Vec_st gam_vec = Vec_st::Zero(STNUM,1);
  for(int i =  0; i < STNUM; i++){
    gam_vec(i) = pz(i)*(bwd[0](i)).real()/pe.real();
    //gam_vec(i) = bwd[0](i)/bwd[0].sum();
    //gam_vec(i) = pabx[0](i);
  }
  return(gam_vec);
}

//明示的インスタンス化
//template class CTM<double>;
template class CTM<std::complex<double>>;
