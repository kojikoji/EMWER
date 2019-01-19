#include<transition_matrix.h>

template<class T>
CTM<T>::CTM(){}

template<class T>
void CTM<T>::init(int snum){
  STNUM = snum;
  MSTNUM = snum -1;
  ISTNUM = snum -2;
  dlt = (double)1.0/STNUM;
  ld.resize(STNUM,1);
  u_mat.resize(STNUM,STNUM);
  ui_mat.resize(STNUM,STNUM);
  //cout<<"substi"<<endl;
  //c mean center
  rc_vec = Vec_st::Zero(1,STNUM);
  //p mean plus
  rp_vec = Vec_st::Zero(1,STNUM);
  //m mean minus
  rm_vec = Vec_st::Zero(1,STNUM);
  //p mean plus
  lrp_vec = Vec_st::Zero(1,STNUM);
  //m mean minus
  lrm_vec = Vec_st::Zero(1,STNUM);
  //c mean center
  rcdneff_vec = Vec_st::Zero(1,STNUM);
  //p mean plus
  lrpdneff_vec = Vec_st::Zero(1,STNUM);
  //m mean minus
  lrmdneff_vec = Vec_st::Zero(1,STNUM);
  //c mean center
  rcdslc_vec = Vec_st::Zero(1,STNUM);
  //p mean plus
  lrpdslc_vec = Vec_st::Zero(1,STNUM);
  //m mean minus
  lrmdslc_vec = Vec_st::Zero(1,STNUM);
  //c mean center
  rcddom_vec = Vec_st::Zero(1,STNUM);
  //p mean plus
  lrpddom_vec = Vec_st::Zero(1,STNUM);
  //m mean minus
  lrmddom_vec = Vec_st::Zero(1,STNUM);
  //cout<<"comp"<<endl;
}
template<class T>
double CTM<T>::f_rq(double x){
  double ans = genpt*x*(1-x)/(4*neff);
  return(ans);
};
template<class T>
double CTM<T>::f_rq_dneff(double x){
  double ans = -genpt*x*(1-x)/(4*neff*neff);
  return(ans);
};
template<class T>
double CTM<T>::f_rq_dslc(double x){
  double ans = 0;
  return(ans);
};
template<class T>
double CTM<T>::f_rq_ddom(double x){
  double ans = 0;
  return(ans);
};
template<class T>
double CTM<T>::f_rq_dslcslc(double x){
  double ans = 0;
  return(ans);
};
template<class T>
double CTM<T>::f_rq_dslcdom(double x){
  double ans = 0;
  return(ans);
};
template<class T>
double CTM<T>::f_rq_ddomdom(double x){
  double ans = 0;
  return(ans);
};
template<class T>
double CTM<T>::f_rqp(double x){
  double a = (1-2*x)/(4*neff);
  double b = slc*x*(1-x)*(x+dom-2*dom*x)/(1+slc*x*(x+2*dom-2*dom*x));
  return(genpt*(a-b));
};
template<class T>
double CTM<T>::f_rqp_dneff(double x){
  double ans = -genpt*(1-2*x)/(4*neff*neff);
  return(ans);
};
template<class T>
double CTM<T>::f_rqp_dslc(double x){
  double bunbo = (1+slc*x*(x+2*dom-2*dom*x));
  double a = -x*(1-x)*(x+dom-2*dom*x)/bunbo;
  double b = slc*x*(1-x)*(x+dom-2*dom*x)*x*(x+2*dom-2*dom*x)/(bunbo*bunbo);
  return(genpt*(a+b));
};
template<class T>
double CTM<T>::f_rqp_ddom(double x){
  double bunbo = (1+slc*x*(x+2*dom-2*dom*x));
  double a = -slc*x*(1-x)*(1-2*x)/bunbo;
  double b = 2*slc*slc*x*(1-x)*(x+dom-2*dom*x)*x*(1-x)/(bunbo*bunbo);
  return(genpt*(a+b));
};
template<class T>
double CTM<T>::f_rqp_dslcslc(double x){
  double b = x*(x+ 2*dom*(1-x));
  double c = x*(1-x)*(x + dom*(1-2*x));
  double var1 = -2*slc*c*b*b*pow(1+slc*b,-3);
  double var2 = 2*c*b*pow(1+slc*b,-2);
  return(var1+var2);
};
template<class T>
double CTM<T>::f_rqp_dslcdom(double x){
  double b = x*(x+ 2*dom*(1-x));
  double bp = 2*x*(1-x);
  double c = x*(1-x)*(x + dom*(1-2*x));
  double cp = x*(1-x)*(1-2*x);
  double var1 = -2*pow(slc,2)*c*b*bp*pow(1+slc*b,-3);  
  double var2 = slc*cp*b*pow(1+slc*b,-2);  
  double var3 = 2*slc*c*bp*pow(1+slc*b,-2);  
  double var4 = -cp*pow(1+slc*b,-1);  
  return(var1+var2+var3+var4);
};
template<class T>
double CTM<T>::f_rqp_ddomdom(double x){
  double b = x*(x+ 2*dom*(1-x));
  double bp = 2*x*(1-x);
  double c = x*(1-x)*(x + dom*(1-2*x));
  double cp = x*(1-x)*(1-2*x);
  double var1 = -2*pow(slc,3)*c*pow(bp,2)*pow(1+slc*b,-3);  
  double var2 = 2*pow(slc,2)*cp*bp*pow(1+slc*b,-2);  
  return(var1+var2);
};
template<class T>
double CTM<T>::f_r(double x){
  double ans = exp(-2*neff*log(1+slc*x*(x+2*dom*(1-x))));
 return(ans);
}
template<class T>
double CTM<T>::f_irdneff(double x){
  double inlog = (1+slc*x*(x+2*dom*(1-x)));
  double inlogdneff = -2*log(inlog);
  double ans = inlogdneff*f_r(x);
  ans = -ans*pow(f_r(x),-2);
  return(ans);
}
template<class T>
double CTM<T>::f_irdslc(double x){
  double inlog = (1+slc*x*(x+2*dom*(1-x)));
  double inlogddom = x*(x + 2*dom*(1-x));
  double ans = -2*neff*(inlogddom/inlog)*f_r(x);
  ans = -ans*pow(f_r(x),-2);
  return(ans);
}
template<class T>
double CTM<T>::f_irddom(double x){
  double inlog = (1+slc*x*(x+2*dom*(1-x)));
  double inlogddom = 2*slc*x*(1-x);
  double ans = -2*neff*(inlogddom/inlog)*f_r(x);
  ans = -ans*pow(f_r(x),-2);
  return(ans);
}
template<class T>
double CTM<T>::f_q(double x){
  double coff = 4*neff/(genpt*x*(1-x));
  double ans = coff/f_r(x);
  return(ans);
}
template<class T>
double CTM<T>::f_iqdneff(double x){
  double coff = 4/(genpt*x*(1-x));
  double ans = coff/f_r(x) + coff*neff*f_irdneff(x);
  ans = -ans*pow(f_q(x),-2);
  return(ans);
}
template<class T>
double CTM<T>::f_iqdslc(double x){
  double coff = 4*neff/(genpt*x*(1-x));
  double ans = coff*f_irdslc(x);
  ans = -ans*pow(f_q(x),-2);
  return(ans);
}
template<class T>
double CTM<T>::f_iqddom(double x){
  double coff = 4*neff/(genpt*x*(1-x));
  double ans = coff*f_irddom(x);
  ans = -ans*pow(f_q(x),-2);
  return(ans);
}

template<class T>
bool CTM<T>::param_refresh(Vec_3d theta){
  neff = theta(0);
  slc = theta(1);
  dom = theta(2);
  //original r mat
  int rp_flag= 0;  
  int rm_flag= 0;
  bool fail = false;
  //pointers
  double *rp_ptr = rp_vec.data();
  double *rm_ptr = rm_vec.data();
  double *lrp_ptr = lrp_vec.data();
  double *lrm_ptr = lrm_vec.data();
  double *rc_ptr = rc_vec.data();
  for(int j = 1;j < STNUM-1;j++){
    double xp = dlt*(j+1);
    double xm = dlt*j;
    double dltj = dlt;
    double dltp = dlt;
    double dltm = dlt;
    rp_ptr[j] = (2*f_rq(xp) - dltp*f_rqp(xp))/(dltj*(dltj+dltp));
    rm_ptr[j] = (2*f_rq(xm) + dltm*f_rqp(xm))/(dltj*(dltj+dltm));
    lrm_ptr[j] = log(rm_ptr[j]);
    lrp_ptr[j] = log(rp_ptr[j]);
    rc_ptr[j] = -rp_ptr[j] - rm_ptr[j];
    if(rp_ptr[j]<0){
      rp_flag--;
    }
    if(rm_ptr[j]<0){
      rm_flag--;
    }
  }
  if(rp_flag < 0){
    fail = true;
  }
  if(rm_flag < 0){
    fail = true;
  }
  //dneff r mat
  //pointer
  double *lrpdneff_ptr = lrpdneff_vec.data();
  double *lrmdneff_ptr = lrmdneff_vec.data();
  double *rcdneff_ptr = rcdneff_vec.data();
  //m mean minus
  for(int j = 1;j < STNUM-1;j++){
    double xp = dlt*(j+1);
    double xm = dlt*j;
    double dltj = dlt;
    double dltp = dlt;
    double dltm = dlt;
    double rmdneff = (2*f_rq_dneff(xm) + dltp*f_rqp_dneff(xm))/(dltj*(dltj+dltm));
    double rpdneff = (2*f_rq_dneff(xp) - dltp*f_rqp_dneff(xp))/(dltj*(dltj+dltp));
    lrmdneff_ptr[j] = rmdneff/rm_vec(j);
    lrpdneff_ptr[j] = rpdneff/rp_vec(j);    
    rcdneff_ptr[j] = -rmdneff -rpdneff;
  }

  //dslc r mat
  //pointer
  double *lrpdslc_ptr = lrpdslc_vec.data();
  double *lrmdslc_ptr = lrmdslc_vec.data();
  double *rcdslc_ptr = rcdslc_vec.data();
  for(int j = 1;j < STNUM-1;j++){
    double xp = dlt*(j+1);
    double xm = dlt*j;
    double dltj = dlt;
    double dltp = dlt;
    double dltm = dlt;
    double rmdslc = (2*f_rq_dslc(xm) + dltp*f_rqp_dslc(xm))/(dltj*(dltj+dltm));
    double rpdslc = (2*f_rq_dslc(xp) - dltp*f_rqp_dslc(xp))/(dltj*(dltj+dltp));
    lrmdslc_ptr[j] = rmdslc/rm_ptr[j];
    lrpdslc_ptr[j] = rpdslc/rp_ptr[j];    
    rcdslc_ptr[j] = -rmdslc -rpdslc;
  }

  //ddom r mat
  double *lrpddom_ptr = lrpddom_vec.data();
  double *lrmddom_ptr = lrmddom_vec.data();
  double *rcddom_ptr = rcddom_vec.data();
  for(int j = 1;j < STNUM-1;j++){
    double xp = dlt*(j+1);
    double xm = dlt*j;
    double dltj = dlt;
    double dltp = dlt;
    double dltm = dlt;
    double rmddom = (2*f_rq_ddom(xm) + dltp*f_rqp_ddom(xm))/(dltj*(dltj+dltm));
    double rpddom = (2*f_rq_ddom(xp) - dltp*f_rqp_ddom(xp))/(dltj*(dltj+dltp));
    lrmddom_ptr[j] = rmddom/rm_ptr[j];
    lrpddom_ptr[j] = rpddom/rp_ptr[j];    
    rcddom_ptr[j] = -rmddom -rpddom;
  }
  return(fail);
}
template<class T>
bool CTM<T>::dd_refresh(Vec_3d theta){
  bool fail = param_refresh(theta);
  //slcslc size set
  //c mean center
  rcdslcslc_vec = Vec_st::Zero(1,STNUM);
  //p mean plus
  lrpdslcslc_vec = Vec_st::Zero(1,STNUM);
  //m mean minus
  lrmdslcslc_vec = Vec_st::Zero(1,STNUM);
  //domdom size set
  //c mean center
  rcddomdom_vec = Vec_st::Zero(1,STNUM);
  //p mean plus
  lrpddomdom_vec = Vec_st::Zero(1,STNUM);
  //m mean minus
  lrmddomdom_vec = Vec_st::Zero(1,STNUM);
  //slcdom size set
  //c mean center
  rcdslcdom_vec = Vec_st::Zero(1,STNUM);
  //p mean plus
  lrpdslcdom_vec = Vec_st::Zero(1,STNUM);
  //m mean minus
  lrmdslcdom_vec = Vec_st::Zero(1,STNUM);

  //dslcslc
  for(int j = 1;j < STNUM-1;j++){
    double xp = dlt*(j+1);
    double xm = dlt*j;
    double dltj = dlt;
    double dltp = dlt;
    double dltm = dlt;
    double rmdslc = (2*f_rq_dslc(xm) + dltp*f_rqp_dslc(xm))/(dltj*(dltj+dltm));
    double rpdslc = (2*f_rq_dslc(xp) - dltp*f_rqp_dslc(xp))/(dltj*(dltj+dltp));
    double rmdslcslc = (2*f_rq_dslcslc(xm) + dltp*f_rqp_dslcslc(xm))/(dltj*(dltj+dltm));
    double rpdslcslc = (2*f_rq_dslcslc(xp) - dltp*f_rqp_dslcslc(xp))/(dltj*(dltj+dltp));
    lrmdslcslc_vec(j) = -pow(rmdslc/rm_vec(j),2) + rmdslcslc/rm_vec(j);
    lrpdslcslc_vec(j) = -pow(rpdslc/rp_vec(j),2) + rpdslcslc/rp_vec(j);
    rcdslcslc_vec(j) = -rmdslcslc -rpdslcslc;
  }

  //ddomdom
  for(int j = 1;j < STNUM-1;j++){
    double xp = dlt*(j+1);
    double xm = dlt*j;
    double dltj = dlt;
    double dltp = dlt;
    double dltm = dlt;
    double rmddom = (2*f_rq_ddom(xm) + dltp*f_rqp_ddom(xm))/(dltj*(dltj+dltm));
    double rpddom = (2*f_rq_ddom(xp) - dltp*f_rqp_ddom(xp))/(dltj*(dltj+dltp));
    double rmddomdom = (2*f_rq_ddomdom(xm) + dltp*f_rqp_ddomdom(xm))/(dltj*(dltj+dltm));
    double rpddomdom = (2*f_rq_ddomdom(xp) - dltp*f_rqp_ddomdom(xp))/(dltj*(dltj+dltp));
    lrmddomdom_vec(j) = -pow(rmddom/rm_vec(j),2) + rmddomdom/rm_vec(j);
    lrpddomdom_vec(j) = -pow(rpddom/rp_vec(j),2) + rpddomdom/rp_vec(j);
    rcddomdom_vec(j) = -rmddomdom -rpddomdom;
  }
  //ddomslc
  for(int j = 1;j < STNUM-1;j++){
    double xp = dlt*(j+1);
    double xm = dlt*j;
    double dltj = dlt;
    double dltp = dlt;
    double dltm = dlt;
    double rmdslc = (2*f_rq_dslc(xm) + dltp*f_rqp_dslc(xm))/(dltj*(dltj+dltm));
    double rpdslc = (2*f_rq_dslc(xp) - dltp*f_rqp_dslc(xp))/(dltj*(dltj+dltp));
    double rmddom = (2*f_rq_ddom(xm) + dltp*f_rqp_ddom(xm))/(dltj*(dltj+dltm));
    double rpddom = (2*f_rq_ddom(xp) - dltp*f_rqp_ddom(xp))/(dltj*(dltj+dltp));
    double rmdslcdom = (2*f_rq_dslcdom(xm) + dltp*f_rqp_dslcdom(xm))/(dltj*(dltj+dltm));
    double rpdslcdom = (2*f_rq_dslcdom(xp) - dltp*f_rqp_dslcdom(xp))/(dltj*(dltj+dltp));
    lrmdslcdom_vec(j) = -rmdslc*rmddom*pow(rm_vec(j),-2) + rmdslcdom/rm_vec(j);
    lrpdslcdom_vec(j) = -rpdslc*rpddom*pow(rp_vec(j),-2) + rpdslcdom/rp_vec(j);
    rcdslcdom_vec(j) = -rmdslcdom -rpdslcdom;
  }

  return(fail);
}
/* DGEEV prototype */
extern void dgeev( char* jobvl, char* jobvr, int* n, double* a,
                int* lda, double* wr, double* wi, double* vl, int* ldvl,
                double* vr, int* ldvr, double* work, int* lwork, int* info );

//ld and u and ui mat refreshed about inner
template<class T>
void CTM<T>::eigen_refresh(){
  Mat_st r_mat = make_r_mat();
  Eigen::EigenSolver<Eigen::MatrixXd> es(r_mat);
  ld = es.eigenvalues().real();
  u_mat = (es.eigenvectors()).real();
  ui_mat = u_mat.inverse();    
}
template<>
void CTM<std::complex<double>>::eigen_refresh(){
  Mat_st r_mat = make_r_mat();
  Eigen::EigenSolver<Eigen::MatrixXd> es(r_mat);
  ld = es.eigenvalues();
  u_mat = (es.eigenvectors());
  ui_mat = u_mat.inverse();    
}
//all element refresh (including umat uimat)
template<class T>
bool CTM<T>::arefresh(Vec_3d theta){
  neff = theta(0);
  slc = theta(1);
  dom = theta(2);
  bool fail = param_refresh(theta);
  eigen_refresh();
  kap_mat = Mat_gen<T>::Zero(STNUM,STNUM);
  return(fail);
}

template<class T>
Mat_st CTM<T>::make_r_mat(){
  Mat_st r_mat = Mat_st::Zero(STNUM,STNUM);
  for(int j = 1;j < STNUM-1;j++){
    r_mat(j-1,j) = rm_vec(j);
    r_mat(j+1,j) = rp_vec(j);
    r_mat(j,j) = -r_mat(j+1,j) - r_mat(j-1,j);
  }
  return(r_mat);
}
      
template<class T>
Mat_st CTM<T>::make_u_ldn_ui(double n){
  Mat_st ldn_mat = Mat_st::Zero(STNUM,STNUM);
  for(int j = 0;j < STNUM;j++){
    ldn_mat(j,j) = pow(ld(j),n);
  }
  return(u_mat*ldn_mat*ui_mat);
}

template<class T>
Mat_st CTM<T>::make_exp_tr(double t){
  Mat_st eld_mat = Mat_st::Zero(STNUM,STNUM);
  for(int j = 0;j < STNUM;j++){
    eld_mat(j,j) = exp(t*ld(j));
  }
  return(u_mat*eld_mat*ui_mat);
}

template<>
Mat_st CTM<std::complex<double>>::make_u_ldn_ui(double n){
  Mat_gen<std::complex<double>> ldn_mat = Mat_gen<std::complex<double>>::Zero(STNUM,STNUM);
  for(int j = 0;j < STNUM;j++){
    ldn_mat(j,j) = pow(ld(j),n);
  }
  return((u_mat*ldn_mat*ui_mat).real());
}

template<>
Mat_st CTM<std::complex<double>>::make_exp_tr(double t){
  Mat_gen<std::complex<double>> eld_mat = Mat_gen<std::complex<double>>::Zero(STNUM,STNUM);
  for(int j = 0;j < STNUM;j++){
    eld_mat(j,j) = exp(t*ld(j));
  }
  return((u_mat*eld_mat*ui_mat).real());
}


template<class T>
Mat_st CTM<T>::make_exp_tr_num(double t,int num){
  return(Mat_st::Zero(STNUM,STNUM));
}
template class CTM<double>;
template class CTM<std::complex<double>>;
