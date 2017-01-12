#include<Matrix.h>
#include<transition_matrix.h>
CTM::CTM(){}
void CTM::init(int snum){
  STNUM = snum;
  MSTNUM = snum -1;
  H = (double)1.0/MSTNUM;
  np_vec.resize(MSTNUM,1);
  nm_vec.resize(MSTNUM,1);
  tf_vec.resize(STNUM,1);
  sym_c.resize(STNUM,1);
  sym_s.resize(MSTNUM,1);
  dg.resize(STNUM,1);
  sqdg.resize(STNUM,1);
  isqdg.resize(STNUM,1);
  ld.resize(STNUM,1);
  r_c.resize(STNUM,1);
  r_p.resize(MSTNUM,1);
  r_m.resize(MSTNUM,1);
  k_uv.resize(STNUM,STNUM);
  u_mat.resize(STNUM,STNUM);
  ui_mat.resize(STNUM,STNUM);
}
double CTM::f_i(int i){
  double ans = H*i*(1-H*i);
  return(ans);
}
double CTM::f_p(Vec_2d theta){
  double ans = pow(theta(0),2)*exp(theta(1)*H)/(2*H*H);
  return(ans);
}
double CTM::f_m(Vec_2d theta){
  double ans = pow(theta(0),2)*exp(-theta(1)*H)/(2*H*H);
  return(ans);
}
double CTM::fc_dsig(Vec_2d theta){
  double ans;
  ans = -theta(0)*(exp(theta(1)*H) + exp(-theta(1)*H))/(H*H);
  return(ans);
}
double CTM::fc_dgam(Vec_2d theta){
  double ans;
  ans = pow(theta(0),2)*(exp(-theta(1)*H) - exp(theta(1)*H))/(2*H);
  return(ans);
}
double CTM::fp_dsig(Vec_2d theta){
  double ans;
  ans = 2.0/theta(0);
  return(ans);
}
double CTM::fp_dgam(Vec_2d theta){
  double ans;
  ans = H;
  return(ans);
}
double CTM::fm_dsig(Vec_2d theta){
  double ans;
  ans = 2.0/theta(0);
  return(ans);
}
double CTM::fm_dgam(Vec_2d theta){
  double ans;
  ans = -H;
  return(ans);
}
//ld and u and ui mat refreshed about inner
void CTM::eigen_refresh(){
  //eigen solver setting
  
    //eigen solver setting
  char jobz[] = "V";
  int ISTNUM = STNUM-2;
  int nrow = ISTNUM;
  Mat_st p_mat = Mat_st::Zero(ISTNUM,ISTNUM);
  int ldz = ISTNUM;
  double work[2*ISTNUM-2];
  //cout <<"sym_s"<<endl<<sym_s.str()<<endl;
  //cout <<"sym_c"<<endl<<sym_c.str()<<endl;
  Vec_st isym_c = Vec_st::Zero(ISTNUM,1);
  Vec_st isym_s = Vec_st::Zero(ISTNUM-1,1);
  for(int i = 0;i < ISTNUM;i++){
    isym_c(i) = sym_c(i+1);
  }
  for(int i = 0;i < ISTNUM-1;i++){
    isym_s(i) = sym_s(i+1);
  }
  int info;
  Vec_st ild = isym_c;
  //solve eigen
  dstev_(jobz,&nrow,ild.get(),isym_s.get(),p_mat.get(),&ldz,work,&info);
  if(info != 0){
    std::cout <<"Error!:eigen refresh failed"<<endl;
  }
  ld = Vec_st::Zero(STNUM,1);
  u_mat = Mat_st::Zero(STNUM,STNUM);
  for(int i = 0; i < ISTNUM;i++){
    ld(i) = ild(i);
    for(int j = 0;j < ISTNUM;j++){
      u_mat(i+1,j) = p_mat(i,j);
    }
  }
  //cout<<u_mat.str();
  ui_mat = u_mat.transpose();
}
//refresh pi_, mi_ and ci_ mat
void  CTM::i_mat_refresh(){
  //make modified uimat for j and umat for i
  pi_mat = Mat_st::Zero(STNUM,STNUM);
  mi_mat = Mat_st::Zero(STNUM,STNUM);
  ci_mat = Mat_st::Zero(STNUM,STNUM);
  Mat_st idu_mat = diag(isqdg)*u_mat;
  Mat_st uid_mat = ui_mat*diag(sqdg);
  ci_mat = uid_mat*diag(fiv)*idu_mat;
  pi_mat = uid_mat.col_cut(1,STNUM-1)*diag(r_p)*idu_mat.row_cut(0,MSTNUM-1);
  mi_mat = uid_mat.col_cut(0,MSTNUM-1)*diag(r_m)*idu_mat.row_cut(1,STNUM-1);
}

  //refresh u_mat, ui_mat, ld and dg
void  CTM::u_mat_refresh(){
  //change dg for u_mat
  dg(0) = dg(1);
  dg(STNUM-1) = dg(MSTNUM-1);
  for(int i = 0;i < dg.size();i++){
    sqdg(i) = sqrt(dg(i));
    isqdg(i) = sqrt(1/dg(i));      
  }
  bool flag = true;
  for(int j = 0;j < STNUM;j++){
    if(ld(j) == 0){
      if(flag){
	u_mat(0,j) = 1;
	flag = false;
      }else{
	u_mat(STNUM-1,j) = 1;
      }
    }else{
      u_mat(0,j) = u_mat(1,j)*r_m(0)/ld(j);
      u_mat(STNUM-1,j) = u_mat(STNUM-2,j)*r_p(MSTNUM-1)/ld(j);
    }
  }
  //cout<<"1ust-1*ui1"<<endl<<(u_mat.row(STNUM-1)*ui_mat.col(1)).str()<<endl;
  for(int j = 1;j < STNUM-1;j++){
    ui_mat(STNUM-2,j) = -u_mat.row(0)*ui_mat.col(j);
    ui_mat(STNUM-1,j) = -u_mat.row(STNUM-1)*ui_mat.col(j);
  }
  ui_mat(STNUM-2,0) = 1;
  ui_mat(STNUM-1,STNUM-1) = 1;
}
//refresh epmx
void  CTM::epmx_refresh(){  
  Mat_st ermt = u_mat*ui_mat;
  for(int i = 0;i < STNUM;i++){
    for(int j = 0;j < STNUM;j++){
      if(i != j&&epmx < fabs(ermt(i,j))){
	epmx = ermt(i,j);
      }
    }
  }
  //epmx = 1.0e-15;
  ldmx = -10e100;
  for(int i = 0;i < STNUM;i++){
    ldmx = max(ldmx,ld(i));
  }
}

//all element refresh (including umat uimat)
void CTM::arefresh(Vec_2d theta){
  sym_s = Vec_mst::Zero(MSTNUM,1);
  sym_c = Vec_st::Zero(STNUM,1);
  dg = Vec_st::Zero(STNUM,1);
  fiv = Vec_st::Zero(STNUM,1);
  double a = 2*theta(1)*H;
  for(int i = 1;i < STNUM-1;i++){
    dg(i) = f_i(i)*exp(-a*i)*pow(theta(0),2)/(2*H*H);
    sym_c(i) = -dg(i)*(exp(a*(i+1.0/2)) + exp(a*(i-1.0/2)));
    sym_s(i-1) = sqrt(dg(i))*sqrt(dg(i-1))*exp(a*(i-1.0/2));
    r_p(i) = f_p(theta)*f_i(i);
    r_m(i-1) = f_m(theta)*f_i(i);
    r_c(i) = - r_p(i) - r_m(i-1);
    fiv(i) = f_i(i);
  }
  //cout<<"s"<<endl<<sym_s.str()<<endl;
  //cout<<"c"<<endl<<sym_c.str()<<endl;
  //cout<<"dg"<<endl<<dg.str()<<endl;
  r_p(0) = 0;
  r_m(MSTNUM-1) = 0;
  r_c(0) = 0;
  r_c(STNUM-1) = 0;
  eigen_refresh();
  u_mat_refresh();
  i_mat_refresh();
  epmx_refresh();
}
