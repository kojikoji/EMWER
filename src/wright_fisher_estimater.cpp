#include<iostream>
#include<q_function.h>
#include<wright_fisher_estimater.h>
#include<lbfgsb.hpp>
#include<data.h>
#include<math.h>
#include <chrono>
bool WFE::refresh(Vec_st ntheta){
  qfun = 0;
  //qfun_d
  qfun_dpop = 0;
  qfun_dslc = 0;
  qfun_ddom = 0;
  bool fail = cqf.refresh(ntheta);
  qfun += cqf.llh();
  //llh_d
  qfun_dpop += cqf.llh_dpop();
  qfun_dslc += cqf.llh_dslc();
  qfun_ddom += cqf.llh_ddom();
  return(fail);
}
WFE::CMR::CMR(WFE &_wfe,double c,vector<bool> opt_flag):upper(_wfe),_c(c),_i(0),_opt_flag(opt_flag){}
int WFE::CMR::operator()(vector<double>& x,double& fn,vector<double>& gr){
  //unnormalized prameter(n,s) theta is caliculated
  //x is changed
  Vec_3d ptheta = upper.theta;
  Vec_3d ntheta = ptheta;
  for(int i = 0; i < x.size(); i++){
    ntheta(i) = x[i];
  }
  upper.refresh(ntheta);
  //compute qfun and qfun_d
  //subsititute for fn and gr
  fn = -upper.qfun;
  gr.assign(x.size(),0);
  if(_opt_flag[0]){
    gr[0] += -upper.qfun_dpop;
  }
  if(_opt_flag[1]){
    gr[1] += -upper.qfun_dslc;
  }
  if(_opt_flag[2]){
    gr[2] += -upper.qfun_ddom;
  }
  return(0);
}
template<class T>
void WFE::steepest_descent(T clh,vector<double> sbound,vector<double> dbound){//compute likelihood
  double pslc = theta(1);
  vector<double> lower = {10,sbound[0],dbound[0]};
  vector<double> upper = {3000,sbound[1],dbound[1]};
  vector<int> nbd = {2,2,2};
  Lbfgsb minimizer;
  minimizer.set_eps(1.0e-4);
  minimizer.set_maxit(30);
  minimizer.set_bounds(lower,upper,nbd);
  vector<double> ntheta0;
  for(int i = 0; i < theta.size(); i++){
    ntheta0.push_back(theta(i));
  }
  minimizer.minimize(ntheta0, clh);
  vector<double> ntheta = minimizer.best_x();
  for(int i = 0; i < ntheta.size(); i++){
    theta(i) = ntheta[i];
  }
}
WFE::WFE(){}
void WFE::load_data(Data& data){
  cqf.load_data(data);
}
void WFE::data_clear(){
  cqf.data_clear();
}
void WFE::optimize(Data& data,Vec_3d ftheta,vector<bool> _opt_flag,double _genpt,int _snum,vector<double> afs,BDR bdr, double beta, double dom_beta,bool avd_flag){
  sign_slc = 1;
  opt_flag = _opt_flag;
  genpt = _genpt;
  snum = _snum;
  //setting for cqf
  double cqf_num = 1;
  lh = 0;
  cqf.init(genpt,snum,afs,beta, dom_beta);
  cqf.load_data(data);
  cqf.afs_opt = opt_flag[3];
  pz_vec = cqf.get_pz();
  theta = ftheta;
  bool calc_fail = cqf.arefresh(theta);
  lh = cqf.get_log_pe();
  init_lh = lh;
  init_snum = snum;
  if(theta(0) == 0){
    theta(0) = 100;
  }
  double est_upper = 40;
  if(opt_flag[0]==true){
    est_upper = 1000;
  }
  int count = 0;
  //cout<<lh<<endl;
  bool end_flag;
  int safe_count = 0;
  double plh;
  Vec_3d ptheta;
  double  sdlt = 1.0e-3*snum;
  do{
    //init snum changer
    bool snum_change = false;
    //refresh parameter
    ptheta = theta;
    plh = lh;
    //snum check
    bool fail = refresh(ptheta);
    if(opt_flag[1]){
      if(fail){
	snum_change = true;
      }
      sdlt = 0.1;
      int sign_ptheta = (qfun_dslc > 0) - (qfun_dslc < 0);
      double qfun_dslc_ptheta = qfun_dslc;
      // when selection become negative, change the allele for caluculation
      if(qfun_dslc < 0 && count==0){
	cqf.change_allele();
	cqf.arefresh(theta);
	fail = refresh(ptheta);
	sign_slc = -1;
	count++;
	continue;
      }
      double prop = 1.0;
      Vec_st sb_theta;
      int search_count = 0;
      // search the boundary for the direction indicated by gradients
      sb_theta = ptheta;
      if(qfun_dslc > 0){
	sb_theta(1) = bdr.get_upper();
      }
      if(qfun_dslc < 0){
	sb_theta(1) = bdr.get_lower();
      }
      // increase state number
      // if gradient direction is same at the boundary for discretization
      refresh(sb_theta);
      int sign_bound = (qfun_dslc > 0) - (qfun_dslc < 0);
      if(sign_bound == sign_ptheta){
	snum_change = true;
	calc_fail = cqf.arefresh(theta);
      }
      //snum change
      if(snum_change){
	bool fail_increase = bdr.increase_state();
	// stop estimation when state number reaches upper limit
	if(fail_increase){
	  theta = sb_theta;
	  calc_fail = cqf.arefresh(theta);
	  lh = cqf.get_log_pe();
	  break;
	}
	snum = bdr.get_snum();
	cqf.init(genpt,snum,afs,beta, dom_beta);
	pz_vec = cqf.get_pz();
	calc_fail = cqf.arefresh(theta);
	lh = cqf.get_log_pe();
	fail = false;
	continue;
      }
    }
    vector<double> sbound = {bdr.get_lower(),bdr.get_upper()};
    vector<double> dbound = {bdr.get_dlower(),bdr.get_dupper()};
    //end bound setting
    //steepest decent
    double c = 2.0;
    CMR cmr((*this),c,opt_flag);
    steepest_descent(cmr,sbound,dbound);
    count++;
    calc_fail = cqf.arefresh(theta);
    lh = cqf.get_log_pe();
    Vec_3d res_vec = theta - ptheta;
    end_flag = fabs(res_vec(0)) < 1.0e-2 && fabs(res_vec(1)) < 1.0e-3 && fabs(res_vec(2)) < 1.0e-2;
    if(opt_flag[3]){
      double error = (cqf.get_pz() - pz_vec).sum();
      end_flag = end_flag && (error < 1.0e-3);
      pz_vec = cqf.get_pz();
    }
    end_flag = end_flag || std::isnan(lh) || lh < plh;
    if(count >= est_upper){
      cout<<"stop est"<<endl;
      calc_fail = true;
    }
  }while(!end_flag && count < est_upper && !calc_fail);
  if(calc_fail){
    cout<<"calc failure"<<endl;
    sign_slc = 0;
  }
}
Vec_3d WFE::ans_get(){
  return(theta);
}
Vec_3d WFE::afs_get(){
  return(pz_vec);
}
// -2log(likelihood ratio)
double WFE::stat_chi2(){
  double log_pe = lh;
  double log_pe0 = init_lh;
  if(snum != init_snum){
    Vec_3d slc0 = theta;
    slc0(1) = 0;
    cqf.arefresh(slc0);
    log_pe0 = cqf.get_log_pe();
  }
  double ans = -2*(log_pe0-log_pe);
  return(ans);
}
// variance of estimation
Vec_3d WFE::estimate_variance(){
  Vec_3d est_var = Vec_3d::Zero(3,1);
  double oq_dpop = cqf.llh_dpop();
  double oq_dslc = cqf.llh_dslc();
  double oq_ddom = cqf.llh_ddom();
  double delta = 1.0e-5;
  //theta change  slc
  Vec_3d theta_slc = theta;
  theta_slc(1) += delta;
  cqf.arefresh(theta_slc);
  double ldss;
  //dslcslc
  double q_dslc = cqf.llh_dslc();
  double num_dslcslc = (1/delta)*(q_dslc - oq_dslc);
  if(opt_flag[1] && opt_flag[2]){
    //dslcdom
    double q_ddom = cqf.llh_ddom();
    double num_dslcdom = (1/delta)*(q_ddom - oq_ddom);
    //theta change dom
    Vec_3d theta_dom = theta;
    theta_dom(2) += delta;
    cqf.arefresh(theta_dom);
    //ddomdom
    q_ddom = cqf.llh_ddom();
    double num_ddomdom = (1/delta)*(q_ddom - oq_ddom);
    Eigen::Matrix2d info_mat;
    info_mat <<num_dslcslc,num_dslcdom,
      num_dslcdom,num_ddomdom;
    Eigen::Matrix2d info_mat_inv = info_mat.inverse();
    est_var(0,0) = -info_mat_inv(0,0);
    est_var(1,0) = -info_mat_inv(1,1);
  }else{
    est_var(0,0) = -1/num_dslcslc;
    if(!opt_flag[1] && opt_flag[0]){
      //theta change  pop
      Vec_3d theta_pop = theta;
      theta_pop(0) += delta;
      cqf.arefresh(theta_pop);
      //dpoppop
      double q_dpop = cqf.llh_dpop();
      double num_dpoppop = (1/delta)*(q_dpop - oq_dpop);
      est_var(0,0) = -1/num_dpoppop;
    }      
  }
  est_var(2,0) = sign_slc;
  return(est_var);
}
