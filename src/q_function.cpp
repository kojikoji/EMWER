#include<iostream>
#include<q_function.h>
#include <chrono>
//#define TIME(PROCESS,TAG)  start2 = chrono::system_clock::now();;PROCESS;end2 = std::chrono::system_clock::now();dur2 = end2 - start2;msec2 = std::chrono::duration_cast<std::chrono::microseconds>(dur2).count();std::cout<<TAG <<"\t"<< msec2 << " milli sec \n";
#ifndef TIME
#define TIME(PROCESS,TAG)  PROCESS;
#endif
auto start2 = chrono::system_clock::now();
auto end2 = std::chrono::system_clock::now();       // 計測終了時刻を保存
auto dur2 = end2 - start2;        // 要した時間を計算
auto msec2 = std::chrono::duration_cast<std::chrono::milliseconds>(dur2).count();
CQF::CQF():data(){
}
void CQF::data_clear(){
  data.clear();
}
Mat_st element_each_double(double el,int N,int M){
  Mat_st ans = Mat_st::Zero(N,M);
  for(int i = 0;i < N;i++){
    for(int j = 0;j < M;j++){
      ans(i,j) = el;
    }
  }
  return(ans);
}
void CQF::init(double _genpt,int snum,vector<double> afs, double _beta, double _dom_beta){
  STNUM = snum;
  MSTNUM = snum -1;
  H = (double)1.0/MSTNUM;
  ctmr.init(snum);
  genpt = _genpt;
  save_flag = false;
  pz = prior_make(afs);
  beta = _beta;
  dom_beta = _dom_beta;
  //cout<<dom_beta<<endl;
  data.pabx_refresh(snum);
}
void CQF::change_allele(){
  data.change_allele();
}
void CQF::load_data(Data& _data){
  data = _data;
  data.pabx_refresh(STNUM);
}
  //f~ and ctm.and ctm. is refreshed for R new pararmeter. 
bool CQF::refresh(Vec theta){
  slc = theta(1);
  dom = theta(2);
  bool fail = ctmr.param_refresh(theta);
  return(fail);
}
double get_mass(vector<double> afs,int ind, int length){
  double prop = (double)ind/(length-1);
  double dlt = 1.0/(length-1);
  double prop_max = prop + dlt/2;
  double prop_min = prop - dlt/2;
  int afs_size = afs.size();
  double afs_max = afs_size*prop_max;
  double afs_min = afs_size*prop_min;
  int afs_min_int = ceil(afs_min);
  int afs_max_int = floor(afs_max);
  double ans = 0;
  for(int i = std::max((double)afs_min_int,(double)0);i < afs_max_int;i++){
    ans += afs[i];
  }
  if(!(afs_max_int > afs_size-1)){
    ans += afs[afs_max_int]*(afs_max - std::max((double)afs_max_int,(double)afs_min_int));
  }
  if(!(afs_min_int < 0)){
    ans += afs[afs_min_int-1]*(afs_min_int - afs_min);
  }
  return(ans);
}
  
Vec_st CQF::prior_make(vector<double> afs){
  Vec_st prior(STNUM,1);
  if(afs.size()==0){
    double val = 1.0/(STNUM);
    for(int z = 0;z < STNUM;z++){
      prior(z) = val;
    }
  }else{    
    for(int z = 0;z < STNUM;z++){
      prior(z) = get_mass(afs,z,STNUM);
    }
  }
  prior = (1.0/prior.sum())*prior;
  return(prior);
}
//calcurate ln(exp(lna) + exp(lnb))
double log_sum_exp(double lna,double lnb){
  double ans;
  if(lna > lnb){
    ans = lna + log(exp(lnb-lna) + 1);
  }else{
    ans = lnb + log(exp(lna-lnb) + 1);
  }
  return(ans);
}
  //all element refresh (including umat uimat)
bool CQF::arefresh(Vec theta){
  refresh(theta);
  bool fail = false;
  auto aref_process = [&](auto &ctm){
    TIME(ctm.arefresh(theta),"CTM AREF")
    //parameter reset
    log_pe = 0;
    double slog_pe = 0;
    //cout<<data.size()<<endl;  
    double std_gen = 30;
    Vec_st one_vec = element_each_double(1.0,STNUM,1);
    Vec_st gam_vec = Vec_st::Zero(STNUM,1);
    double est_snum = (one_vec.transpose()*ctm.make_exp_tr(std_gen)*one_vec)(0,0);
    save_flag = (abs(est_snum - STNUM) > 1);
    //cout<<"data size " <<data.size()<<endl;
    for(int d = 0;d < data.size();d++){
      double pe;
      if(data.tpl(d).size() > 1){
	//pre pab
	TIME(ctm.lpe_refresh(data.pabl(d),data.tpl(d),pz),"LPE");
	pe = ctm.pe_get();
	if(pe < 0 || pe > 1.0){fail = true;cout<<"pe"<<pe<<endl;}
	//cout<<"likelihood "<<pe<<endl;
      }
      TIME(
	   log_pe += log(pe);
	   ctm.kapd_refresh();
	   if(afs_opt){
	     gam_vec = gam_vec + ctm.make_gam_vec(pz);
	   },"KAPD");
    }
    //slc prior
    log_pe += -dom_beta*slc*slc/2;
    TIME(
	 if(afs_opt){
	   cout<<pz<<endl;
	   gam_vec = gam_vec +(beta-1)*one_vec;
	   pz = (1/gam_vec.sum())*gam_vec;
	 }
	 ctm.tkap_mat_refresh(); 
	 fs_vec = ctm.make_fs_vec(); 
	 np_vec = ctm.make_np_vec(); 
	 nm_vec = ctm.make_nm_vec();
	 ,"FSVEC")
    if(fs_vec.real().minCoeff() < 0){
      if(abs(fs_vec.real().minCoeff())/fs_vec.real().maxCoeff() > 0.01){
	fail = true;cout<<"fs"<<endl;
	cout<<fs_vec.real().minCoeff()<<endl;
      }
    }
    if(np_vec.real().minCoeff() < 0){
      if(abs(np_vec.real().minCoeff())/np_vec.real().maxCoeff() > 0.01){
	fail = true;cout<<"np"<<endl;
	cout<<np_vec.real().minCoeff()<<endl;
      }
    }
    if(nm_vec.real().minCoeff() < 0){
      if(abs(nm_vec.real().minCoeff())/nm_vec.real().maxCoeff() > 0.01){
	fail = true;cout<<"nm"<<endl;
	cout<<nm_vec.real().minCoeff()<<endl;
      }
    }
    //cout<<"test"<<endl;
    //cout<<dom<<endl;
    //cout<<dom_beta<<endl;   
    //cout<<(dom_beta-1) * log(dom) + (dom_beta - 1) * log(1-dom)<<endl;
  };
  aref_process(ctmr);
  return(fail);
}
//caluculate ll h 
double CQF::llh(){
  Vec_st q = Vec_st::Zero(1,1);
  auto process = [&](auto &ctm){
    q += ctm.rc_vec*fs_vec;
    q += ctm.lrp_vec*np_vec; 
    q += ctm.lrm_vec*nm_vec;
  };
  process(ctmr);
  return(q(0,0) -slc*slc*dom_beta/2);
}
double CQF::llh_dpop(){
  Vec_st q = Vec_st::Zero(1,1); 
  //parameter cdonversion
  auto process = [&](auto &ctm){
    q += ctm.rcdneff_vec*fs_vec;
    q += ctm.lrpdneff_vec*np_vec;
    q += ctm.lrmdneff_vec*nm_vec;
  };
  process(ctmr);
  return(q(0,0));
}
double CQF::llh_dslc(){
  Vec_st q = Vec_st::Zero(1,1);
  auto process = [&](auto &ctm){
    //parameter cdonversion
    q += ctm.rcdslc_vec*fs_vec;
    q += ctm.lrpdslc_vec*np_vec;
    q += ctm.lrmdslc_vec*nm_vec;
  };
  process(ctmr);
  //slc prior
  return(q(0,0) -slc*dom_beta);
}
double CQF::llh_ddom(){
  Vec_st q = Vec_st::Zero(1,1);
  auto process = [&](auto &ctm){  
    //parameter cdonversion
    q += ctm.rcddom_vec*fs_vec;
    q += ctm.lrpddom_vec*np_vec;
    q += ctm.lrmddom_vec*nm_vec;
  };
  process(ctmr);
  return(q(0,0));
}
Mat_st CQF::make_info_mat(Vec theta){
  ctmr.arefresh(theta);
  ctmr.dd_refresh(theta);
  double dslcslc = 0;
  dslcslc += (ctmr.rcdslcslc_vec*fs_vec)(0,0);
  dslcslc += (ctmr.lrpdslcslc_vec*np_vec)(0,0);
  dslcslc += (ctmr.lrmdslcslc_vec*nm_vec)(0,0);
  double ddomdom = 0;
  ddomdom += (ctmr.rcddomdom_vec*fs_vec)(0,0);
  ddomdom += (ctmr.lrpddomdom_vec*np_vec)(0,0);
  ddomdom += (ctmr.lrmddomdom_vec*nm_vec)(0,0);
  double dslcdom = 0;
  dslcdom += (ctmr.rcdslcdom_vec*fs_vec)(0,0);
  dslcdom += (ctmr.lrpdslcdom_vec*np_vec)(0,0);
  dslcdom += (ctmr.lrmdslcdom_vec*nm_vec)(0,0);
  Mat_st info_mat = Mat_st::Zero(2,2);
  info_mat(0,0) = dslcslc;
  info_mat(0,1) = dslcdom;
  info_mat(1,0) = dslcdom;
  info_mat(1,1) = ddomdom;
  return(info_mat);
}
double CQF::next_slc(){
  /*
  double a = pop*H*(nsp - nsm);
  double b = fs*genpt/H;
  double c = sqrt(pow(a,2)+pow(b,2));
  double d;
  if(b > 0){
    d = (a + c)/b;
  }else{
    d = (a - c)/b;
  }
  double ans = log(d)/(pop*H);
  */
  double ans = 0;
  return(ans);
}
double CQF::next_pop(){
  /*
  double ans = genpt*fs/(H*H*(nsp+nsm));
  */
  double ans = 1.0;
  return(ans);
}
double CQF::get_log_pe(){
  return(log_pe);
}
Vec_st CQF::get_pz(){
  return(pz);
}
void CQF::print_loaded_data(){
  for(int d = 0; d < data.size(); d++){
    cout<<"data number: "<<d<<endl;
    cout<<"apl"<<endl;
    auto apl = data.apl(d);
    for(auto aitr = apl.begin(); aitr != apl.end(); aitr++){
      cout<<"\t"<<*aitr;
    }
    cout<<endl;
    cout<<"bpl"<<endl;
    auto bpl = data.bpl(d);
    for(auto bitr = bpl.begin(); bitr != bpl.end(); bitr++){
      cout<<"\t"<<*bitr;
    }
    cout<<endl;
    cout<<"tpl"<<endl;
    auto tpl = data.tpl(d);
    for(auto titr = tpl.begin(); titr != tpl.end(); titr++){
      cout<<"\t"<<*titr;
    }
    cout<<endl;
  }
}
//need :
// ctm.nc_pget(i) ctm.fcc_dgam(sigam)
// fsc nspc
//ctm.drefresh(data.epl(d),data.ecpl(d),data.tpl(d));
