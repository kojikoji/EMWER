#include<iostream>
#include<transition_matrix.h>
#include<q_function.h>
CQF::CQF():data(),datac(){
  sigam.resize(2,1);
  sigamc.resize(2,1);
}
void CQF::data_add(vector<int> apl,vector<int> bpl,vector<double> tpl,vector<int> acpl,vector<int> bcpl,vector<double> tcpl){
  data.add(apl,bpl,tpl);
  datac.add(acpl,bcpl,tcpl);
}
void CQF::data_clear(){
  data.clear();
  datac.clear();
}
void CQF::to_sigam(){
   //sigma
  sigam(0) = sqrt(genpt/pop);
  //gamma
  sigam(1) = slc*pop;
  //for control
  //sigma
  sigamc(0) = sqrt(genpt/pop);
  //gamma
  sigamc(1) = 0;
 }
void CQF::init(double _genpt,int snum,vector<double> afs){
  STNUM = snum;
  MSTNUM = snum -1;
  H = (double)1.0/MSTNUM;
  ctm.init(snum);
  genpt = _genpt;
  save_flag = false;
  pz = prior_make(afs);
}
  //f~ and ctm.and ctm. is refreshed for R new pararmeter. 
void CQF::refresh(double _pop,double _slc){
  pop = _pop;
  slc = _slc;
  to_sigam();
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
      prior(z) = afs[z];
    }
  }
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
void CQF::arefresh(double _pop,double _slc){
  //refresh pop slc sigam
  refresh(_pop,_slc);
  ctm.arefresh(sigam);
  fs = 0;
  nsp = 0;
  nsm = 0;
  fsc = 0;
  nspc = 0;
  nsmc = 0;
  log_pe = 0;
  double slog_pe = 0;
  for(int d = 0;d < data.size();d++){
    double dpe = 1;
    if(data.tpl(d).size() > 1){
      if(save_flag){
	//cout<<"In safe mode"<<endl;
	dpe = ctm.safe_lpe_refresh(data.apl(d),data.bpl(d),data.tpl(d),pz);
      }else{
	dpe =  ctm.lpe_refresh(data.apl(d),data.bpl(d),data.tpl(d),pz);
      }
    }
    //cout<<"d"<<d<<endl;
    if(dpe>0){
      //wcout<<"dpe"<<dpe<<endl;
      //cout<<"pz_vec"<<lpz_vec.str()<<endl;
      //cout<<"pzc_vec"<<lpzc_vec.str()<<endl;
    }
    log_pe += dpe;
    //slog_pe += sdpe;
    double fd = 0;
    double ndp = 0;
    double ndm = 0;
    fs += ctm.tf_get(0);
    nsp += ctm.n_pget(0);
    nsm += ctm.n_mget(0);
    /*
    cout <<"fs: "<<fs;
    cout <<" nsp: "<<nsp;
    cout <<" nsm: "<<nsm<<endl;
    */
  }
  //cout<<" same logpe? "<<log_pe<<"\t"<<slog_pe<<endl;
  //cout <<"fs: "<<fs;
  //cout <<" nsp: "<<nsp;
  //cout <<" nsm: "<<nsm<<endl;
}
//caluculate ll h 
double CQF::llh(){
  double q = 0; 
  q += -(ctm.f_p(sigam) + ctm.f_m(sigam))*fs;
  q += log(ctm.f_p(sigam))*nsp;
  q += log(ctm.f_m(sigam))*nsm;
  return(q);
}
double CQF::llh_dsig(){
  double q_sigam = 0;
  q_sigam += ctm.fc_dsig(sigam)*fs;
  q_sigam += ctm.fp_dsig(sigam)*nsp;
  q_sigam += ctm.fm_dsig(sigam)*nsm;
  return(q_sigam);
}
double CQF::llh_dgam(){
  double q_sigam = 0;
  q_sigam += ctm.fc_dgam(sigam)*fs;
  q_sigam += ctm.fp_dgam(sigam)*nsp;
  q_sigam += ctm.fm_dgam(sigam)*nsm;
  return(q_sigam);
}
double CQF::llh_dpop(){
  //parameter cdonversion
  double sig_dpop = -3*sqrt(genpt)*pow(pop,-1.5)/2;
  double gam_dpop = -2*slc/pow(pop,2);
  double llh_dpop = llh_dsig()*sig_dpop + llh_dgam()*gam_dpop;
  return(llh_dpop);
}
double CQF::llh_dslc(){
  //parameter cdonversion
  double sig_dslc = 0;
  double gam_dslc = 1/pop;
  double llh_dslc = llh_dsig()*sig_dslc + llh_dgam()*gam_dslc;
  return(llh_dslc);
}
double CQF::next_slc(){
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
  return(ans);
}
double CQF::next_pop(){
  double ans = genpt*fs/(H*H*(nsp+nsm));
  return(ans);
}
double CQF::get_log_pe(){
  return(log_pe);
}
//need :
// ctm.nc_pget(i) ctm.fcc_dgam(sigam)
// fsc nspc
//ctm.drefresh(data.epl(d),data.ecpl(d),data.tpl(d));
