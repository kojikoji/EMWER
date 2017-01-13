//extern STNUM Vec_2d
#include<iostream>
#include<q_function.h>
#include<wright_fisher_estimater.h>
#include<lbfgsb.hpp>
#include<math.h>
void WFE::refresh(){
   cqf.refresh(pop,slc);
  llh = cqf.llh();
  //parameter cdonversion
  double pop_dnpop = popr;
  double slc_dnslc = slcr;
  //llh_d
  double llh_dpop = cqf.llh_dpop();
  double llh_dslc = cqf.llh_dslc();
  llh_dnpop = llh_dpop*pop_dnpop;
  llh_dnslc = llh_dslc*slc_dnslc;
}
  
WFE::CMR::CMR(WFE &_wfe,double c,vector<bool> opt_flag):upper(_wfe),_c(c),_i(0),_opt_flag(opt_flag){}
int WFE::CMR::operator()(const vector<double>& x,double& fn,vector<double>& gr){
  //unnormalized prameter(n,s) theta is caliculated
  //x is changed
  upper.pop = x[0]*upper.popr;
  upper.slc = x[1]*upper.slcr;
  //compute llh and llh_d
  upper.refresh();
  //subsititute for fn and gr
  fn = -upper.llh;
  gr.assign(x.size(),0);
  if(_opt_flag[0]){
    gr[0] += -upper.llh_dnpop;
  }
  if(_opt_flag[1]){
    gr[1] += -upper.llh_dnslc;
  }
  return(0);
}
template<class T>
void WFE::steepest_descent(T clh){//compute likelhood
  Lbfgsb minimizer;
  minimizer.set_eps(1.0e-12);
  minimizer.set_maxit(30);
  popr = 100;//pop/10;
  slcr = 1.0e-4;//slc/10;
  vector<double> ntheta0;
  ntheta0.push_back(pop/popr);
  ntheta0.push_back(slc/slcr);  
  minimizer.minimize(ntheta0, clh);
  vector<double> ntheta = minimizer.best_x();
  pop = ntheta[0]*popr;
  slc = ntheta[1]*slcr;
}

WFE::WFE(){}
void WFE::data_add(vector<int> apl,vector<int> bpl,vector<double> tpl,vector<int> acpl,vector<int> bcpl,vector<double> tcpl){
  cqf.data_add(apl,bpl,tpl,acpl,bcpl,tcpl);
}
void WFE::data_clear(){
  cqf.data_clear();
}
void WFE::optimize(Vec_2d ftheta,vector<bool> _opt_flag,double _genpt,int snum,vector<double> afs){
  opt_flag = _opt_flag;
  genpt = _genpt;
  cqf.init(genpt,snum,afs);
  pop = 2*ftheta(0);
  slc = ftheta(1)/2;
  if(pop == 0){
    pop = 1;
  }
  double est_upper = 40;
  if(opt_flag[0]==true){
    slc = 0;
    est_upper = 1000;
  }
  double ppop = pop;
  double pslc = slc;
  int count = 0;
  double llh,pllh;
  cqf.arefresh(pop,slc);
  llh = cqf.get_log_pe();
  pllh = llh;
  // limitation of slc
  double ih = slc;
  double il = slc;
  double oh = 0.5;//100/(2*pop);
  double ol = -0.5;//-100/(2*pop);
  bool end_flag;
  int safe_count = 0;
  do{
    //refresh parameter
    ppop = pop;
    pslc = slc;
    pllh = llh;
    //cqf.arefresh(pop,slc);
    //pllh = cqf.get_log_pe();
    if(!opt_flag[0] && opt_flag[1]){
      slc = cqf.next_slc();
    }else if(opt_flag[0] && !opt_flag[1]){
      pop = cqf.next_pop();
    }else{
      //M step: optimize theta
      double c = 2.0;
      CMR cmr((*this),c,opt_flag);
      steepest_descent(cmr);
    }
    int arefct = 0;
    cqf.arefresh(pop,slc);
    llh=cqf.get_log_pe();
    //llh = cqf.get_log_pe();
    if(slc > oh){
      break;
    }
    if(slc < ol){
      break;
    }
    end_flag = (fabs(pop - ppop) < 1.0e-4 && fabs(slc  - pslc) < 1.0e-5);
    end_flag = end_flag || std::isnan(llh) || llh < pllh;
    if(llh < pllh){
      pop = ppop;
      slc = pslc;
      llh = pllh;
    }
    if(end_flag && opt_flag[1]){      
      pslc = slc;
      pllh = llh;
      double slcpl = pslc + 0.05;
      cqf.arefresh(pop,slcpl);
      double llhpl = cqf.get_log_pe();
      double slcmn = pslc - 0.05;
      cqf.arefresh(pop,slcmn);
      double llhmn = cqf.get_log_pe();
      llh = max(llh,cqf.get_log_pe());
      if(max(llhmn,llhpl) > pllh){
	cout<<"Trans to safe mode in"<<endl;
	end_flag = false;
	if(llhpl > llhmn){
	  slc = slcpl;
	}else{
	  slc = slcmn;
	}
	cqf.save_flag = true;
	cqf.arefresh(pop,slc);
      }else{
	cout<<"OK estimation is completed!"<<endl;
      }
    }
    if(cqf.save_flag){
      safe_count++;
    }
    if(safe_count>10){
      cout<<"too much safe mode"<<endl;
      break;
    }
    if(std::isnan(llh) || llh < pllh){
      slc = pslc;
      break;
    }
    count++;
    
   }while(!end_flag && count < est_upper);
}
Vec_2d WFE::ans_get(){
  Vec_2d ans(2,1);
  ans(0) = pop/2;
  ans(1) = slc*2;
  return(ans);
}
// -2log(likelihood ratio)
double WFE::stat_chi2(){
  cqf.arefresh(pop,slc);
  double log_pe = cqf.get_log_pe();
  cqf.arefresh(pop,0);
  double log_pe0 = cqf.get_log_pe();
  double ans = -2*(log_pe0-log_pe);
  return(ans);
}
