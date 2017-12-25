#include<boundary.hpp>
BDR::BDR(Vec origin_theta,vector<int> _snum_list,bool dom_opt,double dom_min,double dom_max){
  snum_list = _snum_list;
  Vec theta = origin_theta;
  cout<<"in"<<endl;
  dom_lower = theta(2);
  dom_upper = theta(2);
  if(dom_opt){
    dom_lower = dom_min;
    dom_upper = dom_max;
  }
  for(int snum : snum_list){
    //新しい状態数に変更
    cout<<"snum"<<snum<<endl;
    ctmr.init(snum);
    double upper = 1;
    double lower = -1;
    for(double n_dom = dom_lower;n_dom <= dom_upper; n_dom+= unit_dom){
      //ドミナンス変更
      theta(2) = n_dom;
      //上限を決める
      double n_upper = 0;
      for(double n_slc = 0; n_slc < 1; n_slc += unit){
	theta(1) = n_slc;
	n_upper = theta(1) - 2*unit;
	//失敗したらその2step前が境界
	if(ctmr.param_refresh(theta)){
	  break;
	}
      }
      //cout<<"up slc:"<<n_upper<<"\t dom:"<<n_dom<<endl;
      if(n_upper < upper){
	upper = n_upper;
      }
      //下限を決める
      double n_lower = 0;
      for(double n_slc = 0; n_slc > -1; n_slc -= unit){
	theta(1) = n_slc;
	n_lower = theta(1) + 2*unit;
	//失敗したらその2step前が境界
	if(ctmr.param_refresh(theta)){
	  break;
	}
      }
      //cout<<"low slc:"<<n_lower<<"\t dom:"<<n_dom<<endl;
      if(n_lower > lower){
	lower = n_lower;
      }
    }
    upper_list.push_back(upper);
    lower_list.push_back(lower);
  }
}
bool BDR::increase_state(){
  bool success = snum_index < snum_list.size() - 1;
  if(success){
    snum_index++;
  }
  return(!success);
}
bool BDR::decrease_state(){
  bool success = snum_index > 0;
  if(success){
    snum_index--;
  }
  return(!success);
}
  
int BDR::get_snum(){
  return(snum_list[snum_index]);
}

double BDR::get_lower(){
  return(lower_list[snum_index]);
}

double BDR::get_upper(){
  return(upper_list[snum_index]);
}
double BDR::get_dlower(){
  return(dom_lower);
}

double BDR::get_dupper(){
  return(dom_upper);
}
