//extern STNUM Vec_2d
#include<iostream>
#include<q_function.h>
#include<lbfgsb.hpp>
//optimizer for Start End Data class
class WFE{
  vector<bool> opt_flag;
  //parameter population and slection 
  double pop,slc;
  CQF cqf;
  double llh;
  double llh_dnpop;
  double llh_dnslc;
  //rate of gamma and sigma to 10
  double popr,slcr;
  //generation per time
  double genpt;
  void drefresh();
  //refresh llh and llh_d
  void refresh();
  //minimizer of llh for r 
  class CMR{
  public:
    CMR(WFE &_wfe,double c,vector<bool> opt_flag);
    int operator()(const vector<double>& x,double& fn,vector<double>& gr);
    double _c;
    int _i;
    vector<bool> _opt_flag;
    WFE &upper;
  };
  /*
  class DMR{
  public:
    DMR(WFE &_wfe,double c,vector<bool> opt_flag);
    int operator()(const vector<double>& x,double& fn,vector<double>& gr);
    double _c;
    int _i;
    vector<bool> _opt_flag;
    WFE &upper;
  };
  */
  //operate clh and minimize the lh of clh 
  template<class T>
  void steepest_descent(T clh);
public:
  WFE();
  //give first value of parameters
  void ftheta_give(Vec_2d _theta);
  //give data which have start and end and time
  void data_add(vector<int> apl,vector<int> bpl,vector<double> tpl,vector<int> acpl,vector<int> bcpl,vector<double> tcpl);
  void data_clear();
  void optimize(Vec_2d theta,vector<bool> _opt_flag,double _genpt,int snum,vector<double> afs);
  //fix error of discrete and continuous model
  //void cont_dis_fix(Vec_2d ftheta,vector<bool> _opt_flag,double _genpt);
  //get optimized prameters
  Vec_2d ans_get();
  //accuracy of optimized
  //vector<double> acc_get();
  //fixed error between contiuous and discrete model
  //vector<double> fixedans_get();
  // -2log(likelihood ratio)
  double stat_chi2();
};
