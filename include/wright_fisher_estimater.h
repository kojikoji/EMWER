//extern STNUM Vec_2d
#include<iostream>
#include<q_function.h>
#include <boundary.hpp>
#include<lbfgsb.hpp>
#include<data.h>
typedef Eigen::MatrixXd Vec_3d;
//optimizer for Start End Data class
class WFE{
  double eps_dom = 1.0e-15;
  vector<bool> opt_flag;
  //parameter population and slection 
  double pop,slc;
  Vec_3d theta;
  Vec_3d pz_vec;
  CQF cqf;
  int snum;
  double lh;
  double init_lh;
  double init_snum;
  double qfun;
  double qfun_dnpop;
  double qfun_dnslc;
  double qfun_dpop;
  double qfun_dslc;
  double qfun_ddom;
  int sign_slc;
  //rate of gamma and sigma to 10
  double popr,slcr;
  //generation per time
  double genpt;
  void drefresh();
  //refresh llh and llh_d
  bool refresh(Vec_st ntheta);
  //minimizer of llh for r 
  class CMR{
  public:
    CMR(WFE &_wfe,double c,vector<bool> opt_flag);
    int operator()(vector<double>& x,double& fn,vector<double>& gr);
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
    void steepest_descent(T clh,vector<double> sbound,vector<double> dbound);
public:
  WFE();
  //give first value of parameters
  void ftheta_give(Vec_3d _theta);
  void load_data(Data& data);
  //give data which have start and end and time
  void data_clear();
  void optimize(Data& data,Vec_3d theta,vector<bool> _opt_flag,double _genpt,int snum,vector<double> afs, BDR bdr, double beta=2, double dom_beta = 2, bool avd_flag = false);
  //fix error of discrete and continuous model
  //void cont_dis_fix(Vec_3d ftheta,vector<bool> _opt_flag,double _genpt);
  //get optimized prameters
  Vec_3d ans_get();
  Vec_3d afs_get();
  //accuracy of optimized
  //vector<double> acc_get();
  //fixed error between contiuous and discrete model
  //vector<double> fixedans_get();
  // -2log(likelihood ratio)
  double stat_chi2();
  Vec_3d estimate_variance();
};
