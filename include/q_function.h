#ifndef QF
#define QF
#include<iostream>
#include<data.h>
#include<transition_matrix.h>

typedef Eigen::MatrixXd Mat_st;
//typedef Matrix<std::complex<double>> Mat_stc;
typedef Eigen::MatrixXd Vec_st;
//typedef Matrix<std::complex<double>> Vec_stc;
typedef Eigen::MatrixXd TVec_st;
typedef Eigen::MatrixXd Vec_mst;
typedef Eigen::MatrixXd Vec_3d;
typedef Eigen::MatrixXd Vec;

#define EIGEN_NO_DEBUG 
#define EIGEN_DONT_PARALLELIZE 
#define EIGEN_MPL2_ONLY 

#include <Eigen/Core>
#include <Eigen/Eigenvalues> 

//caluculate q function
class CQF{
  // number of discretization of AF
  int STNUM;
  //STNUM-1
  int MSTNUM;
  double H;
  //CTM<double> ctmr;
  CTM<std::complex<double>> ctmr;
  Data data;
  double dom_beta;
  double beta;
  double genpt;
  Vec_ist fs_vec;
  Vec_ist np_vec;
  Vec_ist nm_vec;
  double log_pe;
  double vdlh,vdlh_dpop,vdlh_dslc;
  double dom,slc;
  Vec_st pz;
 public:
  bool save_flag = false;
  bool afs_opt = false;
  CQF();
  //set genpt
  void init(double _genpt,int snum,vector<double> afs, double _beta=2, double _dom_beta=0);
  //change allele 1 - x 
  void change_allele();
  //give data which have stactm. and end and time
  void load_data(Data& _data);
  void data_clear();
  Vec_st prior_make(vector<double> afs);
  //ctm.and ctm. is refershed foctm.new pactm.metectm. 
  bool refresh(Vec theta);
  //all element refresh (including umat uimat)
  bool arefresh(Vec theta);
  //caluculate llh
  double llh();
  double llh_dpop();
  double llh_dslc();
  double llh_ddom();
  Mat_st make_info_mat(Vec theta);
  double next_slc();
  double next_pop();
  //get log likelihood
  double get_log_pe();
  Vec_st get_pz();
  double set_sflag(bool _save_flag);
  void print_loaded_data();
};
#endif
