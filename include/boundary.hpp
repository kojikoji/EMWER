#ifndef BD
#define BD
#include<iostream>
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
//#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

//caluculate q function
class BDR{
  // number of discretization of AF
  vector<int> snum_list;
  vector<double> lower_list;
  vector<double> upper_list;
  CTM<std::complex<double>> ctmr;
  int snum_index = 0;
  double unit = 0.01;
  double unit_dom = 0.1;
  double dom_lower;
  double dom_upper;
 public:
  BDR(Vec theta,vector<int> _snum_lsit,bool dom_opt=false,double dom_min=0,double dom_max=1.0);
  bool increase_state();
  bool decrease_state();
  vector<double> get_boundary();
  double get_lower();
  double get_upper();
  double get_dlower();
  double get_dupper();
  int get_snum();
};
#endif
