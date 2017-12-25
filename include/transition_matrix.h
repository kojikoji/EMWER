
#ifndef TM
#define TM
//#include<Matrix.h>
#define EIGEN_NO_DEBUG // コード内のassertを無効化．
#define EIGEN_DONT_PARALLELIZE // 並列を無効化．
#define EIGEN_MPL2_ONLY // LGPLライセンスのコードを使わない．

#include <Eigen/Core>
//#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 

#include <chrono>
#include <vector>
using  Mat_st = Eigen::MatrixXd;
using  Mat_stc = Eigen::MatrixXcd;
using  Vec_st = Eigen::MatrixXd;
using Vec_ist = Eigen::MatrixXd;
//typedef Matrix<std::complex<double>> Vec_stc;
using TVec_st =  Eigen::MatrixXd;
using Vec_mst =  Eigen::MatrixXd;
using Vec_3d = Eigen::MatrixXd;

template<class T>
using Mat_gen = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;

using namespace std;
template<class T>
class CTM;
// computation of transition matrix
template<class T>
class CTM{
  // number of discretization of AF
  int STNUM;
  //STNUM-1
  int MSTNUM;
  //STNUM-2
  int ISTNUM;
  //discretization resolution
  double dlt;
  //time discretization num band used in safe_nf_refresh
  double safe_disc = 10;
  //parameters
  double slc;
  double neff;
  double dom;
  double genpt = 1;
  // p(ad,bd|d=1..dsize)
  T pe;
  //p(alpha,beta|endpoint)
  vector<Vec_st> pabx;
  //probability of transition from z to e after t for z = 0:stnum-1
  //forward probability P(a(1:k),b(1:k),xk=i)
  vector<Mat_gen<T> > fwd;
  //backward probability P(a(k:D),b(k:D),xk=i)
  vector<Mat_gen<T> > bwd;
  //d pse /drij
  //sym*dg colwise is r
  // dependecy of x
  Mat_gen<T> u_mat;
  Mat_gen<T> ui_mat;
  Mat_gen<T> kapd_mat;
  Mat_gen<T> kap_mat;
  Mat_gen<T> tkap_mat;
  Mat_gen<T> ld;
  //data
  //number of endpoint data
  int dsize;
  //time point
  vector<double> tpl;
  //u*sum(td*kd)
  double to_xi(int i);
  Mat_gen<T> exp_eg(double t);
  Mat_gen<T> etr_make(double t);
  Mat_gen<T> diag(Mat_gen<T> vec);
  bool test_flag = true;
public:
  //c mean center
  Vec_st rc_vec;
  //p mean plus
  Vec_mst rp_vec;
  //m mean minus
  Vec_mst rm_vec;
  //p mean plus
  Vec_mst lrp_vec;
  //m mean minus
  Vec_mst lrm_vec;
  //c mean center
  Vec_st rcdneff_vec;
  //p mean plus
  Vec_mst lrpdneff_vec;
  //m mean minus
  Vec_mst lrmdneff_vec;
  //c mean center
  Vec_st rcdslc_vec;
  //p mean plus
  Vec_mst lrpdslc_vec;
  //m mean minus
  Vec_mst lrmdslc_vec;
  //c mean center
  Vec_st rcddom_vec;
  //p mean plus
  Vec_mst lrpddom_vec;
  //m mean minus
  Vec_mst lrmddom_vec;
  //c mean center
  Vec_st rcdslcslc_vec;
  //p mean plus
  Vec_mst lrpdslcslc_vec;
  //m mean minus
  Vec_mst lrmdslcslc_vec;
  //c mean center
  Vec_st rcdslcdom_vec;
  //p mean plus
  Vec_mst lrpdslcdom_vec;
  //m mean minus
  Vec_mst lrmdslcdom_vec;
  //c mean center
  Vec_st rcddomdom_vec;
  //p mean plus
  Vec_mst lrpddomdom_vec;
  //m mean minus
  Vec_mst lrmddomdom_vec;
  CTM();
  void init(int snum);
  //r and rd is refreshed for new parameter  
  bool param_refresh(Vec_3d theta);
  //secoudary differencitation
  bool dd_refresh(Vec_3d theta);
  //eigen deconposition done
  void eigen_refresh();
  //refresh pi_, mi_ and ci_ mat
  void i_mat_refresh();
  //refresh epmx
  void epmx_refresh();
  //all element refresh (including umat uimat)
  bool arefresh(Vec_3d theta);
  void tkap_mat_refresh();
  //refresh k_uv and pse
  //refresh pze and return vector of pzx1*pzx2*... for z
  bool lpe_refresh(vector<Vec_st> _pabl,vector<double> _tpl,Vec_st pz);
  double safe_lpe_refresh(vector<Vec_st> _pabl,vector<double> _tpl,Vec_st pz);
  //
  void kapd_refresh();
  //refresh N and F ,argument is qz = pz*px1*px2*...*pc1*pc2*...
  void nf_refresh();
  //in dangerous selection, sufficient statistics are safly estimated
  void safe_nf_refresh();
  //get parameters
  double pe_get();  
  //part of r and dr which depend on parameters
  double f_rq(double x);
  double f_rq_dneff(double x);
  double f_rq_dslc(double x);
  double f_rq_ddom(double x);
  double f_rq_dslcslc(double x);
  double f_rq_ddomdom(double x);
  double f_rq_dslcdom(double x);
  double f_rqp(double x);
  double f_rqp_dneff(double x);
  double f_rqp_dslc(double x);
  double f_rqp_ddom(double x);
  double f_rqp_dslcslc(double x);
  double f_rqp_dslcdom(double x);
  double f_rqp_ddomdom(double x);
  double f_r(double x);
  double f_q(double x);
  double f_irdneff(double x);
  double f_iqdneff(double x);
  double f_irdslc(double x);
  double f_iqdslc(double x);
  double f_irddom(double x);
  double f_iqddom(double x);
  Mat_st make_r_mat();
  Mat_st make_u_ldn_ui(double n);
  Mat_st make_exp_tr(double t);
  Mat_st make_exp_tr_num(double t,int num);
  Vec_st make_fs_vec();
  Vec_st make_np_vec();
  Vec_st make_nm_vec();
  Vec_st make_gam_vec(Vec_st &pz);
};
template <class T>
void print_vec(vector<T> vec);
template <class T>
Mat_gen<T> element_each(T el,int N,int M);

double no_minus(double val);
Mat_st no_minus(Mat_st mat);
double binom(double p,int a,int b);
#endif
