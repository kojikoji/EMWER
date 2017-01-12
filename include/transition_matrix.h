#ifndef TM
#define TM
#include<Matrix.h>
typedef Matrix<double> Mat_st;
typedef Matrix<double> Vec_st;
typedef Matrix<double> TVec_st;
typedef Matrix<double> Vec_mst;
typedef Matrix<double> Vec_2d;
// computation of transition matrix
class CTM{
  // number of discretization of AF
  int STNUM;
  //STNUM-1
  int MSTNUM;
  //discretization resolution
  double H;
  //time discretization num band used in safe_nf_refresh
  double safe_disc = 10;
  //pz * pze1 * pze2 * ...
  Vec_st q_vec;
  // sum z of qz
  double qe;
  // p(ad,bd|d=1..dsize)
  double pe;
  //sufficient statistics
  double fd,ndp,ndm;
  //error
  double epmx,ldmx;
  //p(alpha,beta|endpoint)
  vector<Vec_st> pabx;
  //probability of transition from z to e after t for z = 0:stnum-1
  //forward probability P(a(1:k),b(1:k),xk=i)
  vector<Vec_st > fwd;
  //backward probability P(a(k:D),b(k:D),xk=i)
  vector<Vec_st > bwd;
  //d pse /drij
  Vec_mst np_vec,nm_vec;
  Vec_st tf_vec;
   //matrix r differenciated by parameter
  vector<Vec_st> rdsig_c;
  vector<Vec_st> rdgam_c;
  //symetry matrix is this
  //of center
  Vec_st sym_c;
  //of side
  Vec_mst sym_s;
  //sym*dg colwise is r
  Vec_st dg;
  Vec_st sqdg;
  Vec_st isqdg;
  //c mean center
  Vec_st r_c;
  //p mean plus
  Vec_mst r_p;
  //m mean minus
  Vec_mst r_m;
  // dependecy of x
  Vec_st fiv;
  Mat_st k_uv;
  Mat_st u_mat;
  Mat_st ui_mat;
  Vec_st ld;
  //data
  //number of endpoint data
  int dsize;
  //alpha nuber of major fold
  vector<int> apl;
  //beta number of minor fold
  vector<int> bpl;
  //time point
  vector<double> tpl;
  //u*sum(td*kd)
  //use for clcuration of sufficient statistics
  Mat_st pi_mat;
  Mat_st mi_mat;
  Mat_st ci_mat;
  double to_xi(int i);
  Vec_st exp_eg(double t);
  void utk_refresh();
  double tpze_drc(int i);
  double tpze_drp(int i);
  double tpze_drm(int i);
  Mat_st etr_make(double t);
  template <class T>
  Matrix<T> diag(Matrix<T> vec);
public:
  CTM();
  void init(int snum);
  //r and rd is refreshed for new parameter  
  void refresh(Vec_2d theta);
  //eigen deconposition done
  void eigen_refresh();
  //refresh u_mat and ui_mat and dg
  void u_mat_refresh();
  //refresh pi_, mi_ and ci_ mat
  void i_mat_refresh();
  //refresh epmx
  void epmx_refresh();
  //all element refresh (including umat uimat)
  void arefresh(Vec_2d theta);
  void pabx_refresh();
  //refresh k_uv and pse
  //refresh pze and return vector of pzx1*pzx2*... for z
  double lpe_refresh(vector<int> _apl,vector<int> _bpl,vector<double> _tpl,Vec_st pz);
  double safe_lpe_refresh(vector<int> _apl,vector<int> _bpl,vector<double> _tpl,Vec_st pz);
  //refresh N and F ,argument is qz = pz*px1*px2*...*pc1*pc2*...
  void nf_refresh();
  //in dangerous selection, sufficient statistics are safly estimated
  void safe_nf_refresh();
  double tf_get(int i);
  double n_pget(int i);
  double n_mget(int i);
  double pse_get();  
  //part of r and dr which depend on parameters
  double f_i(int i);
  double f_p(Vec_2d theta);
  double f_m(Vec_2d theta);
  double fp_dsig(Vec_2d theta);
  double fm_dsig(Vec_2d theta);
  double fc_dsig(Vec_2d theta);
  double fc_dgam(Vec_2d theta);
  double fp_dgam(Vec_2d theta);
  double fm_dgam(Vec_2d theta);
};
#endif
