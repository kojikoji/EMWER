#ifndef QF
#define QF
#include<iostream>
#include<data.h>
#include<transition_matrix.h>
typedef Matrix<double> Mat_st;
typedef Matrix<double> Vec_st;
typedef Matrix<double> TVec_st;
typedef Matrix<double> Vec_mst;
typedef Matrix<double> Vec_2d;
//caluculate q function
class CQF{
  // number of discretization of AF
  int STNUM;
  //STNUM-1
  int MSTNUM;
  double H;
  CTM ctm,ctmc;
  Data data,datac;
  Vec_2d sigam;
  Vec_2d sigamc;
  double pop,slc,genpt;
  double fs;
  double nsp;
  double nsm;
  double fsc;
  double nspc;
  double nsmc;
  double log_pe;
  double vdlh,vdlh_dpop,vdlh_dslc;
  Vec_st pz;
 public:
  bool save_flag = false;
  CQF();
  //set genpt
  void init(double _genpt,int snum,vector<double> afs);
  //give data which have stactm. and end and time
  void data_add(vector<int> apl,vector<int> bpl,vector<double> tpl,vector<int> acpl,vector<int> bcpl,vector<double> tcpl);
  void data_clear();
  //refresh sigam for pop slc
  void to_sigam();
  Vec_st prior_make(vector<double> afs);
  //ctm.and ctm. is refershed foctm.new pactm.metectm. 
  void refresh(double _pop,double _slc);
  //all element refresh (including umat uimat)
  void arefresh(double _pop,double _slc);
  //caluculate llh
  double llh();
  double llh_dsig();
  double llh_dgam();
  double llh_dpop();
  double llh_dslc();
  double next_slc();
  double next_pop();
  /*
  //for continuous discrete fix
  void cont_dis_set(double _pop,double _slc);
  double lpdis(int i, int j);
  double lpdis_dpop(int i, int j);
  double lpdis_dslc(int i, int j);
  void drefresh(double _pop,double _slc);
  double dlh();
  double dlh_dpop();
  double dlh_dslc();
  */
  //get log likelihood
  double get_log_pe();
  double set_sflag(bool _save_flag);
};
#endif
