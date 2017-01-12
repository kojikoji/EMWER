#include <vector>
void sample_make(std::vector<double> theta,double time,int data_num);
void data_setter(double l,double u);
void data_adder(double s,double t,double time);
std::vector<double> em(std::vector<double> ftheta,std::vector<bool> opt_flag,double genpt);
void data_clear();
void sample_add(double s,double rand,double time);
void sample_condition(double t,double lower,double upper,std::vector<double> smpl_theta);
double hypo_test(std::vector<double> null_theta,std::vector<double> sample_theta,double alpha,int try_num);
double estimate_variance(std::vector<double> null_theta,int try_num);
void sample_make_i(std::vector<double> theta,double time,int data_num,std::vector<double> init_list);
