#include <iostream>
#include <math.h>
#include <vector>
#include <complex>
#include <fstream>
#include <random>
#include <string>
#include "cmdline.h"
using namespace std;
static random_device rd;
static mt19937 gen(rd());
template<typename T_n>
T_n uniform_distribution(mt19937 &_mt,T_n min,T_n max){//[min,max]
    return((T_n)((max - min)*(_mt() - _mt.min())/(_mt.max() - _mt.min())) + min);
}
double rand_normal( double mu, double sigma ){
  double u1 = (double)uniform_distribution<double>(gen,0,1);
  double u2 = (double)uniform_distribution<double>(gen,0,1);
  double z=sqrt( -2.0*log(u1) ) * sin( 2.0*M_PI*u2 );
  return mu + sigma*z;
}
double next_x(double x,int psize,double s){
  double eta = x/(1 - s*(1-x));
  double var = psize*eta*(1-eta);
  double mean = psize*eta;
  double apbi = rand_normal(mean,sqrt(var));       // approximatio of binominal  
  double nx = apbi/psize;
  if(nx < 0){			// exception process
    nx = 0;
  }else if(nx >1){
    nx= 1;
  }
  return(nx);
}
int poisson_sampler(int n,double p){
  int ct = 0;
  for(int i = 0;i < n;i++){
    double x = (double)uniform_distribution<double>(gen,0,1);
    if(x <= p){
      ct++;
    }
  }
  return(ct);
}
double one_sampler(double slc,int all_num,int gen_num,double &x,int fld,double flc,ofstream& ofs,int pgen = 0){
  //double x = 0.5;
  //int rgen_num = uniform_distribution<double>(gen,0.2*gen_num,1.8*gen_num);      
  int rgen_num = gen_num;
  for(int g = 0;g < rgen_num;g++){
    x = next_x(x,all_num,slc);
  }
  
  double rate;
  if(flc > 0){
    rate = (double)uniform_distribution<double>(gen,1-flc,1+flc);
  }else{
    rate = 1;
  }
  double n = fld*rate;
  int alpha = poisson_sampler(n,x);
  int beta = n - alpha;
  ofs <<alpha<<","<<beta<<","<<rgen_num+pgen;
  return(x);
  //cout <<slc<<"alpha "<<alpha<<"\t"<<"beta "<<beta<<"\t"<<"x "<<x<<endl;
}

void wright_fisher_sampling(int dnum,double slc,int all_num,vector<int> genl,int fld,double flc,double lower,double upper,ofstream &ofs){
  if(!(lower <= upper && lower >= 0 && upper <= 1)){
    cout<<"ERROR initial boundi is wrong"<<endl;
  }
  double t=1;
  vector<double> x_list;
  x_list.assign(2,2);
  double init = uniform_distribution<double>(gen,lower,upper);
  for(int i= 0;i <dnum;i++){
    //double init = uniform_distribution<double>(gen,0.8,0.9);
    //double init = 0.5;//uniform_distribution<double>(gen,0.1,0.9);
    double x = init;
    x = one_sampler(slc,all_num,0,x,fld,flc,ofs);
    for(int k = 1;k < genl.size();k++){
      int gen_num = genl[k]-genl[k-1];
      if(k != 0){
	ofs <<",";
      }
      x = one_sampler(slc,all_num,gen_num,x,fld,flc,ofs,genl[k-1]);
    }
    if(i != dnum-1){
      ofs <<"\t";
    }
  }
  ofs<<endl;
}
vector<string> split(const string &str, const string &delim){
  vector<string> res;
  size_t current = 0, found, delimlen = delim.size();
  while((found = str.find(delim, current)) != string::npos){
    res.push_back(string(str, current, found - current));
    current = found + delimlen;
  }
  res.push_back(string(str, current, str.size() - current));
  return res;
}
vector<int> to_int(vector<string> strl){
  vector<int> intl;
  for(auto strlit = strl.begin();strlit != strl.end();strlit++){
    intl.push_back(atoi((*strlit).c_str()));
  }
  return(intl);
}
vector<double> to_double(vector<string> strl){
  vector<double> doublel;
  for(auto strlit = strl.begin();strlit != strl.end();strlit++){
    doublel.push_back(atof((*strlit).c_str()));
  }
  return(doublel);
}
void afs_make(int disc,double lower,double upper,string afsf){
  if(afsf != ""){
    ofstream oafsf(afsf);
    double h = (double)1.0/(disc-1);
    int lwind = (lower + h/2)/h;
    int upind = (upper + h/2)/h;
    double uniprob = (double)1.0/(upind - lwind + 1);
    for(int i = 0;i < disc;i++){
      if(i >= lwind && i <= upind){
	oafsf<<uniprob<<endl;
      }else{
	oafsf<<0<<endl;
      }
    }
  }
}
int main(int argc,char* argv[]){
  cmdline::parser p;
  //p.add("hoge", 'h', "hoge flag with no value");
  p.add<int>("replicate", 'r', "number of replicate");
  p.add<int>("number", 'n', "number of run");
  p.add<int>("population", 'p', "population of system");
  p.add<string>("generation", 'g', "list of number of generations for experiments");
  p.add<double>("selection", 's', "selection coefficient");
  p.add<int>("folds", 'f', "folds number");
  p.add<double>("fluctuate", 'u', "ratio of fluctuate of depth",false,0);
  p.add<string>("bound", 'b', "boud for initial AF",false,"0,1");
  p.add<string>("afs", 'a', "af spectrum file",false,"");
  p.add<int>("disc", 'd', "discretization of afs",false,50);
  p.add<string>("output", 'o', "output file");
  p.add<string>("control", 'c', "output file of control",false);
  p.add("help", 0, "print help");
  // add other flags...

  if (!p.parse(argc, argv)||p.exist("help")){
    std::cout<<p.error_full()<<p.usage();
    return 0;
  }
  string test = "1,2,3,4";
  split(test,",");
  int spnum = p.get<int>("number");//number of sample is 100
  int dnum = p.get<int>("replicate");
  int psize = p.get<int>("population")*2;
  vector<string> genl_str = split(p.get<string>("generation"),",");
  vector<int> genl = to_int(genl_str);
  double slc = p.get<double>("selection")/2;
  int fld = p.get<int>("folds");
  double flc = p.get<double>("fluctuate");
  double lower = to_double(split(p.get<string>("bound"),","))[0];
  double upper = to_double(split(p.get<string>("bound"),","))[1];
  string ofname = p.get<string>("output");
  string ofnamec = p.get<string>("control");
  ofstream ofs(ofname);
  for(int i = 0;i < spnum;i++){
    wright_fisher_sampling(dnum,slc,psize,genl,fld,flc,lower,upper,ofs);
  }
  string afsf = p.get<string>("afs");
  int disc = p.get<int>("disc");
  afs_make(disc,lower,upper,afsf);
  return(0);
}



