#include <vector>
#include <fstream>
#include <omp.h>
#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <algorithm>
#include <chrono>
#include <wright_fisher_estimater.h>
#include <boundary.hpp>
#include <cmdline.h>
#define TIME(PROCESS,TAG)  start5 = chrono::system_clock::now();;PROCESS;end5 = std::chrono::system_clock::now();dur5 = end5 - start5;msec5 = std::chrono::duration_cast<std::chrono::microseconds>(dur5).count();std::cout<<TAG <<"\t"<< msec5 << " milli sec \n";
#ifndef TIME
#define TIME(PROCESS,TAG)  PROCESS;
#endif
auto start5 = chrono::system_clock::now();
auto end5 = std::chrono::system_clock::now();       // 計測終了時刻を保存
auto dur5 = end5 - start5;        // 要した時間を計算
auto msec5 = std::chrono::duration_cast<std::chrono::milliseconds>(dur5).count();
using namespace std;
static double t_smpl = 0.1;
static double gamma_smpl = 5;
static double sigma_smpl = 1;
template<typename T_n>
T_n uniform_distribution(mt19937 &_mt,T_n min,T_n max){//[min,max]
    return((T_n)((max - min)*(_mt() - _mt.min())/(_mt.max() - _mt.min())) + min);
}
random_device rd;
mt19937 gen(rd());
vector<double> afs_reader(string afs_file){
  ifstream ifs(afs_file);
  if(!ifs ) {
    cout << "Error:Input afs file not found" << endl;
    exit(-1);
  }
  vector<double> afs;
  int count = 0;
  while(!ifs.eof()){
    double val;
    ifs >> val;
    afs.push_back(val);
  }
  if(afs.size()<2){
    cout << "Error:Input afs is too small" << endl;
    exit(-1);
  }    
  return(afs);
}
void output_process(ofstream &ofs,Vec_3d theta, double chisq ,double msec){
  //cout <<theta.str()<<endl;
    ofs<<theta(0)<<"\t"<<theta(1)<<"\t"<<theta(2);
    //cout <<"Population size"<<":"<<"\t";
    //cout <<theta(0)<<endl;
    //cout <<"Selection coefficient"<<":"<<"\t";
    //cout <<theta(1)<<endl;
    //cout <<"Dominance"<<":"<<"\t";
    //cout <<theta(2)<<endl;
    //cout<<"Likelihood ratio:"<<"\t"<<chisq<<endl;
    ofs<<"\t"<<chisq;
    //cout<<"Run time:"<<"\t"<<msec<<endl;
    ofs<<"\t"<<msec;
    ofs <<endl;
}
string estimate_process(vector<bool> opt_flag,double psize,double slc, double dom,double genpt,int snum, double beta, double dom_beta,vector<double> afs,Data& data, bool avd_flag,BDR bdr, string oafs_file = ""){
    auto start = chrono::system_clock::now();
    // 計測スタート時刻を保存
    WFE wfe;
    Vec_3d ftheta(3,1);
    ftheta(0) = psize;
    ftheta(1) = slc;
    ftheta(2) = dom;
    wfe.optimize(data,ftheta,opt_flag,genpt,snum,afs,bdr,beta, dom_beta,avd_flag);
    Vec_3d theta = wfe.ans_get();
    auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
    auto dur = end - start;        // 要した時間を計算
    auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
      // 要した時間をミリ秒（1/1000秒）に変換して表示
    //std::cout << msec << " milli sec \n";
    Vec_3d est_var = wfe.estimate_variance();
    double chisq = wfe.stat_chi2();
    stringstream sfs;
    sfs<<theta(0)<<"\t"<<theta(1)<<"\t"<<theta(2);
    //cout <<"Population size"<<":"<<"\t";
    //cout <<theta(0)<<endl;
    //cout <<"Selection coefficient"<<":"<<"\t";
    //cout <<theta(1)<<endl;
    //cout <<"Dominance"<<":"<<"\t";
    //cout <<theta(2)<<endl;
    //cout<<"Likelihood ratio:"<<"\t"<<chisq<<endl;
    sfs<<"\t"<<est_var(0)<<"\t"<<est_var(1);
    sfs<<"\t"<<chisq;
    //cout<<"Run time:"<<"\t"<<msec<<endl;
    sfs<<"\t"<<est_var(2);
    sfs <<endl;
    //afs estimation
    if(oafs_file != ""){
      ofstream ofs_afs;
      ofs_afs.open(oafs_file);
      wfe.afs_get();
      ofs_afs<<wfe.afs_get();
      ofs_afs.close();
      //cout<<wfe.afs_get().str();
    }
    return(sfs.str());
}
   
// args: psize generations gamma sigma
int main(int argc,char* argv[]){
  //parsecommand line and get population, generations, selection, and filename
  cmdline::parser p;
  cout<<"begin"<<endl;
  p.add<int>("discretization", 'd', "number of discretization of AF",false,20);
  p.add<int>("threads", 't', "number of threads",false,1);
  p.add<double>("beta", 'b', "hyper parameter of Dirichlet",false,1.5);
  p.add<double>("dom_beta", 'c', "hyper parameter of Beta dist of dominance",false,0);
  p.add<int>("population", 'p', "population of system",false,250);
  p.add<double>("selection", 's', "selection coefficient",false,0);
  p.add<double>("dominance", 'n', "dominance",false,0.5);
  p.add<double>("dommin", 'y', "dominance",false,-0.5);
  p.add<double>("dommax", 'x', "dominance",false,1.5);
  p.add<string>("generation", 'g', "list of number of generations for experiments",false,"");
  p.add<string>("output", 'o', "output file",false,"/Users/kojimayasuhiro/Projects/wfeed/tmp/validation/rlt_sample.est");
  p.add<string>("input", 'i', "input file  for case",false,"/Users/kojimayasuhiro/Projects/wfeed/tmp/validation/test_real.dat");
  p.add<string>("afs", 'q', "file of allele frequency spectrum",false);
  p.add<string>("oafs", 'r', "output file of allele frequency spectrum",false,"");
  p.add<string>("fopt", 'f', "flag for optimization p: population s: selection d:dominance a: allele frequency spectrum",false,"s");
  p.add("allfile", 'a', "estimate one parameter from file of all");
  p.add("fix", 'l', "fix discretization");
  p.add("help", 0, "print help");
  if (!p.parse(argc, argv)||p.exist("help")){
    std::cout<<p.error_full()<<p.usage();
    return 0;
  }
  int snum = p.get<int>("discretization");
  double beta = p.get<double>("beta");
  double dom_beta = p.get<double>("dom_beta");
  int psize = p.get<int>("population");
  int num_thread = p.get<int>("threads");
  double slc = p.get<double>("selection");
  double dom = p.get<double>("dominance");
  double dom_min = p.get<double>("dommin");
  double dom_max = p.get<double>("dommax");
  string fname = p.get<string>("input");
  string ofname = p.get<string>("output");
  string genl = p.get<string>("generation");
  bool cflag = true;
  string fnamec;
  vector<double> afs;
  if(p.exist("afs")){
    afs = afs_reader(p.get<string>("afs"));
  }
  ofstream ofs;
  ofs.open(ofname);
  string opt_code = p.get<string>("fopt");
  vector<bool> opt_flag = {false,false,false,false};
  if(opt_code.find("p") != string::npos){
    opt_flag[0] = true;
  }
  if(opt_code.find("s") != string::npos){
    opt_flag[1] = true;
  }  
  if(opt_code.find("d") != string::npos){
    opt_flag[2] = true;
  }  
  if(opt_code.find("a") != string::npos){
    opt_flag[3] = true;
  }
  bool avd_flag = false;
  if(opt_code.find("v") != string::npos){
    avd_flag = true;
  }
  Data data;
  data.read_file(fname);
  int itr_num = data.line_num;
  double genpt = 1;
  Vec_3d ftheta(3,1);
  ftheta(0) = psize;
  ftheta(1) = slc;
  ftheta(2) = dom;
  
  vector<int> snum_list;
  int init_snum = 20;
  if(p.exist("fix") || opt_flag[0]){
    cout<<"snum fixed"<<endl;
    init_snum = snum;
  } 
  for(int add_snum = init_snum; add_snum < snum+1;add_snum += 10){
    snum_list.push_back(add_snum);
  }
  BDR bdr(ftheta,snum_list,opt_flag[2],dom_min,dom_max);
  if(p.exist("allfile")){
    data.load_all();
    string rlt = estimate_process(opt_flag, psize, slc, dom, genpt, init_snum, beta, dom_beta, afs, data,avd_flag,bdr,p.get<string>("oafs")); //パラメータ推定
    ofs << rlt;
  }else{
    vector<Data> data_vec;
    for(int i = 0; i < num_thread; i++){
      data_vec.push_back(data);
    }
    vector<string> rlt_vec;
    rlt_vec.resize(itr_num);
#pragma omp parallel for num_threads(num_thread)
    for(int d = 0;d < itr_num;d++){//data
      int t_num = omp_get_thread_num();
      data_vec[t_num].load_line(d);
      rlt_vec[d] = estimate_process(opt_flag, psize, slc, dom, genpt, init_snum, beta, dom_beta, afs, data_vec[t_num], avd_flag,bdr); //パラメータ推定
    }
    string rlt = accumulate(rlt_vec.begin(),rlt_vec.end(),string(""));
    ofs << rlt;
  }
  ofs.close();
  return(0);
}
//-p 1200 -g 60 -s 0.005 -i tmp/sim/simg_p1200_gen60_slc0.005.data -o tmp/o.txt
