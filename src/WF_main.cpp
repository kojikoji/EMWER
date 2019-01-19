

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
    cerr << "Error:Input afs file not found" << endl;
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
    cerr << "Error:Input afs is too small" << endl;
    exit(-1);
  }    
  return(afs);
}
void output_process(ofstream &ofs,Vec_3d theta, double chisq ,double msec){
    ofs<<theta(0)<<"\t"<<theta(1)<<"\t"<<theta(2);
    ofs<<"\t"<<chisq;
    ofs<<"\t"<<msec;
    ofs <<endl;
}
string estimate_process(vector<bool> opt_flag,double psize,double slc, double dom,double genpt,int snum, double beta, double dom_beta, double lh_thresh, vector<double> afs,Data& data, bool avd_flag,BDR bdr, string oafs_file = "", bool variant_flag=false){
    WFE wfe;
    Vec_3d ftheta(3,1);
    ftheta(0) = psize;
    ftheta(1) = slc;
    ftheta(2) = dom;
    if(opt_flag[0]){
      wfe.optimize_pop(data,ftheta,opt_flag,genpt,snum,afs,bdr,beta, dom_beta, avd_flag);
    }else{
      wfe.optimize(data,ftheta,opt_flag,genpt,snum,afs,bdr,beta, dom_beta, avd_flag, variant_flag);
    }
    Vec_3d theta = wfe.ans_get();
    Vec_3d est_var = wfe.estimate_variance();
    double chisq = wfe.stat_chi2();
    stringstream sfs;
    sfs<<theta(0)<<"\t"<<theta(1)<<"\t"<<theta(2);
    sfs<<"\t"<<est_var(0)<<"\t"<<est_var(1);
    sfs<<"\t"<<chisq;
    sfs<<"\t"<<est_var(2);
    sfs <<endl;
    //afs estimation
    if(oafs_file != ""){
      ofstream ofs_afs;
      ofs_afs.open(oafs_file);
      wfe.afs_get();
      ofs_afs<<wfe.afs_get();
      ofs_afs.close();
    }
    return(sfs.str());
}

vector<int> csv2int_vec(string str){
  vector<string> sval_vec = split(str,',');
  vector<int> ival_vec;
  for(string sval : sval_vec){
    ival_vec.push_back(stoi(sval));
  }
  return(ival_vec);
}
   
// args: psize generations gamma sigma
int main(int argc,char* argv[]){
  //parsecommand line and get population, generations, selection, and filename
  cmdline::parser p;
  p.add<int>("discretization", 'd', "number of discretization of AF",false,20);
  p.add<int>("threads", 't', "number of threads",false,1);
  p.add<double>("beta", 'b', "hyper parameter of Dirichlet",false,1.5);
  p.add<double>("dom_beta", 'c', "hyper parameter of Beta dist of dominance",false,0);
  p.add<double>("lh_thresh", 'z', "threshold criterion for likelihood convergence",false,10);
  p.add<int>("population", 'p', "population of system",false,250);
  p.add<double>("selection", 's', "selection coefficient",false,0);
  p.add<double>("dominance", 'n', "dominance",false,0.5);
  p.add<double>("dommin", 'y', "dominance",false,-0.5);
  p.add<double>("dommax", 'x', "dominance",false,1.5);
  p.add<string>("generation", 'g', "list of number of generations for experiments",false,"");
  p.add<string>("sync_repl", '\0', "list of replicate indexes for each sync",false,"");
  p.add<string>("sync_genl", '\0', "list of generations for experiments",false,"");
  p.add<string>("sync", '\0', "file name of sync file",false,"");
  p.add<string>("output", 'o', "output file",false,"/Users/kojimayasuhiro/Projects/wfeed/log.txt");
  p.add<string>("input", 'i', "input file  for case",false,"/Users/kojimayasuhiro/Projects/wfeed/tmp/validation/test_real.dat");
  p.add<string>("afs", 'q', "file of allele frequency spectrum",false);
  p.add<string>("oafs", 'r', "output file of allele frequency spectrum",false,"");
  p.add<string>("fopt", 'f', "flag for optimization p: population s: selection d:dominance a: allele frequency spectrum",false,"s");
  p.add("allfile", 'a', "estimate one parameter from file of all");
  p.add("fix", 'l', "fix discretization");
  p.add("variant", 'v', "define the selection value for variant allele");
  p.add("help", 0, "print help");
  if (!p.parse(argc, argv)||p.exist("help")){
    std::cout<<p.error_full()<<p.usage();
    return 0;
  }
  int snum = p.get<int>("discretization");
  double beta = p.get<double>("beta");
  double dom_beta = p.get<double>("dom_beta");
  double lh_thresh = p.get<double>("lh_thresh");
  int psize = p.get<int>("population");
  int num_thread = p.get<int>("threads");
  double slc = p.get<double>("selection");
  double dom = p.get<double>("dominance");
  double dom_min = p.get<double>("dommin");
  double dom_max = p.get<double>("dommax");
  string fname = p.get<string>("input");
  string ofname = p.get<string>("output");
  string genl = p.get<string>("generation");
  string sync_fname = p.get<string>("sync");
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
  if(p.exist("sync")){
    vector<int> sync_genl = csv2int_vec(p.get<string>("sync_genl"));
    vector<int> sync_repl = csv2int_vec(p.get<string>("sync_repl"));
    if(sync_genl.size() != sync_repl.size()){
      cerr<<"Error: Inconsistent size of sync_repl and sync_genl"<<endl;
      exit(-1);
    }
    data.read_sync(sync_fname, sync_repl, sync_genl);
  }else{
    data.read_file(fname);
  }
  int itr_num = data.line_num;
  double genpt = 1;
  Vec_3d ftheta(3,1);
  ftheta(0) = psize;
  ftheta(1) = slc;
  ftheta(2) = dom;
  
  int init_snum = 20;
  if(p.exist("fix") || opt_flag[0] || ftheta(1) != 0){
    cout<<"snum fixed"<<endl;
    init_snum = snum;
  }
  vector<int> snum_array = {20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250};
  vector<int> snum_list;
  for(int add_snum : snum_array){
    if(add_snum > snum){ break;}
    snum_list.push_back(add_snum);
  }
  BDR bdr(ftheta,snum_list,opt_flag[2],dom_min,dom_max);
  if(p.exist("allfile")){
    data.load_all();
    string rlt = estimate_process(opt_flag, psize, slc, dom, genpt, init_snum, beta, dom_beta, lh_thresh, afs, data,avd_flag,bdr,p.get<string>("oafs"));
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
      rlt_vec[d] = estimate_process(opt_flag, psize, slc, dom, genpt, init_snum, beta, dom_beta, lh_thresh, afs, data_vec[t_num], avd_flag,bdr, "", p.exist("variant"));
    }
    string rlt = accumulate(rlt_vec.begin(),rlt_vec.end(),string(""));
    ofs << rlt;
  }
  ofs.close();
  return(0);
}
