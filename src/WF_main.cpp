#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <algorithm>
#include <chrono>
#include <wright_fisher_estimater.h>
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
vector<string> split(const string &str, char delim){
  vector<string> res;
  size_t current = 0, found;
  while((found = str.find_first_of(delim, current)) != string::npos){
    res.push_back(string(str, current, found - current));
    current = found + 1;
  }
  res.push_back(string(str, current, str.size() - current));
  return res;
}
vector<vector<vector<double > > > read_sim_file(string filename,bool aflag=false){
  ifstream ifs( filename );
  if(!ifs ) {
    cout << "Error:Input data file not found" << endl;
    exit(-1);
  }
  vector<vector<vector<double> > > data;
  int n = 0;
  string str;
  while( getline( ifs, str ) ){
    if(str.size()<1){
      cout<<"Input is empty"<<endl;
      exit(-1);
    }
    vector<vector<double> > one_record;
    string token;
    istringstream stream( str );
    int m = 0;
    while( getline( stream, token, '\t' ) ) {
      vector<double> same_init;
      string subtoken;
      istringstream substream( token );
      int k = 0;
      while( getline( substream, subtoken, ',' ) ) {
	stringstream subss;
	double val;
	subss << subtoken;
	subss >>val;
	same_init.push_back(val);
	k++;
      }
      //cout<<"k"<<k<<endl;
      one_record.push_back(same_init);
      m++;
    }
    //cout<<"m"<<m<<endl;
    //if not all data for parameter divede data
    if(!aflag || data.size()== 0){
      data.push_back(one_record);
    }else{
      data[0].insert(data[0].end(),one_record.begin(),one_record.end());
    }
    n++;
  }
  return(data);
}
void sample_register(WFE& wfe,vector<vector<double> > one_record,vector<vector<double> > one_recordc,int genpt,bool add_f=false){
  for(int i = 0;i < one_record.size();i++){
    vector<double> data = one_record[i];
    vector<double> datac;
    if(one_recordc.size() > 0){
       datac = one_recordc[i];
    }else{
      datac = {};
    }
    vector<int> apl,bpl,acpl,bcpl;
    vector<double> tpl,tcpl;
    for(int d = 0;d < data.size()/3;d++){
      double alpha = (double)data[3*d];
      apl.push_back(alpha);
      double beta = (double)data[3*d+1];
      bpl.push_back(beta);
      double time = (double)data[3*d+2]/genpt;
      tpl.push_back(time);
    }
    for(int d = 0;d < datac.size()/3;d++){
      double alpha = (double)datac[3*d];
      acpl.push_back(alpha);
      double beta = (double)datac[3*d+1];
      bcpl.push_back(beta);
      double time = (double)datac[3*d+2]/genpt;
      //cout<<time<<endl;
      tcpl.push_back(time);
    }
    wfe.data_add(apl,bpl,tpl,acpl,bcpl,tcpl);
  }
}
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
Vec_2d estimate_process(WFE &wfe,vector<bool> opt_flag,double psize,double slc,double genpt,int snum,vector<double> afs){
    int hy_count = 0;
    Vec_2d ftheta(2,1);
    ftheta(0) = psize;
    ftheta(1) = slc;
    Vec_2d theta(2,1);
    wfe.optimize(ftheta,opt_flag,genpt,snum,afs);
    theta = wfe.ans_get();
    return(theta);
}
void output_process(WFE &wfe,ofstream &ofs,Vec_2d theta,double psize,string genl,double msec){
    cout <<theta.str()<<endl;
    ofs <<1;
    ofs <<"\t"<<psize;
    ofs <<"\t"<<genl;
    for(int k = 0;k <theta.size();k++){
      ofs<<"\t"<<theta(k);
    }
    cout <<"Selection coefficient"<<":"<<"\t";
    cout <<theta(0)<<endl;
    cout <<"Population size"<<":"<<"\t";
    cout <<theta(1)<<endl;
    double chisq = wfe.stat_chi2();
    cout<<"Likelihood ratio:"<<"\t"<<chisq<<endl;
    ofs<<"\t"<<chisq;
    cout<<"Run time:"<<"\t"<<msec<<endl;
    ofs<<"\t"<<msec;
    ofs <<endl;
}
   
// args: psize generations gamma sigma
int main(int argc,char* argv[]){
  //parsecommand line and get population, generations, selection, and filename
  cmdline::parser p;
  p.add<int>("discretization", 'd', "number of discretization of AF",false,50);
  p.add<int>("population", 'p', "population of system",false,500);
  p.add<double>("selection", 's', "selection coefficient",false,0);
  p.add<string>("generation", 'g', "list of number of generations for experiments",false,"");
  p.add<string>("output", 'o', "output file");
  p.add<string>("input", 'i', "input file  for case");
  p.add<string>("control", 'c', "input file of control",false);
  p.add<string>("afs", 'q', "file of allele frequency spectrum",false);
  p.add("fslc", 'f', "flag for slection with no value");
  p.add("fpop", 'l', "flag for slection");
  p.add("allfile", 'a', "estimate one parameter from file of all");
  p.add("help", 0, "print help");
  if (!p.parse(argc, argv)||p.exist("help")){
    std::cout<<p.error_full()<<p.usage();
    return 0;
  }
  int snum = p.get<int>("discretization");
  int psize = p.get<int>("population");
  double slc = p.get<double>("selection");
  string fname = p.get<string>("input");
  string ofname = p.get<string>("output");
  string genl = p.get<string>("generation");
  bool cflag = true;
  string fnamec;
  if(p.exist("control")){
    fnamec = p.get<string>("control");
  }else{
    cflag = false;
  }
  vector<double> afs;
  if(p.exist("afs")){
    afs = afs_reader(p.get<string>("afs"));
  }
  ofstream ofs;
  ofs.open(ofname);
  vector<bool> opt_flag;
  if(p.exist("fpop")){
    opt_flag.push_back(true);
  }else{
    opt_flag.push_back(false);
  }
  if(p.exist("fslc")){
    opt_flag.push_back(true);
  }else{
    opt_flag.push_back(false);
  }  
  
  vector<double> null_theta = {1,0};
  double genpt = (double)psize*2;
  vector<vector<vector<double> > > data = read_sim_file(fname,p.exist("allfile"));
  vector<vector<vector<double> > > datac;
  if(cflag){
    datac  = read_sim_file(fnamec);
  }
  for(int d = 0;d < data.size();d++){//data
    WFE wfe;
    cout<<data[d].size()<<endl;
    //register each row
    vector<vector<double > > one_record = data[d];
    vector<vector<double > > one_recordc = {};
    sample_register(wfe,one_record,one_recordc,genpt,true);
    auto start = chrono::system_clock::now();      // 計測スタート時刻を保存
    Vec_2d theta = estimate_process(wfe,opt_flag,psize,slc,genpt,snum,afs); //パラメータ推定
    auto end = std::chrono::system_clock::now();       // 計測終了時刻を保存
    auto dur = end - start;        // 要した時間を計算
    auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    // 要した時間をミリ秒（1/1000秒）に変換して表示
    std::cout << msec << " milli sec \n";
    output_process(wfe,ofs,theta,psize,genl,msec);
  }
  return(0);
}
//-p 1200 -g 60 -s 0.005 -i tmp/sim/simg_p1200_gen60_slc0.005.data -o tmp/o.txt
