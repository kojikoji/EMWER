#include<iostream>
#include<sstream>
#include <fstream>
#include<string>
#include<vector>
#include<data.h>
vector<string> split(const string &str, char sep)
{
    vector<std::string> v;
    std::stringstream ss(str);
    std::string buffer;
    while( std::getline(ss, buffer, sep) ) {
        v.push_back(buffer);
    }
    return v;
}
void Data::dnum_check(int d){
  
  try{
    if(d > dsize-1 || d < 0){
      throw "Data index is out of scope";
    }
  }
  catch(char *str){
    std::cout <<str<<std::endl;
  }
}
Data::Data(){
  dsize = 0;
  line_num = 0;
}
void Data::clear(){
  dsize = 0;
  tpll.clear();
  apll.clear();
  bpll.clear();
}
void Data::add(vector<int> _apl,vector<int> _bpl,vector<double> _tpl){
  dsize++;
  tpll.push_back(_tpl);
  //classify repl
  apll.push_back(_apl);
  bpll.push_back(_bpl);
}
double binom(double p,int a,int b){
  int n = a+b;
  //combination and p
  double ans = 1; 
  for(int i = 0;i < a;i++){
    ans *= ((double)(n-i)/(a-i))*p;
  }
  for(int i = 0;i < b;i++){
    ans *= 1-p;
  }  
  if(p == 1){
    if(b==0){
      ans = 1;
    }else{
      ans = 0;
    }
  }
  if(p == 0){
    if(a == 0){
      ans = 1;
    }else{
      ans = 0;
    }
  }
  return(ans);
}
void Data::pabx_refresh(int _snum){
  if(snum != _snum){
    pabll.clear();
    snum = _snum;
    double dlt = 1.0/snum;
    for(int d = 0; d < apll.size(); d++){
      vector<Vec_st> pabl;
      for(int t = 0;t < apll[d].size();t++){
	Vec_st qx=  Vec_st::Zero(snum,1);;
	for(int i = 1;i < snum-1;i++){
	  double xi = 0.5*(dlt*i + dlt*(i+1));
	  qx(i) = binom(xi,apll[d][t],bpll[d][t]);
	  if(qx(i) < 1.0e-5){
	    qx(i) = 0;
	  }
	}
	qx(0) = binom(0,apll[d][t],bpll[d][t]);
	qx(snum-1) = binom(1.0,apll[d][t],bpll[d][t]);
	pabl.push_back(move(qx));
      }
      pabll.push_back(pabl);
    }
  }
}   
int Data::size(){
  return(dsize);
}
vector<int> Data::apl(int d){
  dnum_check(d);
  return(apll[d]);
}
vector<int> Data::bpl(int d){
  dnum_check(d);
  return(bpll[d]);
}
vector<double> Data::tpl(int d){
  dnum_check(d);
  return(tpll[d]);
}
vector<Vec_st> Data::pabl(int d){
  dnum_check(d);
  return(pabll[d]);
}
void Data::read_file(std::string filename){
  ifstream ifs( filename );
  if(!ifs ) {
    cout << "Error:Input data file not found" << endl;
    exit(-1);
  }
  vector<vector<vector<double> > > data;
  int n = 0;
  string str;
  while( getline( ifs, str ) ){
    line_list.push_back(str);
    line_num++;
  }
}

void Data::load_line(int n, bool add_flag){
  if(!add_flag){
    clear();
  }
  string line = line_list[n];
  vector<string> data_list = split(line,'\t');
  for(auto itr = data_list.begin(); itr != data_list.end(); itr++){
    vector<string> val_list = split(*itr,',');
    vector<int> apl;
    vector<int> bpl;
    vector<double> tpl;
    for(int d = 0;d < val_list.size()/3;d++){
      double alpha = std::stoi(val_list[3*d]);
      apl.push_back(alpha);
      double beta = std::stoi(val_list[3*d+1]);
      bpl.push_back(beta);
      double time = std::stoi(val_list[3*d+2]);
      tpl.push_back(time);
    }
    add(apl,bpl,tpl);
  }
}

void Data::load_all(){
  for(int n = 0; n < line_num; n++){
    load_line(n,true);
  }
}
void Data::change_allele(){
  vector<vector<int>> tmpll;
  tmpll = apll;
  apll = bpll;
  bpll = tmpll;
  int tsnum = snum;
  snum = -1;
  pabx_refresh(tsnum);
}
