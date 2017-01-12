#include<iostream>
#include<vector>
#include<data.h>
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
