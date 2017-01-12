#ifndef WDT
#define WDT
#include<vector>
using namespace std;
class Data{
  //number of data set (a,b,t)
  int dsize;
  //alpha list
  vector<vector<int>> apll;
  //beta list
  vector<vector<int>> bpll;
  //list of time of process
  vector<vector<double>> tpll;
  //discretized state number of x
  //confirm if refferd data nuber exist
  void dnum_check(int d);
public:
  Data();
  void clear();
  void add(vector<int> _apl,vector<int> _bpl,vector<double> _tpl);
  int size();
  vector<int> apl(int d);
  vector<int> bpl(int d);
  vector<double> tpl(int d);
};
#endif
