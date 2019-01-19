#ifndef WDT
#define WDT
#include<vector>
//#include<Matrix.h>
#define EIGEN_NO_DEBUG // コード内のassertを無効化．
#define EIGEN_DONT_PARALLELIZE // 並列を無効化．
#define EIGEN_MPL2_ONLY // LGPLライセンスのコードを使わない．

#include <Eigen/Core>
//#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 
#include <map>
typedef Eigen::MatrixXd Vec_st;

using namespace std;
class Data{
  //number of data set (a,b,t)
  int dsize;
  //lines of dat file
  vector<string> line_list;
  //alpha list
  vector<vector<int>> apll;
  //beta list
  vector<vector<int>> bpll;
  //list of time of process
  vector<vector<double>> tpll;
  // replicate idx
  vector<int> unique_ridx_l;
  // sync col idx for each replicate
  vector<vector<int>> col_idx_ll;
  // generation of each column in sync file
  vector<int> sync_t_l;
  // base idx map in sync file field
  std::map<std::string, int> sync_base_idx_dict;
  // snum
  int snum = -1;
  //list of binom dist
  vector<vector<Vec_st>> pabll;
  bool sync_mode;
  //discretized state number of x
  //confirm if refferd data nuber exist
  void dnum_check(int d);
  void load_sync_line(int n);
  void load_dat_line(int n);
public:
  int line_num;
  Data();
  void clear();
  void read_file(std::string file_name);
  void read_sync(std::string filename, vector<int> sync_ridx_l, vector<int> _sync_t_l);
  void load_line(int n, bool add_flag=false);
  void load_all();
  void change_allele();
  void add(vector<int> _apl,vector<int> _bpl,vector<double> _tpl);
  void pabx_refresh(int _snum);
  int size();
  vector<int> apl(int d);
  vector<int> bpl(int d);
  vector<double> tpl(int d);
  vector<Vec_st> pabl(int d);
};
vector<string> split(const string &str, char sep);
#endif
