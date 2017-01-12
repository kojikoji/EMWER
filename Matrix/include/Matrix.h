#ifndef MAT
#define MAT
#include <sstream>
#include <vector>
#include <memory>
#include <typeinfo>
#include <blas.h>
#include <lapacke.h>
#include <iostream>
#include <string>
#include <execinfo.h>
// __cxa_demangle() はこのヘッダ
#include <cxxabi.h>

using namespace std;
static size_t upper_size = 1.0e8;
static size_t all_size = 0;
template <class T>
class Matrix;
template <class T>
class Matrix{
//error print
// バックトレース情報を文字列で取得する
vector<string> split(const string &str, const string &delim) const{
  vector<string> res;
  size_t current = 0, found, delimlen = delim.size();
  while((found = str.find(delim, current)) != string::npos){
    res.push_back(string(str, current, found - current));
    current = found + delimlen;
  }
  res.push_back(string(str, current, str.size() - current));
  return res;
}
vector<string> get_backtrace10() const{
    // backtrace()でバックトレース情報を取得できる
    const int trace_size = 10;
    void* trace[trace_size];
    int size = backtrace(trace, trace_size);

    // backtrace_symbols()でバックトレース情報を文字列に変換できる
    char** symbols = backtrace_symbols(trace, size);

    vector<string> result(symbols, symbols + size);

    // symbolsをfreeする必要がある
    free(symbols);

    return result;
}

// 取得したバックトレース文字列から関数名部分を切り出す
string cut_function_name_part(const string& raw_text) const{
  auto splited_text = split(raw_text," ");
    // '(' ～ '+' までを切り出せばいいみたい
  string func_name;
  if(splited_text.size() > 2){
    func_name = *(splited_text.end() - 3);
  }else{
    func_name = "";
  }      
  return(func_name);
}
// 関数名をデマングルする
string demangle_function_name(const string& mangled) const{
    // __cxa_demangle()でデマングルできる
    int status = 0;
    char* demangled = abi::__cxa_demangle(mangled.c_str(), 0, 0, &status);

    // demangledをfreeする必要がある
    string result = demangled ? demangled : mangled;
    free(demangled);

    if (status != 0) {
        std::ostringstream oss;
        oss << " [error:status=" << status << "]";
        result += oss.str();
    }

    return result;
}
void exception_exit(string erstr) const{
  void *array[10];
  size_t size;
  // get void*'s for all entries on the stack
  cout<<"ERROR!:"<<erstr<<endl;
  vector<string> backtrace_texts = get_backtrace10();

  std::cout << "Backtrace:\n";

  for (vector<string>::const_iterator
	 it = backtrace_texts.begin(), end = backtrace_texts.end();
       it != end;
       ++it) {
    string mangled_func_name = cut_function_name_part(*it);

    std::cout
      << ((!mangled_func_name.empty())
	  ? demangle_function_name(mangled_func_name)
	  : "???")
      << "\n";
  }
  exit(1);
}
  //list of memory poll
  static vector<vector<unique_ptr<T[]> > > pool_list;
  static vector<int> pool_size_list;
  int pool_index;
  unique_ptr<T[]> val_ptr;
  T val;
  int nrow; // mと同じ値
  int ncol; // mと同じ値
  int _size;
  //for blas and lapack 
  int one = 1; // mと同じ値
  double done = 1.0; // mと同じ値
  // for inverse 
  int lda; // mと同じ値
  int info; // 計算が成功すれば0を返す
  static int ipiv[10000]; // 要素数はm,nのうち小さい方とする
  int lwork; // nと同じ値
  static T work[10000]; // 要素数はlworkと同じ値
  // get matrix region
  int demmand_pool();
  // get matrix region
  void get_reg();
  // return matrix region for pool
  void return_reg();
  //some parameters depend on sizeis changed
  void size_set(int n,int m);
  //if size is 0, inform error
  void size_checker() const;
public:
  Matrix();
  Matrix(int n,int m);
  // initialize by arg_mat
  Matrix(const Matrix<T>& arg_mat);
  Matrix(Matrix<T>&& arg_mat);
  ~Matrix();
  T *begin() const;
  T *end() const;
  //if this have region, reg() return true, else false
  bool reg() const;
  //return pointer
  unique_ptr<T[]> give_ptr();
  T* get() const;
  //parameter change and get region
  void resize(int n,int m);
  //if lvalue, copy arg_mat value to this region  
  void operator=(const Matrix<T>& arg_mat);
  //if rvalue, val_ptr indicate arg_mat ptr
  void operator=(Matrix<T>&& arg_mat);
  // (i,j) return refference of Matrix(i,j)
  T &operator()(int i,int j) const;
  // (i) return refference of i th element on memory
  T &operator()(int i) const;
  //Mat+Mat
  Matrix<T> operator+(Matrix<T> &&arg_mat);
  // if arg_mat is lval, make temporaly 
  Matrix<T> operator+(const Matrix<T> &arg_mat);
  //if single element, this can be summed with T 
  T operator+(T a);
  // -Mat
  Matrix<T> operator-() const;
  // Mat - Mat
  Matrix<T> operator-(Matrix<T> &&arg_mat);
  // if arg_mat is lval, make temporaly 
  Matrix<T> operator-(const Matrix<T> &arg_mat);
  //if single element, this can be substracted with T 
  T operator-(T a);
  //Mat*a
  Matrix<T> operator*(double a);
  // Mat * Mat 
  Matrix<T> operator*(const Matrix<T>& arg_mat);
  Matrix<T> operator*(Matrix<T>&& arg_mat);
  // return inverse matrix
  Matrix<T> inverse();  
  //return transpose
  Matrix<T> transpose() const;
  //return Matrix whose all element = 0
  static Matrix<T> Zero(int n,int m);
  void print() const;
  // cast to string ,example cout <<Matrix<<endl;
  string str() const;
  //_size N*M
  int size() const;
  // number of rows
  int rows() const;
  //number of cols
  int cols() const;
  //take nth row
  Matrix<T> row(int n) const;
  //take mth col
  Matrix<T> col(int m) const;
  //cut from i th element to j th element
  Matrix<T> cut(int i,int j) const;
  //cut from i th element to j th col
  Matrix<T> col_cut(int i,int j) const;
  //cut from i th element to j th row
  Matrix<T> row_cut(int i,int j) const;
  //cast for double when nrow = 1 ncol = 1
  operator T();
};


// static
template<class T>
  Matrix<T> operator*(double a,const Matrix<T>& arg_mat);
  
static int ncount = 0;
//demmand index of pool for size n*m
template<class T>
int Matrix<T>::demmand_pool(){
  if(!(nrow > 0 && ncol > 0)){
    return(-1);
  }
  int ans;
  for(auto itr = pool_size_list.begin();itr != pool_size_list.end();++itr){
    if(*itr == nrow*ncol){
      pool_index = distance(pool_size_list.begin(),itr);
      return(pool_index);
    }
  }
  // there is no pool whose size is just, then make the pool
  vector<unique_ptr<T[] > > mempool;
  pool_list.push_back(move(mempool));
  pool_size_list.push_back(nrow*ncol);
  pool_index = pool_size_list.size() -1;
  return(pool_index);
}
// get matrix region
template<class T> /*check*/ void Matrix<T>::get_reg(){
  //select mempool whose size is  same as this matrix
  vector<unique_ptr<T[]> > &mempool = pool_list[pool_index];
  //cout <<"get"<<mempool.size()<<endl;
  //cout<<"index"<<pool_index<<endl;
  //if mempool is not empty, get region from pool
  if(mempool.size()>0){
    val_ptr = move(mempool.back());
    mempool.pop_back();
  }else{
    // if mempool is empty, generate new region and take it
    //cout <<"new"<<endl;
    all_size += nrow*ncol*sizeof(T);
    if(all_size > upper_size){
      //cout << "Error!: too many memory occupied"<<endl;
    }
    val_ptr.reset(new T[nrow*ncol]());
    ncount++;
  }
  //cout <<"get"<<mempool.size()<<endl;
    
}
// return matrix region for pool
template<class T> /*check*/ void Matrix<T>::return_reg(){
  if(reg()){
    //select mempool whose size is  same as this matrix
    vector<unique_ptr<T[] > > &mempool = pool_list[pool_index];
    mempool.push_back(move(val_ptr));
    //cout <<"return"<<mempool.size()<<endl;
    //cout<<"index"<<pool_index<<endl;
  }    
}
//set parameters for size change and pool change
template<class T>
void Matrix<T>::size_set(int n,int m){
  nrow = n;
  ncol = m;
  lda = n;
  lwork = m;
  _size = nrow*ncol;
  if(_size>0){
    //return presize area
    return_reg();
    //pool_index is set
    demmand_pool();
  }
}  
//if size is 0, inform error
template <class T>
void Matrix<T>::size_checker() const{
  if(!(nrow > 0 && ncol > 0)){
    exception_exit("This matrix is unavairable size");
  }
}
template<class T> /*check*/ Matrix<T>::Matrix(){
  size_set(0,0);
}
template<class T> /*check*/ Matrix<T>::Matrix(int n,int m){
  size_set(n,m);
  get_reg();
}
// initialize by arg_mat
template<class T> /*check*/ Matrix<T>::Matrix(const Matrix<T>& arg_mat){
  size_set(0,0);
  //Matrix();
  (*this) = arg_mat;
  //cout <<*oneptr<<endl;
}
template<class T> /*check*/ Matrix<T>::Matrix(Matrix<T>&& arg_mat){
  size_set(0,0);
  //Matrix();
  //copy arg_mat region to this
  (*this) = move(arg_mat);
}
template<class T> /*check*/ Matrix<T>::~Matrix(){
  return_reg();
}
template<class T> /*check*/ T * Matrix<T>::begin() const{
  return(&val_ptr[0]);
}
template<class T> /*check*/ T * Matrix<T>::end() const{
  return(&val_ptr[nrow*ncol]);
}
//if this have region, reg() return true, else false
template<class T> /*check*/ bool Matrix<T>::reg() const{
  bool flag;
  if(val_ptr){
    flag = true;
  }else{
    flag = false;
  }
  return(flag);
}
//  operator Matrix&&()=delete;
//return pointer and this matrix get new
template<class T>
unique_ptr<T[]> Matrix<T>::give_ptr(){
  unique_ptr<T[]> tmp_ptr = move(val_ptr);
    
  return(tmp_ptr);
}
template<class T>
T* Matrix<T>::get() const{
  return(val_ptr.get());
}
template<class T>
void Matrix<T>::resize(int n,int m){
  size_set(n,m);
  get_reg();
}
//if lvalue, copy arg_mat value to this region  
template<class T> /*check*/ void Matrix<T>::operator=(const Matrix<T>& arg_mat){
  //if size is 0, parmeter is reset
  if(_size == 0){
    size_set(arg_mat.rows(),arg_mat.cols());
  }
  if((*this).cols()!=arg_mat.cols() || (*this).rows()!=arg_mat.rows()){
    //different size!
    exception_exit("Diffrent size matrix is substituted");
  }
  //if this dont have region, get new region from mempool
  if(!this->reg()){
    get_reg();
  }
  //copy arg_mat region to this
  copy(arg_mat.begin(), arg_mat.end(), begin());
}
//if rvalue, val_ptr indicate arg_mat ptr
template<class T> /*check*/ void Matrix<T>::operator=(Matrix<T>&& arg_mat){
  //if size is 0, parmeter is reset
  if(_size == 0){
    size_set(arg_mat.rows(),arg_mat.cols());
  }
  if((*this).cols()!=arg_mat.cols() || (*this).rows()!=arg_mat.rows()){
    //different size!
    exception_exit("Diffrent size matrix is substituted");
  }
  //if this have region, return region to pool
  if(this->reg()){
    return_reg();
  }
  //move arg_mat ptr to this    
  val_ptr = arg_mat.give_ptr();
    
}
// (i,j) return refference of Matrix(i,j)
template<class T>
T & Matrix<T>::operator()(int i,int j) const{
  size_checker();
  if(i < 0 || i > nrow-1 || j < 0 || j > ncol -1){
    exception_exit("Out of matrix is selected ");
  }
  return(val_ptr[i+j*nrow]);
}
// (i) return ith element col  first
template<class T>
T & Matrix<T>::operator()(int i) const{
  size_checker();
  int n = i/ncol;
  int m = i%ncol;
  return((*this)(n,m));
}

// Mat + Mat
template<class T>
Matrix<T> Matrix<T>::operator+(Matrix<T> &&arg_mat){
  size_checker();
  if(nrow != arg_mat.rows() || ncol != arg_mat.cols()){
    exception_exit("diffrent size matirces cannot be summed");
  }
  //cout <<"in r"<<endl;
  //cout <<_size<<endl;
  int size = _size;
  double donetmp = done;
  int onetmp = one;
  daxpy_(&size,&donetmp,get(),&onetmp,arg_mat.get(),&onetmp);
  //cout <<arg_mat.str()<<endl;
  // return as raval ref
  return(move(arg_mat));
}
// if arg_mat is lval, make temporaly
template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &arg_mat){
  size_checker();
  //cout <<arg_mat.str()<<endl;
   Matrix<T> arg_tmp(arg_mat);
   //cout <<arg_tmp.str()<<endl;
  //daxpy_(&_size,&done,get(),oneptr.get(),arg_tmp.get(),oneptr.get());
  // cause + switch rval-ref 
  return((*this)+move(arg_tmp));
}
//if single element, this can be summed with T 
template<class T> /*check*/
T Matrix<T>::operator+(T a){
  if(nrow != 1 || ncol != 1){
    exception_exit("This is not single element");
  }
  return((*this)(0,0)+a);
}
// -Mat
template<class T>
Matrix<T> Matrix<T>::operator-() const{
  size_checker();
  Matrix<T> mat(nrow,ncol);
  for(int i = 0;i < nrow;i++){
    for(int j = 0;j < ncol;j++){
      mat(i,j) = - (*this)(i,j);
    }
  }
  return(mat);
}
// Mat - Mat
template<class T>
Matrix<T> Matrix<T>::operator-(Matrix<T> &&arg_mat){
  size_checker();
  return((*this) + -arg_mat);
}

// if arg_mat is lval, make temporaly 
template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &arg_mat){
  size_checker();
  return((*this) + -arg_mat);
}
//if single element, this can be substracted with T 
template<class T> /*check*/
T Matrix<T>::operator-(T a){
  if(nrow != 1 || ncol != 1){
    exception_exit("This is not single element");
  }
  return((*this)(0,0)-a);
}
template<class T> /*check*/ Matrix<T> Matrix<T>::operator*(double a){
  Matrix<T> mat((*this));
  dscal_(&_size, &a, mat.get(), &one);
  return(mat);
}
// Mat * Mat 
template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& arg_mat){
  size_checker();
  if(ncol != arg_mat.rows()){
    exception_exit("In multiplication, size is not fitted");
  }
  Matrix<T> mat(Matrix<T>::Zero(nrow,arg_mat.cols()));
  unique_ptr<int> kptr(new int(arg_mat.cols()));
  char Nstr[] = "N";
  //cout <<"in *"<<endl;
  int nrowtmp = nrow;
  dgemm_(Nstr,Nstr, &nrow, kptr.get(), &ncol,&done,get(),&nrow,arg_mat.get(),&ncol,&done, mat.get(),&nrow);
  return(mat);
}
template<class T>
Matrix<T> Matrix<T>::operator*(Matrix<T>&& arg_mat){
  size_checker();
  return((*this)*arg_mat);
}
//return inverse matrix
template<class T>
Matrix<T> Matrix<T>::inverse(){
  size_checker();
  Matrix<T> mat_tmp((*this));
  // LAPACKのdgetrfサブルーチンを呼んで、行列AをLU分解
  // 引数は全て参照渡し
  int n = nrow;
  int m = ncol;
  dgetrf_( &n, &m, mat_tmp.get(), &lda, ipiv, &info);
  if(info!=0){
    cout <<"error lu"<<endl;
  }
  // LU分解後の行列から逆行列を求める
  // 逆行列は元の配列Aに入る
  dgetri_( &ncol,mat_tmp.get() , &lda, ipiv, work, &lwork, &info);
  if(info!=0){
    cout <<"error inv"<<endl;
  }
  return(mat_tmp);
}
//return transpose matrix of this
template<class T>
Matrix<T> Matrix<T>::transpose() const{
  size_checker();
  Matrix<T> mat_tmp(ncol,nrow);
  for(int i = 0;i < nrow;i++){
    for(int j = 0;j <ncol;j++){
      mat_tmp(j,i) = (*this)(i,j);
    }
  }
  return(mat_tmp);
}
//template<int K>
//  Matrix operator*(Matrix<T>&& arg_mat);
template<class T> /*check*/  Matrix<T> Matrix<T>::Zero(int n,int m){
  Matrix<T> zmat(n,m);
  for(int i = 0;i < n;i++){
    for(int j = 0;j < m;j++){
      zmat(i,j) = 0;
    }
  }
  return(zmat);
}
template<class T> /*check*/ void Matrix<T>::print() const{
  /*
  for(int i = 0;i < N;i++){
    for(int j= 0;j <M;j++){
      cout<<(*this)(i,j)<<"\t";
    }
    cout <<endl;
  }
  */
}
// cast to string ,example cout <<Matrix<<endl;
template<class T> /*check*/
string Matrix<T>::str() const{
  size_checker();
  stringstream ss;
  for(int i = 0;i < nrow;i++){
    for(int j= 0;j <ncol;j++){
      ss.setf(ios::right,ios::adjustfield);
      ss.width(15);
      ss.precision(5);
      ss<<(*this)(i,j);
    }
    ss <<endl;
  }
  return(ss.str());
}
//_size of Matrix
template<class T> /*check*/
int Matrix<T>::size() const{
  return(_size);
}
// number of rows
template<class T> /*check*/  int Matrix<T>::rows() const{
  return(nrow);
}
//number of cols
template<class T> /*check*/  int Matrix<T>::cols() const{
  return(ncol);
}  
//take nth row
template<class T>
Matrix<T> Matrix<T>::row(int n) const{
  size_checker();
  if(n > ncol || n < 0){
    exception_exit("no row");
  }
  Matrix<T> nthrow(1,ncol);
  for(int j = 0;j <ncol;j++){
    nthrow(0,j) = (*this)(n,j);
  }
  return(nthrow);
}
//take mth col
template<class T>
Matrix<T> Matrix<T>::col(int m) const{
  size_checker();
  if(m > ncol || m < 0){
    exception_exit("no column");
  }
  Matrix<T> mthcol(nrow,1);
  for(int i = 0;i <nrow;i++){
    mthcol(i,0) = (*this)(i,m);
  }
  return(mthcol);
}
template<class T>
Matrix<T> Matrix<T>::cut(int i,int j) const{
  size_checker();
  //if large to small rflag on
  bool rflag = false;
  if(i > j){
    rflag = true;
  }
  if(max(i,j) > _size-1 && min(i,j) < 0){
    exception_exit("there is no element of this number");
  }
  Matrix<T> subvec(abs(i-j)+1,1);
  for(int k = 0;k < abs(i-j)+1;k++){
    double val;
    if(rflag){
      val = (*this)(i - k);
    }else{
      val = (*this)(i+k);
    }
    subvec(k) = val;
  }
  return(subvec);
}
template<class T>
Matrix<T> Matrix<T>::col_cut(int i,int j) const{
  size_checker();
  if(max(i,j) > ncol-1 && min(i,j) < 0){
    exception_exit("there is no col of this number");
  }
  //if large to small rflag on
  bool rflag = false;
  if(i > j){
    rflag = true;
  }
  Matrix<T> onecol(nrow,1);
  Matrix<T> submat(nrow,abs(i-j)+ 1);
  for(int k = 0;k < abs(i-j)+1;k++){
    if(rflag){
      onecol = (*this).col(i - k);
    }else{
      onecol = (*this).col(i+k);
    }
    for(int n = 0;n < nrow;n++){
      submat(n,k) = onecol(n);
    }
  }
  return(submat);
}
template<class T>
Matrix<T> Matrix<T>::row_cut(int i,int j) const{
  size_checker();
  if(max(i,j) > nrow-1 && min(i,j) < 0){
    exception_exit("there is no row of this number");
  }
  //if large to small rflag on
  bool rflag = false;
  if(i > j){
    rflag = true;
  }
  Matrix<T> onerow(1,ncol);
  Matrix<T> submat(abs(i-j)+ 1,ncol);
  for(int k = 0;k < abs(i-j)+1;k++){
    if(rflag){
      onerow = (*this).row(i - k);
    }else{
      onerow = (*this).row(i+k);
    }
    for(int m = 0;m < ncol;m++){
      submat(k,m) = onerow(0,m);
    }
  }
  return(submat);
}
    
template<class T>
Matrix<T>::operator T(){
  size_checker();
  if(nrow != 1 || ncol != 1){
    exception_exit("This is not single element");
  }
  return((*this)(0,0));
}

template<class T>
vector<vector<unique_ptr<T[]> > > Matrix<T>::pool_list;
template<class T>
vector<int> Matrix<T>::pool_size_list;
template<class T>
Matrix<T> operator*(T a,Matrix<T>& arg_mat){
  return(arg_mat*a);
}
template<class T>
T operator+(T a,Matrix<T>& arg_mat){
  return(arg_mat+a);
}
template<class T>
Matrix<T> operator*(T a,Matrix<T>&& arg_mat){
  return(arg_mat*a);
}
template<class T>
T operator+(T a,Matrix<T>&& arg_mat){
  return(arg_mat+a);
}
/*
template<class T>
void T::operator=(Matrix<T>& arg_mat){
  if(arg_mat.rows() != 1 || arg_mat.cols() != 1){
    exception_exit("This is not single element");
  }
  (*this) = arg_mat;
}
*/
template<class T>
int Matrix<T>::ipiv[10000];
template<class T>
T Matrix<T>::work[10000];
#endif
