#include <iostream>
#include <vector>
#include "Matrix.h"
using namespace std;
int main(){
  Matrix<double> mat = Matrix<double>::Zero(3,3);
  //cout<<mat.str()<<endl;
  Matrix<double> vec(3,1);
  for(int i = 0;i <3;i++){
    for(int j = 0;j < 3;j++){
      mat(i,j) = i - j;
    }
    vec(i) = i;
  }
  cout<<mat.str();
  cout<<mat.row_cut(1,0).str();  
  for(int i = 0;i < 3;i++){
    //vec = vec + vec;
  }
  mat = mat*mat*mat;
  //mat = vec;
  //mat(10,2);
  //cout<<mat.str()<<endl;
  //cout<<mat.str()<<endl;
  //cout<<vec.str()<<endl;
  double a = (a*vec.transpose()*vec*a); 
  //cout<<a<<endl;
  //cout<<mat.str()<<endl;
  /*
    Vector<double,3> vec;
    for(int i = 0;i <3;i++){
    vec(i)=i;
    }
    mat = vec;
    vec = mat;
    Vector<double,3> vec2(vec);
    //double v = vec*vec;
    string vstr = vec.str();
    string mstr = mat.str();
    cout <<vec2.str()<<endl;
    cout <<mstr<<endl;
  */
  return(0);
}
    
