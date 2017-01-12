#include<Matrix.h>
#include<transition_matrix.h>
#include<random>
double no_minus(double val){
  if(val > 0){
    return(val);
  }else{
    return(0);
  }
}  
template <class T>
Matrix<T> multi_each(Matrix<T> mat1,Matrix<T> mat2){
  Matrix<T> ans = Matrix<T>::Zero(mat1.rows(),mat1.cols());
  for(int i = 0;i < mat1.rows();i++){
    for(int j = 0;j < mat1.cols();j++){
      ans(i,j) = mat1(i,j)*mat2(i,j);
    }
  }
  return(ans);
}
template <class T>
Matrix<T> no_minus(Matrix<T> mat){
  Matrix<T> ans = Matrix<T>::Zero(mat.rows(),mat.cols());
  for(int i = 0;i < mat.rows();i++){
    for(int j = 0;j < mat.cols();j++){
      if(mat(i,j) > 0){
	ans(i,j) = mat(i,j);
      }else{
	ans(i,j) = 0;
      }
    }
  }
  return(ans);
}
template <class T>
Matrix<T> element_each(T el,int N,int M){
  Matrix<T> ans = Matrix<T>::Zero(N,M);
  for(int i = 0;i < N;i++){
    for(int j = 0;j < M;j++){
      ans(i,j) = el;
    }
  }
  return(ans);
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
template <class T>
Matrix<T> CTM::diag(Matrix<T> vec){
  Matrix<T> ans = Matrix<T>::Zero(vec.rows(),vec.rows());
  for(int i = 0;i < vec.rows();i++){
    ans(i,i) = vec(i,0);
  }
  return(ans);
}
double CTM::to_xi(int i){
  double h = (double)1/MSTNUM;
  return(h*i);
}
Vec_st CTM::exp_eg(double t){
  Vec_st exp_eigens(STNUM,1);
  //time setting
  for(int i = 0;i < STNUM;i++){
    exp_eigens(i) = exp(t*ld(i));
  }
  return(exp_eigens);
}
void CTM::pabx_refresh(){
  pabx.clear();
  //cout<<"in"<<endl;
  for(int d = 0;d < dsize;d++){
    Vec_st qx(STNUM,1);
    for(int i = 0;i < STNUM;i++){
      qx(i) = binom(to_xi(i),apl[d],bpl[d]);
    }
    //cout <<apl[d]<<" "<<bpl[d]<<endl;
    //cout<<qx.str()<<endl;
    pabx.push_back(move(qx));
  }
}
//calucu etr
Mat_st CTM::etr_make(double t){
  //cout<<exp(t*ldmx)*epmx<<endl;
  Vec_st exp_eigens = exp_eg(t);
  Mat_st er_mat = element_each<double>(exp(t*ldmx)*epmx,STNUM,STNUM);
  Mat_st etq = u_mat*diag(exp_eigens)*ui_mat;
  //error correct from biased zooming dg
  etq = no_minus(etq - er_mat);
  Mat_st etr = diag(isqdg)*etq*diag(sqdg);
  return(etr);
}
//caluculate data wise 
double CTM::lpe_refresh(vector<int> _apl,vector<int> _bpl,vector<double> _tpl,Vec_st pz){
  apl = _apl;
  bpl = _bpl;
  tpl = _tpl;
  Vec_st tone_vec = element_each((double)1,1,STNUM);
  //Mat_st er_mat =  element_each<double>(epmx,STNUM,STNUM);
  //case
  dsize = tpl.size();
  //calucurate p(alpha,beta|endpoint)
  pabx_refresh();
  //caluculate  etr
  vector<Mat_st> etr(dsize);
  for(int d = 0;d < dsize-1;d++){
    double t = tpl[d+1]-tpl[d];
    etr[d] = etr_make(t);
    //cout<<etr[d].str()<<endl;
    //cout<<t<<endl;
    //cout<<etr[d].str()<<endl;

  }
  //forward
  fwd.resize(dsize);
  fwd[0] = multi_each(pabx[0],pz);
  //cout<<fwd[0].str()<<endl;
  for(int d = 0;d < dsize-1;d++){
    Vec_st nfwd = etr[d]*fwd[d];
    //cout<<etr[d].str()<<endl;
    //cout<<pabx[d+1].str()<<endl;
    fwd[d+1] = multi_each(pabx[d+1],nfwd);
  }
  //cout<<fwd[dsize-1].str()<<endl;
  //forward
  bwd.resize(dsize);
  bwd[dsize-1] = pabx[dsize-1];
  for(int d = 1;d < dsize;d++){
    int rd = dsize - d;
    Vec_st pbwd = (etr[rd-1].transpose())*bwd[rd];
    bwd[rd-1] = multi_each(pabx[rd-1],pbwd);
  }
  if(dsize == 1){
    pe = 1;
  }else{
    //get pe form onevec*fwd(D) into log
    pe = element_each((double)1,1,STNUM)*fwd[dsize-1];
    nf_refresh();
  }
  //cout<<pe<<endl;
  return(log(pe));
}
void CTM::nf_refresh(){
  //modified vec for j and for i
  fd=0;
  ndp=0;
  ndm=0;
  //epmx = 0;
  /*
  double umx = 0,uimx = 0;
  for(int i = 0;i < STNUM;i++){
    for(int j = 0;j < STNUM;j++){
      umx = max(fabs(u_mat(i,j)),umx);
      uimx = max(fabs(ui_mat(i,j)),uimx);
    }
  }
  double u3mx = pow(max(umx,uimx),3);
  double rmmx = 0,rpmx  = 0,fimx = 0;
  for(int i = 0;i < MSTNUM;i++){
    rmmx = max(fabs(r_m(i)),rmmx);
    rpmx = max(fabs(r_p(i)),rpmx);
    fimx = max(fabs(f_i(i)),fimx);
  }
  */
  for(int d = 0;d <dsize-1;d++){
    double t = tpl[d+1] - tpl[d];
    //Mat_st er_mat = element_each<double>(exp(t*ldmx)*epmx,STNUM,STNUM);
    Vec_st exp_eigens = exp_eg(t);
    //caluculate l_vec
    Vec_st r_vec = multi_each(sqdg,fwd[d]);
    Vec_st l_vec = multi_each(bwd[d+1],isqdg);
    //cout<<pabx[d].str()<<endl;
    //cout<<"r_vec : "<<r_vec.str()<<endl;
    //caluculate mat vu_mat
    //caluculate u_eu*ui_ui*u_jv*u_vs*k
    Mat_st tkp = Mat_st::Zero(STNUM,STNUM);
    Mat_st tkm = Mat_st::Zero(STNUM,STNUM);
    Mat_st tkc = Mat_st::Zero(STNUM,STNUM);
    if(t != 0){ // if t = 0, utk is not added
      Mat_st uv_mat(STNUM,STNUM);
      for(int u = 0;u < STNUM;u++){
	for(int v = 0;v < STNUM;v++){
	  double k;
	  if(ld(u)- ld(v) < 10e-200){///MEMO
	    k = exp_eigens(u);
	  }else{
	    k = (exp_eigens(u) - exp_eigens(v))/(t*(ld(u) - ld(v)));
	  }
	  tkp(u,v) = t*pi_mat(u,v)*k;
	  tkm(u,v) = t*mi_mat(u,v)*k;
	  tkc(u,v) = t*ci_mat(u,v)*k;
	}
      }
      fd += l_vec.transpose()*u_mat*tkc*ui_mat*r_vec/pe;
      ndp += l_vec.transpose()*u_mat*tkp*ui_mat*r_vec/pe;
      ndm += l_vec.transpose()*u_mat*tkm*ui_mat*r_vec/pe;
      //sum t*vumat 
    }
  }
  //cout<<"utkrow1 : "<<utk.row(3).str()<<endl;
}
double CTM::safe_lpe_refresh(vector<int> _apl,vector<int> _bpl,vector<double> _tpl,Vec_st pz){
  apl = _apl;
  bpl = _bpl;
  tpl = _tpl;
  Vec_st tone_vec = element_each((double)1,1,STNUM);
  //Mat_st er_mat =  element_each<double>(epmx,STNUM,STNUM);
  //case
  dsize = tpl.size();
  //calucurate p(alpha,beta|endpoint)
  pabx_refresh();
  //caluculate  etr
  vector<Mat_st> edtr(dsize,Mat_st::Zero(STNUM,STNUM));
  for(int d = 0;d < dsize-1;d++){
    double t = tpl[d+1]-tpl[d];
    double dt = t/safe_disc;
    edtr[d] = etr_make(dt);
  }
  //forward
  fwd.assign((dsize-1)*safe_disc+1,Mat_st::Zero(STNUM,1));
  fwd[0] = multi_each(pabx[0],pz);
  for(int d = 0;d < dsize-1;d++){
    double dt = (double)1/safe_disc;
    for(int pt = 0;pt < safe_disc;pt++){
      fwd[d*safe_disc + pt + 1] = edtr[d]*fwd[d*safe_disc + pt];
    }
    fwd[(d+1)*safe_disc] = multi_each(pabx[d+1],fwd[(d+1)*safe_disc]);
  }
  //backward
  bwd.assign((dsize-1)*safe_disc+1,Mat_st::Zero(STNUM,1));
  bwd[(dsize-1)*safe_disc] = element_each((double)1,STNUM,1);
  for(int d = 1;d < dsize;d++){
    int rd = dsize - d;
    double dt = (double)1/safe_disc;
    bwd[(rd)*safe_disc-1] = (edtr[rd-1]).transpose()*multi_each(pabx[rd],bwd[(rd)*safe_disc]);
    for(int pt = 1;pt < safe_disc;pt++){
      bwd[rd*safe_disc - pt - 1] = (edtr[rd-1]).transpose()*bwd[rd*safe_disc - pt];
      int loc = rd*safe_disc - pt - 1;
      //cout<<log(bwd[loc].transpose()*fwd[loc])<<endl;
    }
  }
  //get pe form onevec*fwd(D) into log
  if(dsize == 1){
    pe = 1;
  }else{
    pe = (element_each((double)1,1,STNUM)*fwd[(dsize-1)*safe_disc])(0,0);
    safe_nf_refresh();
  }
  //cout<<"ture"<<log(pe)<<endl;
  return(log(pe));
}
void CTM::safe_nf_refresh(){
  //modified vec for j and for i
  fd=0;
  ndp=0;
  ndm=0; 
  //epmx = 0;
  int loc;
  double t;
  double dt = (double)1/safe_disc;
  for(int d = 0;d <dsize-1;d++){
    for(int pt = 0;pt < safe_disc;pt++){
      t = tpl[d+1] - tpl[d];
      loc = d*safe_disc + pt;
      //cout<<loc<<endl;
      ndp += bwd[loc].cut(1,STNUM-1).transpose()*multi_each(r_p,fwd[loc].cut(0,MSTNUM-1))*dt*t/pe;
      ndm += bwd[loc].cut(0,MSTNUM-1).transpose()*multi_each(r_m,fwd[loc].cut(1,STNUM-1))*dt*t/pe;
      fd += bwd[loc].transpose()*multi_each(fiv,fwd[loc])*dt*t/pe;
    }
  }
  t = tpl[dsize-1] - tpl[dsize-2];
  loc = (dsize-1)*safe_disc;
  //cout<<bwd[loc].size()<<endl;
  //cout<<fwd[loc].size()<<endl;
  ndp += bwd[loc].cut(1,STNUM-1).transpose()*multi_each(r_p,fwd[loc].cut(0,MSTNUM-1))*dt*t/pe;
  ndm += bwd[loc].cut(0,MSTNUM-1).transpose()*multi_each(r_m,fwd[loc].cut(1,STNUM-1))*dt*t/pe;
  fd += bwd[loc].transpose()*multi_each(fiv,fwd[loc])*dt*t/pe;
}

double CTM::n_pget(int i){
  double ans = ndp;
  return(ans);
}
double CTM::n_mget(int i){
  double ans = ndm;
  return(ans);
}
double CTM::tf_get(int i){
  double ans = fd;
  return(ans);
}
