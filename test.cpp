#include <iostream>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

//[[Rcpp::depends(RcppArmadillo)]]

//The whole parameters/variables are {$r,z,c,b,v,m,k,\mu,\sigma,\lambda,\rho,\gamma$}
//functions to calculate harmonic number;
//[[Rcpp::export]]
double har(int n){
  if(n==1)
    return 1;
  double result = har(n-1)+1.0/n;
  return result;
}

//functions to compute Khatri-Rao product;
//[[Rcpp::export]]
mat KR(mat A,mat B){
  int ncol = A.n_cols;
  int nrowa = A.n_rows;
  int nrowb = B.n_rows; 
  mat result(nrowa*nrowb,ncol);
  for(int i=0;i<ncol;i++)
    result.col(i) = kron(A.col(i),B.col(i));
  return result;
}


//functions to calculate pz;
//[[Rcpp::export]]
double log_pz(vec a,mat b1,mat b2,vec c,vec d1,vec d2){
  vec t1 = b1*a+d1;
  vec t2 = b2*a+d2;
  double result = sum((t1%(c==-1)+t2%(c==1))-log(exp(t1)+1+exp(t2)));
  return result;
}
//[[Rcpp::export]]
double log_pz3(vec a1,vec a2,mat b,vec c,vec d1,vec d2){
  vec t1 = b*a1+d1;
  vec t2 = b*a2+d2;
  double result = sum((t1%(c==-1)+t2%(c==1))-log(exp(t1)+1+exp(t2)));
  return result;
}

//functions to calculate px;
//[[Rcpp::export]]
double log_px(double x,double v1,double v2,double mu,double sigma2,double p_1,double p0,double p1){
  double result=log((1/v1)*p_1*((x<mu)&&(x>mu-v1))+R::dnorm(x,mu,sigma2,false)*p0+log(1/v2)*((x>mu)&&(x<mu+v2))*p1);
  return result;
}

//functions to calculate pz_unit;
//[[Rcpp::export]]
double log_pz_unit(int i0,int i1,int i2,int z,mat c1,mat c2,mat c3,mat l1,mat l2,mat b1, mat b2){
  int r = c1.n_cols;
  double r1=0,r2=0,result=0;
  for(int i=0;i<r;i++){
    r1 += c1(i0,i)*c2(i1,i)*l1(i2,i)*(c3(i2,i)==-1);
    r2 += c1(i0,i)*c2(i1,i)*l2(i2,i)*(c3(i2,i)==1);
  }
  r1 = r1+b1(i1,i2);
  r2 = r2+b2(i1,i2);
  result = r1*(z==-1)+r2*(z==1)-log(exp(r1)+1+exp(r2));
  return result;
}

//functions to calculate pz_b;
//[[Rcpp::export]]
double log_pz_b1(int i1,int i2,double b,vec d,mat c1,mat c2,mat c3,mat l1,mat l2,mat b2,cube z){
  int r = c1.n_cols;
  vec r1(d[0],fill::zeros),r2(d[0],fill::zeros);
  for(int i0=0;i0<d[0];i0++){
    for(int ind=0;ind<r;ind++){
      r1(i0) += c1(i0,ind)*c2(i1,ind)*l1(i2,ind)*(c3(i2,ind)==-1);
      r2(i0) +=c1(i0,ind)*c2(i1,ind)*l2(i2,ind)*(c3(i2,ind)==1);
    }
    r1(i0) += b;
    r2(i0) += b2(i1,i2);
  }
  double result;
  vec zv = vectorise(z.subcube(0,i1,i2,d[0]-1,i1,i2));
  result = sum((r1%(zv==-1)+r2%(zv==1))-log(exp(r1)+1+exp(r2)));
  return result;
}

//[[Rcpp::export]]
double log_pz_b2(int i1,int i2,double b,vec d,mat c1,mat c2,mat c3,mat l1,mat l2,mat b1,cube z){
  int r = c1.n_cols;
  vec r1(d[0],fill::zeros),r2(d[0],fill::zeros);
  for(int i0=0;i0<d[0];i0++){
    for(int ind=0;ind<r;ind++){
      r1(i0) += c1(i0,ind)*c2(i1,ind)*l1(i2,ind)*(c3(i2,ind)==-1);
      r2(i0) +=c1(i0,ind)*c2(i1,ind)*l2(i2,ind)*(c3(i2,ind)==1);
    }
    r1(i0) += b1(i1,i2);
    r2(i0) += b;
  }
  double result;
  vec zv = vectorise(z.subcube(0,i1,i2,d[0]-1,i1,i2));
  result = sum((r1%(zv==-1)+r2%(zv==1))-log(exp(r1)+1+exp(r2)));
  return result;
}
//function to update variables and parameters

//function to update c1
//[[Rcpp::export]]
List update_c1(mat c1,mat c2,mat c3,mat lambda1,mat lambda2,vec d,cube z,mat b1,mat b2,vec gamma,double m,double rho,double al,double bl)
{
  List re(5);
  for(int i=0;i<d[0];i++){
    int r = c1.n_cols,r_star = R::rpois(m/d[0]);
    mat m1 = KR(lambda1%(c3==-1),c2);
    mat m2 = KR(lambda2%(c3==1),c2);
    vec index = zeros<vec>(r);
    vec a = vectorise(c1.row(i));
    vec zvec = vectorise(z.row(i));
    for(int j=0;j<r;j++){
      
      double x = sum(c1.col(j))-c1(i,j);
      if((x>0)&&x<(d[0]-1)){
        index(j) = 1;
        a(j) =1;
        double log_p1 = log(x)+log_pz(a,m1,m2,zvec,vectorise(b1.as_col()),vectorise(b2.as_col()));
        a(j) =0;
        double log_p0 = log(d[0]-x)+log_pz(a,m1,m2,zvec,vectorise(b1.as_col()),vectorise(b2.as_col()));
        double p0 = 1/(exp(log_p1-log_p0)+1);
        c1(i,j) = (randu()>p0);a(j) = c1(i,j);
      }
    }
    if((sum(index) != 0) && (sum(index) != r)) {
      c1 = c1.cols(find(index));  c2 = c2.cols(find(index)); c3 = c3.cols(find(index)); lambda1 = lambda1.cols(find(index));lambda2 = lambda2.cols(find(index));
      a = a(find(index));  m1 = m1.cols(find(index));  m2 = m2.cols(find(index));
    }
    
    if(r_star >0){
      //propose new feature
      
      mat c1t(d[0],r_star,fill::zeros),c2t(d[1],r_star,fill::randu),c3t(d[2],r_star,fill::randu),lambda1t=randg<mat>(d[2],r_star, distr_param(al,bl)),lambda2t=randg<mat>(d[2],r_star, distr_param(al,bl));
      c2t = ones(size(c2t))%(c2t<rho);c3t = ones(size(c3t))%(c3t>(1-gamma[2]))-ones(size(c3t))%(c3t<gamma[0]);
      c1t.row(i) = ones<rowvec>(r_star);
      mat c1_new=join_rows(c1,c1t),c2_new=join_rows(c2,c2t),c3_new=join_rows(c3,c3t),l1_new=join_rows(lambda1,lambda1t),l2_new=join_rows(lambda2,lambda2t);
      mat m1_new = KR(l1_new%(c3_new==-1),c2_new),m2_new = KR(l2_new%(c3_new==1),c2_new);
      vec a_new = vectorise(c1_new.row(i));
      
      double log_p = log_pz(a_new,m1_new,m2_new,zvec,vectorise(b1),vectorise(b2))-log_pz(a,m1,m2,zvec,vectorise(b1),vectorise(b2));
      double p = exp(log_p);
      double u = randu();
      if(u<p){
        c1 = c1_new;c2 = c2_new;c3=c3_new;lambda1 = l1_new;lambda2 = l2_new;
      }
    }
  }
  re[0]=c1;re[1]=c2;re[2]=c3;re[3]=lambda1;re[4]=lambda2;
  return re;
}

//function to update c2
//[[Rcpp::export]]
mat update_c2(mat c1,mat c2,mat c3,mat lambda1,mat lambda2,vec d,cube z,mat b1,mat b2,double rho)
{
  int r = c1.n_cols;
  mat m1 = KR(lambda1%(c3==-1),c1);
  mat m2 = KR(lambda2%(c3==1),c1);
  for(int i=0;i<d[1];i++){
    vec a = vectorise(c2.row(i));
    vec zvec = vectorise(z.col(i));
    NumericVector b1r = NumericVector(b1.row(i).begin(), b1.row(i).end());
    NumericVector b2r = NumericVector(b2.row(i).begin(), b2.row(i).end());
    for(int j=0;j<r;j++){
      
      a(j) =1;
      double log_p1 = log(rho)+log_pz(a,m1,m2,zvec,rep_each(b1r,int(d[0])),rep_each(b2r,int(d[0])));
      a(j) =0;
      double log_p0 = log(1-rho)+log_pz(a,m1,m2,zvec,rep_each(b1r,int(d[0])),rep_each(b2r,int(d[0])));
      double p0 = 1/(exp(log_p1-log_p0)+1);
      c2(i,j) = (randu()>p0);a(j) = c2(i,j);
    }
  }
  return c2;
}

//function to update c3 and lambda1,lambda2
//[[Rcpp::export]]
List update_c3(mat c1,mat c2,mat c3,mat lambda1,mat lambda2,vec d,cube z,mat b1,mat b2,vec gamma,double sigmal,double al,double bl)
{
  int r = c1.n_cols;
  List re(3);
  mat m12 = KR(c2,c1);
  for(int i=0;i<d[2];i++){
    vec a = vectorise(c3.row(i));
    vec w1 = vectorise(lambda1.row(i));
    vec w2 = vectorise(lambda2.row(i));
    vec zvec = vectorise(z.slice(i));
    NumericVector b1l = NumericVector(b1.col(i).begin(), b1.col(i).end());
    NumericVector b2l = NumericVector(b2.col(i).begin(), b2.col(i).end());
    for(int j=0;j<r;j++){
      
      a(j) =-1;
      vec a1 = (a==-1)%w1;
      vec a2 = (a==1)%w2;
      double log_p_1 = log(gamma[0])+log_pz3(a1,a2,m12,zvec,rep_each(b1l,d[0]),rep_each(b2l,d[0]));
      a(j) = 0;
      a1 = (a==-1)%w1;
      a2 = (a==1)%w2;
      double log_p0 = log(gamma[1])+log_pz3(a1,a2,m12,zvec,rep_each(b1l,d[0]),rep_each(b2l,d[0]));
      a(j) = 1;
      a1 = (a==-1)%w1;
      a2 = (a==1)%w2;
      double log_p1 = log(gamma[2])+log_pz3(a1,a2,m12,zvec,rep_each(b1l,d[0]),rep_each(b2l,d[0]));
      vec log_pr={log_p_1,log_p0,log_p1},pr = exp(log_pr-max(log_pr));pr = pr/sum(pr);
      double u = randu();
      c3(i,j)=0;
      if(u>(1-pr(2))) c3(i,j)=1;
      if(u<pr(0)) c3(i,j)=-1;
      a(j) = c3(i,j);
      
      double temp1;
      do temp1 = sigmal*randn()+lambda1(i,j); while(temp1<0);
      vec tempv = w1;
      tempv(j) = temp1;
      vec a1t = tempv%vectorise(c3.row(i)==-1);
      //vec a2t = w2%vectorise(c3.row(i)==1);
      a1 = w1%vectorise(c3.row(i)==-1);
      a2 = w2%vectorise(c3.row(i)==1);
      double log_p = R::dgamma(temp1,al,bl,true)+log_pz3(a1t,a2,m12,zvec,rep_each(b1l,d[0]),rep_each(b2l,d[0]))-R::dgamma(lambda1(i,j),al,bl,true)-log_pz3(a1,a2,m12,zvec,rep(b1l,d[0]),rep(b2l,d[0]));
      double p = exp(log_p);
      u = randu();
      if(u<p){
        lambda1(i,j)=temp1;w1 = tempv;
      }
        
      
      double temp2;
      do temp2 = sigmal*randn()+lambda2(i,j); while(temp2<0);
      tempv = w2;
      tempv(j) = temp2;
      //a1t = tempv%vectorise(c3.row(i)==-1);
      vec a2t = tempv%vectorise(c3.row(i)==1);
      a1 = w1%vectorise(c3.row(i)==-1);
      a2 = w2%vectorise(c3.row(i)==1);
      log_p = R::dgamma(temp2,al,bl,true)+log_pz3(a1,a2t,m12,zvec,rep_each(b1l,d[0]),rep_each(b2l,d[0]))-R::dgamma(lambda2(i,j),al,bl,true)-log_pz3(a1,a2,m12,zvec,rep(b1l,d[0]),rep(b2l,d[0]));
      p = exp(log_p);
      u = randu();
      if(u<p){
        lambda2(i,j)=temp2;w2 = tempv;
      }
    }
  }
  re[0]=c3;re[1]=lambda1;re[2]=lambda2;
  return re;
}


//function to update z and v
//[[Rcpp::export]]
List update_zv(mat c1,mat c2,mat c3,mat lambda1,mat lambda2,vec d,cube z,mat b1,mat b2,mat v1,mat v2,double av,double bv,cube y,mat mu,mat sigma2)
{
  //int r = c1.n_cols;
  double k0=5;
  List re(3);
  for(int i2=0;i2<d[2];i2++){
    for(int i1=0;i1<d[1];i1++){
      vec ztemp(d[0],fill::zeros);
      double temp1 = randg( distr_param(av,bv) );
      double temp2 = randg( distr_param(av,bv) );
      if(temp1 < sqrt(sigma2(i1,i2))*k0) temp1=sqrt(sigma2(i1,i2))*k0;
      if(temp2 < sqrt(sigma2(i1,i2))*k0) temp2=sqrt(sigma2(i1,i2))*k0;
      double log_p=0;
      for(int i0=0;i0<d[0];i0++){
        //double log_p_1 = log_pz_unit(i0,i1,i2,-1,c1,c2,c3,lambda1,lambda2,b1,b2)+log_px(y(i0,i1,i2),-1,v1(i1,i2),v2(i1,i2),mu(i1,i2),sigma2(i1,i2));
        //double log_p0 = log_pz_unit(i0,i1,i2,0,c1,c2,c3,lambda1,lambda2,b1,b2)+log_px(y(i0,i1,i2),0,v1(i1,i2),v2(i1,i2),mu(i1,i2),sigma2(i1,i2));
        //double log_p1 = log_pz_unit(i0,i1,i2,1,c1,c2,c3,lambda1,lambda2,b1,b2)+log_px(y(i0,i1,i2),1,v1(i1,i2),v2(i1,i2),mu(i1,i2),sigma2(i1,i2));
        //double p_1 = exp(log_p_1), p0 = exp(log_p0),p1 = exp(log_p1);
        //double u = randu();
        //ztemp(i0) = (u>(1-p1/(p_1+p0+p1)))-(u<(p_1/(p_1+p0+p1)));
        double t0 = log_pz_unit(i0,i1,i2,0,c1,c2,c3,lambda1,lambda2,b1,b2);
        double t_1 = log_pz_unit(i0,i1,i2,-1,c1,c2,c3,lambda1,lambda2,b1,b2);
        double t1 = log_pz_unit(i0,i1,i2,1,c1,c2,c3,lambda1,lambda2,b1,b2);
        double log_p0 = t0+R::dnorm(y(i0,i1,i2),mu(i1,i2),sqrt(sigma2(i1,i2)),true);
        double log_p1 = t1-log(temp2);
        double log_p_1 = t_1-log(temp1);
        if((y(i0,i1,i2)>=mu(i1,i2))  && (y(i0,i1,i2) < (mu(i1,i2) + temp2)))  ztemp(i0) = ((log_p1 - log_p0) > log(1 / randu() - 1));
        if((y(i0,i1,i2)<mu(i1,i2)) && (y(i0,i1,i2) > (mu(i1,i2) - temp1)))  ztemp(i0) = -((log_p_1 - log_p0) > log(1 / randu() - 1));
        log_p += log_px(y(i0,i1,i2),temp1,temp2,mu(i1,i2),sigma2(i1,i2),exp(t_1),exp(t0),exp(t1))-log_px(y(i0,i1,i2),v1(i1,i2),v2(i1,i2),mu(i1,i2),sigma2(i1,i2),exp(t_1),exp(t0),exp(t1));
      }
      //z.subcube(0,i1,i2,d[0]-1,i1,i2) = ztemp;
      double p = exp(log_p),u = randu();
      if(u<p){
        v1(i1,i2) = temp1;v2(i1,i2) = temp2;z.subcube(0,i1,i2,d[0]-1,i1,i2) = ztemp;
      }
    }
  }
  re[0]=z;re[1]=v1;re[2]=v2;
  return re;
}

//function to update mu and sigma2
//[[Rcpp::export]]
mat update_norm1(cube z,cube y,mat mu,mat v1,mat v2)
{
  double as=1e-3,bs=1e-3,k0=5;
  int nrow = mu.n_rows,ncol = mu.n_cols,fir = z.n_rows;mat result(nrow,ncol);
  for(int i =0;i<nrow;i++){
    for(int j =0;j<ncol;j++){
      vec uk={v1(i,j) / k0,  v2(i,j) / k0};
      vec zt = vectorise(z.subcube(0,i,j,fir-1,i,j)),con = ones(fir),yt = vectorise(y.subcube(0,i,j,fir-1,i,j));
      result(i,j) = 1 / randg(distr_param(as+0.5*sum(con%(zt==0)), 1 / (bs+0.5*sum((zt==0)%pow(yt-mu(i,j),2)))));
      if(result(i,j) > pow(min(uk), 2)) result(i,j) = pow(min(uk), 2);
    }
  }
  return result;   
}

//[[Rcpp::export]]
mat update_norm2(cube z,cube y,double sigmamu2,mat sigma2)
{
  int nrow = sigma2.n_rows,ncol = sigma2.n_cols,fir = z.n_rows;mat result(nrow,ncol);
  for(int i =0;i<nrow;i++){
    for(int j =0;j<ncol;j++){
      vec zt = vectorise(z.subcube(0,i,j,fir-1,i,j)),con = ones(fir),yt = vectorise(y.subcube(0,i,j,fir-1,i,j));
      double st = 1/(1/sigmamu2+sum(con%(zt==0))/sigma2(i,j));
      result(i,j) = sqrt(st)*randn()+st*sum((zt==0)%yt)/sigma2(i,j);
    }
  }
  return result;   
}
//function to update b1 and b2
//[[Rcpp::export]]
List update_b(vec d,double mub,double sigmab,mat &b1,mat &b2,mat c1,mat c2,mat c3,mat l1,mat l2,cube z)
{
  //int r = c1.n_cols;
  List re(2);
  for(int i2=0;i2<d[2];i2++){
    for(int i1=0;i1<d[1];i1++){
      double temp1 = randn()+b1(i1,i2);
      double log_p = R::dnorm(temp1,mub,sigmab,true)+log_pz_b1(i1,i2,temp1,d,c1,c2,c3,l1,l2,b2,z)-R::dnorm(b1(i1,i2),mub,sigmab,true)-log_pz_b1(i1,i2,b1(i1,i2),d,c1,c2,c3,l1,l2,b2,z);
      double p = exp(log_p);
      double u = randu();
      if(u<p)
        b1(i1,i2)=temp1;
      
      double temp2 = randn()+b2(i1,i2);
      log_p = R::dnorm(temp2,mub,sigmab,true)+log_pz_b2(i1,i2,temp2,d,c1,c2,c3,l1,l2,b1,z)-R::dnorm(b2(i1,i2),mub,sigmab,true)-log_pz_b1(i1,i2,b2(i1,i2),d,c1,c2,c3,l1,l2,b1,z);
      p = exp(log_p);
      u = randu();
      if(u<p)
        b2(i1,i2)=temp2;
    }
  }
  re[0]=b1;re[1]=b2;
  return re;
}
