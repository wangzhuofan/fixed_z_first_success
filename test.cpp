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

//function to generate dirichlet random variable;
//[[Rcpp::export]]
mat rdirichlet_cpp(int num_samples,vec alpha_m) {
  int distribution_size = alpha_m.n_elem;
  // each row will be a draw from a Dirichlet
  arma::mat distribution = arma::zeros(num_samples, distribution_size);
  
  for (int i = 0; i < num_samples; ++i) {
    double sum_term = 0;
    // loop through the distribution and draw Gamma variables
    for (int j = 0; j < distribution_size; ++j) {
      double cur = randg(distr_param(alpha_m[j],1.0));
      distribution(i,j) = cur;
      sum_term += cur;
    }
    // now normalize
    for (int j = 0; j < distribution_size; ++j) {
      distribution(i,j) = distribution(i,j)/sum_term;
    }
  }
  return(distribution);
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
double log_px(double x,double z,double v1,double v2,double mu,double sigma2){
  double result=-log(v1)*(z==-1)+R::dnorm(x,mu,sigma2,true)*(z==0)-log(v2)*(z==1);
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
      if(x>0){
        index(j) = 1;
        a(j) =1;
        double log_p1 = log(x)+log_pz(a,m1,m2,zvec,vectorise(b1.as_col()),vectorise(b2.as_col()));
        a(j) =0;
        double log_p0 = log(d[0]-x)+log_pz(a,m1,m2,zvec,vectorise(b1.as_col()),vectorise(b2.as_col()));
        double p0 = 1/(exp(log_p1-log_p0)+1);
        c1(i,j) = (randu()>p0);a(j) = c1(i,j);
      }
    }
    if((sum(index) != 0) & (sum(index) != r)) {
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

//function to update class r and new feature

//List new_feature(mat c1,mat c2,mat c3,mat lambda1,mat lambda2,vec d,cube z,mat b1,mat b2,vec gamma,double m,double rho,double al,double bl)
//{
  //int r = c1.n_cols;
  //List re(5);
  //int r_star = R::rpois(m/d[0]);
  //if(r_star >0){
    //propose new feature
    
    //mat c1t(d[0],r_star,fill::ones),c2t(d[1],r_star,fill::randu),c3t(d[2],r_star,fill::randu),lambda1t=randg<mat>(d[2],r_star, distr_param(al,bl)),lambda2t=randg<mat>(d[2],r_star, distr_param(al,bl));
    //c2t = ones(size(c2t))%(c2t<rho);c3t = ones(size(c3t))%(c3t>(1-gamma[2]))-ones(size(c3t))%(c3t<gamma[0]);
    //mat c1_new=join_rows(c1,c1t),c2_new=join_rows(c2,c2t),c3_new=join_rows(c3,c3t),l1_new=join_rows(lambda1,lambda1t),l2_new=join_rows(lambda2,lambda2t);
    //mat m1_old = KR(lambda1%(c3==-1),c2),m2_old = KR(lambda2%(c3==1),c2),m1_new = KR(l1_new%(c3_new==-1),c2_new),m2_new = KR(l2_new%(c3_new==1),c2_new);
    
    //double flag =0;
    //for(int i=0;i<d[0];i++){
      //vec a_old = vectorise(c1.row(i));
      //vec a_new = vectorise(c1_new.row(i));
      //vec zvec = vectorise(z.row(i));
      //double log_p = log_pz(a_new,m1_new,m2_new,zvec,vectorise(b1),vectorise(b2))-log_pz(a_old,m1_old,m2_old,zvec,vectorise(b1),vectorise(b2));
      //double p = exp(log_p);
      //double u = randu();
      //if(u>p){
        //flag +=1;
        //c1_new.submat(i,r,i,r_star+r-1) = zeros<mat>(1,r_star);
      //}
    //}
    //if(flag<d[0]){
      //c1 = c1_new;c2 = c2_new;c3 = c3_new;lambda1 = l1_new;lambda2 = l2_new;
    //}
  //}
  //re[0]=c1;re[1]=c2;re[2]=c3;re[3]=lambda1;re[4]=lambda2;
  //return re;
//}

//function to update z and v
//[[Rcpp::export]]
List update_zv(mat c1,mat c2,mat c3,mat lambda1,mat lambda2,vec d,cube z,mat b1,mat b2,mat v1,mat v2,double av,double bv,cube y,mat mu,mat sigma2)
{
  int r = c1.n_cols;double k0=5;
  List re(3);
  for(int i2=0;i2<d[2];i2++)
    for(int i1=0;i1<d[1];i1++){
      vec ztemp(d[0]);
      double temp1 = randg( distr_param(av,bv) );
      double temp2 = randg( distr_param(av,bv) );
      if(temp1 < sqrt(sigma2(i1,i2))*k0) temp1=sqrt(sigma2(i1,i2))*k0;
      if(temp2 < sqrt(sigma2(i1,i2))*k0) temp2=sqrt(sigma2(i1,i2))*k0;
      double log_p=0;
      for(int i0=0;i0<d[0];i0++){
        double log_p_1 = log_pz_unit(i0,i1,i2,-1,c1,c2,c3,lambda1,lambda2,b1,b2)+log_px(y(i0,i1,i2),-1,v1(i1,i2),v2(i1,i2),mu(i1,i2),sigma2(i1,i2));
        double log_p0 = log_pz_unit(i0,i1,i2,0,c1,c2,c3,lambda1,lambda2,b1,b2)+log_px(y(i0,i1,i2),0,v1(i1,i2),v2(i1,i2),mu(i1,i2),sigma2(i1,i2));
        double log_p1 = log_pz_unit(i0,i1,i2,1,c1,c2,c3,lambda1,lambda2,b1,b2)+log_px(y(i0,i1,i2),1,v1(i1,i2),v2(i1,i2),mu(i1,i2),sigma2(i1,i2));
        double p_1 = exp(log_p_1), p0 = exp(log_p0),p1 = exp(log_p1);
        double u = randu();
        ztemp(i0) = (u>(1-p1/(p_1+p0+p1)))-(u<(p_1/(p_1+p0+p1)));
        log_p += log_px(y(i0,i1,i2),z(i0,i1,i2),v1(i1,i2),v2(i1,i2),mu(i1,i2),sigma2(i1,i2))-log_px(y(i0,i1,i2),z(i0,i1,i2),temp1,temp2,mu(i1,i2),sigma2(i1,i2));
      }
      double p = exp(-log_p);
      double u = randu();
      if(u<p){
        v1(i1,i2) = temp1;v2(i1,i2) = temp2;z.subcube(0,i1,i2,d[0]-1,i1,i2) = ztemp;
      }
    }
    re[0]=z;re[1]=v1;re[2]=v2;
    return re;
}

//function to update mu and sigma2
//[[Rcpp::export]]
mat update_norm(cube z,cube y,double as,double bs,double sigmamu,mat mu,double k0,mat v1,mat v2)
{
  as=1e-3;bs=1e-3;k0=5;
  int nrow = mu.n_rows,ncol = mu.n_cols,fir = z.n_rows;mat result(2*nrow,ncol);
  for(int i =0;i<nrow;i++)
    for(int j =0;j<ncol;j++){
      vec uk={v1(i,j) / k0,  v2(i,j) / k0};
      vec zt = vectorise(z.subcube(0,i,j,fir-1,i,j)),con = ones(fir),yt = vectorise(y.subcube(0,i,j,fir-1,i,j));
      result(i,j) = 1 / randg(distr_param(as+0.5*sum(con%(zt==0)), 1 / (bs+0.5*sum((zt==0)%pow(yt-mu(i,j),2)))));
      if(result(i,j) > pow(min(uk), 2)) result(i,j) = pow(min(uk), 2);
      double st = 1/(1/sigmamu+sum(con%(zt==0))/result(i,j));
      result(nrow+i,j) = sqrt(st)*randn()+st*sum((zt==0)%yt)/result(i,j);
    }
   return result;   
}
//function to update b1 and b2
//[[Rcpp::export]]
List update_b(vec d,double mub,double sigmab,mat &b1,mat &b2,mat c1,mat c2,mat c3,mat l1,mat l2,cube z)
{
  int r = c1.n_cols;
  List re(2);
  for(int i2=0;i2<d[2];i2++)
    for(int i1=0;i1<d[1];i1++){
      double temp1 = randn()*sigmab+b1(i1,i2);
      double log_p = R::dnorm(temp1,mub,sigmab,true)+log_pz_b1(i1,i2,temp1,d,c1,c2,c3,l1,l2,b2,z)-R::dnorm(b1(i1,i2),mub,sigmab,true)-log_pz_b1(i1,i2,b1(i1,i2),d,c1,c2,c3,l1,l2,b2,z);
      double p = exp(log_p);
      double u = randu();
      if(u<p)
        b1(i1,i2)=temp1;
      
      double temp2 = randn()*sigmab+b2(i1,i2);
      log_p = R::dnorm(temp2,mub,sigmab,true)+log_pz_b2(i1,i2,temp2,d,c1,c2,c3,l1,l2,b1,z)-R::dnorm(b2(i1,i2),mub,sigmab,true)-log_pz_b1(i1,i2,b2(i1,i2),d,c1,c2,c3,l1,l2,b1,z);
      p = exp(log_p);
      u = randu();
      if(u<p)
        b2(i1,i2)=temp2;
    }
    re[0]=b1;re[1]=b2;
    return re;
}


//List Multi_way_MCMC(cube y,int iters){
  //List matrix_gibbs(iters);
  //initialization
  //fixed parameters
  //double as=1,bs=1,sigmamu=1,al=1,bl=1,av=1,bv=1,am=1, bm=1,sigmal=0.1,mub=0,sigmab=0.1,arho=1,brho=1,psi_1=1/3,psi0=1/3,psi1=1/3;
  //dimension
  //vec d(3);
  //d[0] = y.n_rows;
  //d[1] = y.n_cols;
  //d[2] = y.n_slices;
  //class number
  //int r=1;
  //observation and latent
  //cube z=randi<cube>(d[0],d[1],d[2],distr_param(-1,1));
  //class
  //mat c1 = randi<mat>(d[0], r, distr_param(0,1));
  //mat c2 = randi<mat>(d[1], r, distr_param(0,1));
  //mat c3 = randi<mat>(d[2], r, distr_param(-1,1));
  //mat lambda1(d[2],r,fill::randn);
  //mat lambda2(d[2],r,fill::randn);
  
  //paramters
  //mat mu(d[1],d[2],fill::zeros);
  //mat v1(d[1],d[2],fill::randu);
  //mat v2(d[1],d[2],fill::randu);
  //mat sigma2(d[1],d[2],fill::ones);
  //mat b1(d[1],d[2],fill::randn);
  //mat b2(d[1],d[2],fill::randn);
  //double m=1,rho=0.2;
  //vec gamma={0.1,0.8,0.1};
  //MCMC update
  //for(int it=0;it<iters;it++){
    //update $c_1$ and delete  zero-column and related feature
    //List up = update_c1(c1,c2,c3,lambda1,lambda2,d,z,b1,b2);
    //update $c_2$
    //mat c1 = up[0],c2=up[1],c3=up[2],lambda1 = up[3],lambda2 = up[4];r = c1.n_cols;
    //c2 = update_c2(c1,c2,c3,lambda1,lambda2,d,z,b1,b2,rho);
    //List up2 = update_c3(c1,c2,c3,lambda1,lambda2,d,z,b1,b2,gamma,sigmal,al,bl);
    //mat t0=up2[0],t1=up2[1],t2=up2[2];
    //c3 = t0;lambda1 = t1;lambda2 = t2;
    //List up3 = new_feature(c1,c2,c3,lambda1,lambda2,d,z,b1,b2,gamma, m,rho,al,bl);
    //mat w0 = up3[0],w1=up3[1],w2=up3[2],w3 =up3[3],w4=up3[4];
    //c1=w0;c2=w1;c3=w2;lambda1=w3;lambda2=w4;r=c1.n_cols;
    //List up4 =update_zv(c1,c2,c3,lambda1,lambda2,d,z,b1,b2,v1,v2,av,bv,y,mu,sigma2);
    //cube q1=up4[0];mat q2 = up4[1],q3 = up4[2];
    //z = q1;v1 = q2;v2= q3;
    //m = randg( distr_param(am+r,bm+har(d[0])));
    //rho = R::rbeta(arho+count(c2.begin(),c2.end(),1),brho+count(c2.begin(),c2.end(),0));
    //List up5 = update_b(d,mub,sigmab,b1,b2,c1,c2,c3,lambda1,lambda2,z);
    //mat s1=up5[0],s2=up5[1];
    //b1 = s1;b2=s2;
    //mat up6 = update_norm(z,y,as,bs,sigmamu,mu,5,v1,v2);
    //mu = up6.rows(0,d[1]-1);sigma2 = up6.rows(d[1],2*d[1]-1);
    //matrix_gibbs[it] = join_cols(c1,c2,c3);
  //}
  //return matrix_gibbs;
//}