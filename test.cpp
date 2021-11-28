#include<RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

using namespace std;
//[[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
uvec test(vec zero_col)
{
  uvec a= find_unique(zero_col);
  return a;
}
  
      