/*
 * https://github.com/coatless/Rmath/blob/master/plogis.c
 * 
 * 
 * 
 * 
 * 
 */


#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
NumericMatrix log1pexp(NumericMatrix mat) {
  
  
  
  int m = mat.nrow();
  int n = mat.ncol();
  
  NumericMatrix res(m, n);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++) {
      res(i, j) = log1pexp(mat(i, j)); // log(1 + exp(x))
    }
    
  return res;
} 