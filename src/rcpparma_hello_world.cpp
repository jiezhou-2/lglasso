// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
arma::mat rcpparma_hello_world() {
    arma::mat m1 = arma::eye<arma::mat>(3, 3);
    arma::mat m2 = arma::eye<arma::mat>(3, 3);
	                     
    return m1 + 3 * (m1 + m2);
}


// another simple example: outer product of a vector, 
// returning a matrix
//
// [[Rcpp::export]]
arma::mat rcpparma_outerproduct(const arma::colvec & x) {
    arma::mat m = x * x.t();
    return m;
}

// and the inner product returns a scalar
//
// [[Rcpp::export]]
double rcpparma_innerproduct(const arma::colvec & x) {
    double v = arma::as_scalar(x.t() * x);
    return v;
}


// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
Rcpp::List rcpparma_bothproducts(const arma::colvec & x) {
    arma::mat op = x * x.t();
    double    ip = arma::as_scalar(x.t() * x);
    return Rcpp::List::create(Rcpp::Named("outer")=op,
                              Rcpp::Named("inner")=ip);
}



// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]


arma::mat phifunction(const arma::vec t, double tau, double ty){
  int m=t.n_elem;
  arma::mat M=arma::eye(m,m);
  if (m==1){
    M=1;
  }
  else{
    for (int i=1; i<=m; i++) {
      for (int j=i+1; j<=m; j++){
        M(i-1,j-1)=exp(-tau*pow((abs(t(i-1)-t(j-1))),ty));
        M(j-1,i-1)=M(i-1,j-1);
      }
    }
  
  }
  return M;
}



// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]


arma::mat iss(arma::mat idata, double itau, double ty){
  arma::vec t=idata.col(1);
  int p=idata.n_cols-2;
  int n=idata.n_rows;
  double tau;
  arma::mat inversephi;
  inversephi=arma::inv(phifunction(t=t,tau=itau,ty=ty));
  arma::mat si(p,p);
  si.fill(0);
  arma::mat yy;
  yy=idata.cols(2,p+1);
  arma::rowvec mm;
  mm=arma::mean(yy);
  arma::mat xx;
  xx.ones(n,1);
  arma::mat mu;
  mu=xx*mm;
  arma::mat cc;
  cc=yy-mu;
  arma::rowvec bb1(n);
  arma::rowvec bb2(n);
  for (int i=0; i< n-1; i++){
  bb1= cc.row(i);
    for (int j=0; j< n-1; j++){
  bb2= cc.row(j);
  si+=inversephi(i,j)*(bb1.t()*bb2);
    }
  }
  return  si;
}


