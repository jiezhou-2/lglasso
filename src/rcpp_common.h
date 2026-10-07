// src/rcpp_common.h ---------------------------------------------------------
#ifndef RCPP_COMMON_H
#define RCPP_COMMON_H

/*  *** DO NOT PUT ANY RCPP ATTRIBUTES HERE ***  */

/*  Include the heavy‑weight header once */
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <unordered_map>

/*  Bring the most convenient symbols into the global namespace */
using namespace Rcpp;      // Rcpp::NumericVector, Rcpp::List, …
using arma::mat;           // arma::mat  (matrix of doubles)
using arma::vec;           // arma::vec  (column vector)

/*  Small helper macro – makes the export line shorter */
#define EXPORT_RCPP [[Rcpp::export]]

#endif // RCPP_COMMON_H
