// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

/*
 C++ helper for the longitudinal graphical‑lasso routine.
 Returns a numeric vector of length two

 [1] = -0.5 * p * log(det(A))
 -0.5 * tr(A⁻¹ %*% X %*% B %*% t(X))

 [2] = gradient of the above w.r.t. the correlation decay
 parameter τ (scalar) for the “expFixed” model.

 Arguments
 X   : n × p matrix of measurements for a single subject
 (rows = observations, columns = variables)
 t   : numeric vector of length n – the time points belonging to the rows B   : p × p precision matrix (symmetric positive‑definite)
 tau : positive scalar – decay parameter in the correlation model

 The correlation matrix A has entries
 A_{ab} = exp( -tau * |t_a - t_b| ).

 Its derivative w.r.t. τ is
 dA_{ab}/dτ = -|t_a - t_b| * A_{ab}.
 */

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector lglasso_cpp_obj(const arma::mat& X,          // n × p
                              const NumericVector& t,    // length n
                              const arma::mat& B,        // p × p
                              const double tau) {

  // -------------------------------------------------
  // 1) Build the correlation matrix A (n × n)
  // -------------------------------------------------
  int n = t.size();                              // number of observations
  arma::mat A(n, n, arma::fill::zeros);          // <-- fixed line
  for (int i = 0; i < n; ++i) {
    for (int j = i; j < n; ++j) {
      double d   = std::abs(t[i] - t[j]);
      double val = std::exp(-tau * d);
      A(i, j)    = val;
      A(j, i)    = val;                         // symmetry
    }
  }

  // -------------------------------------------------
  // 2) Cholesky of A  (A = L %*% Lᵀ)
  // -------------------------------------------------
  arma::mat L = arma::chol(A, "lower");
  double logdetA = 2.0 * arma::sum(arma::log(L.diag()));

  // -------------------------------------------------
  // 3) Compute M = X %*% B %*% t(X)   (n × n)
  // -------------------------------------------------
  arma::mat M = X * B * X.t();                 // n × n

  // -------------------------------------------------
  // 4) Inverse of A (still cheap because n is tiny)
  // -------------------------------------------------
  arma::mat Linv = arma::solve(arma::trimatl(L),
                               arma::eye<arma::mat>(n, n));
  arma::mat Ainv = Linv.t() * Linv;            // n × n

  // -------------------------------------------------
  // 5) Trace term tr( A⁻¹ %*% M )
  // -------------------------------------------------
  double trace_term = arma::trace(Ainv * M);

  // -------------------------------------------------
  // 6) Gradient w.r.t. τ
  // -------------------------------------------------
  // dA/dτ matrix (same shape as A)
  arma::mat dA(n, n, arma::fill::zeros);
  for (int i = 0; i < n; ++i) {
    for (int j = i; j < n; ++j) {
      double dist = std::abs(t[i] - t[j]);
      double dAij = -dist * A(i, j);            // ∂A_{ij}/∂τ
      dA(i, j) = dAij;
      dA(j, i) = dAij;                         // symmetry
    }
  }

  // d log|A| / dτ = tr( A⁻¹ dA )
  double dlogdet_dtau = arma::trace(Ainv * dA);

  // d tr( A⁻¹ M ) / dτ = - tr( A⁻¹ dA A⁻¹ M )
  //                  = - tr( (A⁻¹ M A⁻¹) dA )
  arma::mat G = Ainv * M * Ainv;               // n × n
  double dtrace_dtau = - arma::trace(G * dA);

  // -------------------------------------------------
  // 7) Return the (negative) log‑likelihood and its gradient
  // -------------------------------------------------
  int p = X.n_cols;                           // number of variables

  NumericVector out(2);
  out[0] = -0.5 * p * logdetA - 0.5 * trace_term;          // objective value
  out[1] = -0.5 * p * dlogdet_dtau - 0.5 * dtrace_dtau;    // gradient
  return out;
}
