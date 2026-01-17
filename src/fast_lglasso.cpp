// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Helper that returns a numeric vector of length 2:
//   [1] = -0.5 * p * log(det(A))
//   [2] = gradient w.r.t. tau (scalar) for the 'expFixed' correlation model
// Arguments:
//   X   : p × n_i matrix of measurements for a single subject (already transposed)
//   t   : vector of time points (length = n_i)
//   B   : current precision matrix (p × p, SPD)
//   tau : current decay parameter (scalar)
//
// The correlation matrix A has entries exp(-tau * |t_a - t_b|)
// Its derivative w.r.t. tau is -|t_a - t_b| * A.
//
// Because p is usually modest (≤ 200) we can afford a full
// Cholesky decomposition of A each time.
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector lglasso_cpp_obj(const arma::mat& X,
                              const NumericVector& t,
                              const arma::mat& B,
                              const double tau) {

  // 1️⃣ compute A and A⁻¹ efficiently via Cholesky
  //    (A is symmetric positive definite for any tau>0)
  int n = t.size();
  arma::mat A(n, n, arma::fill::zeros);
  for (int i = 0; i < n; ++i) {
    for (int j = i; j < n; ++j) {
      double d = std::abs(t[i] - t[j]);
      double val = std::exp(-tau * d);
      A(i, j) = val;
      A(j, i) = val;  // symmetric
    }
  }

  // Cholesky factor L: A = L * L.t()
  arma::mat L = arma::chol(A, "lower");

  // log(det(A)) = 2 * sum(log(diag(L)))
  double logdetA = 2.0 * arma::sum(arma::log(L.diag()));

  // Solve A⁻¹ * X efficiently: first solve L * Y = X, then L.t() * Z = Y
  arma::mat Y = arma::solve(arma::trimatl(L), X);
  arma::mat AinvX = arma::solve(arma::trimatu(L.t()), Y);

  // trace term = tr( B %*% X.t() %*% A⁻¹ %*% X )
  //   = sum( B %*% (X.t() %*% A⁻¹ %*% X) ) = sum( B %*% (X.t() %*% AinvX) )
  arma::mat X_t_Ainv_X = X.t() * AinvX;       // p × p
  double trace_term = arma::trace(B * X_t_Ainv_X);

  // ----- gradient w.r.t. tau ----------------------------------------
  // dA/dtau = -|t_a - t_b| * A
  // Hence, d log det(A) / d tau = trace( A⁻¹ %*% dA/dtau )
  //   = - sum_{a,b} |t_a - t_b| * (A⁻¹ ∘ A)_{ab}
  // where ∘ is the Hadamard product.
  // Instead of forming A⁻¹ explicitly we exploit the Cholesky factors:
  //   A⁻¹ = (L⁻¹)ᵗ %*% L⁻¹
  arma::mat Linv = arma::solve(arma::trimatl(L), arma::eye<arma::mat>(n, n));
  arma::mat Ainv = Linv.t() * Linv;

  double dlogdet_dtau = 0.0;
  double dtrace_dtau  = 0.0;
  for (int a = 0; a < n; ++a) {
    for (int b = a; b < n; ++b) {
      double dist = std::abs(t[a] - t[b]);
      double weight = -dist * A(a, b); // dA/dtau element
      double a_inv = Ainv(a, b);
      dlogdet_dtau += weight * a_inv;               // symmetric, count once
      // For the trace part we need   tr( B %*% X.t() %*% dA⁻¹ %*% X )
      // dA⁻¹ = -A⁻¹ (dA/dtau) A⁻¹  (matrix calculus)
      // The scalar contribution becomes  - weight * (Xᵗ A⁻¹)_{·a} B (A⁻¹ X)_{b·}
      // We compute it via the already‑available AinvX:
      arma::colvec col_a = AinvX.col(a);
      arma::colvec col_b = AinvX.col(b);
      // B * col_a
      arma::colvec Bcol_a = B * col_a;
      double term = arma::dot(Bcol_a, col_b); // = col_aᵗ B col_b
      dtrace_dtau += weight * term; // again symmetric; we count each pair once
      if (a != b) { // off‑diag appears twice in the full sum
        dlogdet_dtau += weight * a_inv;
        dtrace_dtau  += weight * term;
      }
    }
  }

  // The contribution to the *negative* log‑likelihood is
  //   -0.5 * p * logdetA - 0.5 * trace_term
  // So we return the (positive) objective parts that the R code
  // will multiply by -1.
  NumericVector out(2);
  out[0] = -0.5 * X.n_rows * logdetA - 0.5 * trace_term;   // obj value
  out[1] = -0.5 * X.n_rows * dlogdet_dtau - 0.5 * dtrace_dtau; // d obj / d tau
  return out;
}
