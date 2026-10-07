#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <numeric>          // std::iota
#include <R_ext/Random.h>   // R_unif_index
#include <R_ext/Applic.h>   // lbfgsb(), the optimizer behind optim(method = "L-BFGS-B")
// [[Rcpp::depends(RcppArmadillo)]]   // harmless in a package; needed only for sourceCpp()
//' Exponential correlation matrix
//'
//' Computes \eqn{M_{ij} = \exp(-\tau |t_i - t_j|^{expFix})} with ones on the diagonal.
//'
//' @param t Numeric vector of time points.
//' @param tau Positive scalar decay parameter.
//' @param expFix Exponent applied to the distance (default 1).
//' @return An n x n numeric matrix.
//' @examples
//' phifunction(1:5, tau = 0.3)
//' phifunction(c(0, 2, 5), tau = 0.7, expFix = 2)
//' @export
 // [[Rcpp::export]]
 arma::mat phifunction(const arma::vec& t, double tau, double expFix = 1.0) {
   if (tau <= 0) Rcpp::stop("tau should be positive!");
   const arma::uword n = t.n_elem;
   arma::mat M(n, n, arma::fill::ones);
   for (arma::uword j = 0; j < n; ++j) {
     for (arma::uword i = j + 1; i < n; ++i) {
       double d = std::pow(std::abs(t(i) - t(j)), expFix);
       M(i, j) = M(j, i) = std::exp(-tau * d);
     }
   }
   return M;
 }



// [[Rcpp::depends(RcppArmadillo)]]

// ---- helpers -------------------------------------------------------------

// Returns a numeric matrix: column 0 = time (2nd column of datai),
// columns 1..p = outcomes (columns 3.. of datai).
// The first column (subject ID) is skipped, so it may be character/factor.
static arma::mat dataToMatrix(SEXP x) {
  if (Rf_inherits(x, "data.frame")) {
    Rcpp::DataFrame df(x);
    const int nc = df.size();
    if (nc < 3) Rcpp::stop("datai must have at least 3 columns (id, time, variables)");
    const int nr = df.nrows();
    arma::mat M(nr, nc - 1);
    for (int j = 1; j < nc; ++j) {
      Rcpp::NumericVector col = Rcpp::as<Rcpp::NumericVector>(df[j]);
      M.col(j - 1) = arma::vec(col.begin(), nr, false, true);   // copy
    }
    return M;
  }
  if (Rf_isMatrix(x)) {
    arma::mat A = Rcpp::as<arma::mat>(x);
    if (A.n_cols < 3) Rcpp::stop("datai must have at least 3 columns (id, time, variables)");
    return A.cols(1, A.n_cols - 1);
  }
  Rcpp::stop("datai must be a data.frame or a numeric matrix");
}

// Uses R's own as.character() so the labels match names(wi) exactly as in R
// (e.g. factors -> labels, 1e5 -> "1e+05").
static std::vector<std::string> groupToString(SEXP g) {
  static Rcpp::Function asChar("as.character");
  Rcpp::CharacterVector cv = asChar(g);
  std::vector<std::string> out(cv.size());
  for (R_xlen_t i = 0; i < cv.size(); ++i) {
    if (Rcpp::CharacterVector::is_na(cv[i])) Rcpp::stop("groupi must not contain NA");
    out[i] = Rcpp::as<std::string>(cv[i]);
  }
  return out;
}

// Split row indices by group label, keeping order of first appearance (= unique()).
// lev[k] is the k-th group label, idx[k] the row indices belonging to it.
static void splitByGroup(const std::vector<std::string>& g,
                         std::vector<std::string>& lev,
                         std::vector<std::vector<arma::uword>>& idx) {
  lev.clear();
  idx.clear();
  std::unordered_map<std::string, std::size_t> pos;
  for (arma::uword r = 0; r < g.size(); ++r) {
    auto it = pos.find(g[r]);
    std::size_t k;
    if (it == pos.end()) {
      k = lev.size();
      pos.emplace(g[r], k);
      lev.push_back(g[r]);
      idx.emplace_back();
    } else {
      k = it->second;
    }
    idx[k].push_back(r);
  }
}


// Extract the time column (2nd column) of datai as a numeric vector.
// Works for data.frame / tibble and numeric matrix; the ID column is never touched.
static arma::vec timeColumn(SEXP x) {
  if (Rf_inherits(x, "data.frame")) {
    Rcpp::DataFrame df(x);
    if (df.size() < 2) Rcpp::stop("datai must have at least 2 columns (id, time, ...)");
    Rcpp::NumericVector col = Rcpp::as<Rcpp::NumericVector>(df[1]);
    return arma::vec(col.begin(), col.size());          // copies the data
  }
  if (Rf_isMatrix(x)) {
    arma::mat A = Rcpp::as<arma::mat>(x);
    if (A.n_cols < 2) Rcpp::stop("datai must have at least 2 columns (id, time, ...)");
    return A.col(1);
  }
  Rcpp::stop("datai must be a data.frame or a numeric matrix");
}

// Extract sample (1st column) and weight (2nd column) from the output of importanceSample().
static void sampleAndWeight(SEXP x, arma::vec& s, arma::vec& w) {
  if (Rf_inherits(x, "data.frame")) {
    Rcpp::DataFrame df(x);
    if (df.size() < 2) Rcpp::stop("importancesSample must have 2 columns (sample, weight)");
    Rcpp::NumericVector a = Rcpp::as<Rcpp::NumericVector>(df[0]);
    Rcpp::NumericVector b = Rcpp::as<Rcpp::NumericVector>(df[1]);
    s = arma::vec(a.begin(), a.size());
    w = arma::vec(b.begin(), b.size());
    return;
  }
  if (Rf_isMatrix(x)) {
    arma::mat A = Rcpp::as<arma::mat>(x);
    if (A.n_cols < 2) Rcpp::stop("importancesSample must have 2 columns (sample, weight)");
    s = A.col(0);
    w = A.col(1);
    return;
  }
  Rcpp::stop("importancesSample must be a data.frame or a numeric matrix");
}

// First column of a data.frame / matrix (subject ID); returned as RObject so it stays protected.
static Rcpp::RObject firstColumn(SEXP x) {
  if (Rf_inherits(x, "data.frame")) {
    Rcpp::DataFrame df(x);
    if (df.size() < 1) Rcpp::stop("data has no columns");
    return Rcpp::RObject(df[0]);
  }
  if (Rf_isMatrix(x)) {
    Rcpp::NumericMatrix A(x);
    return Rcpp::RObject(Rcpp::NumericVector(A(Rcpp::_, 0)));
  }
  Rcpp::stop("data must be a data.frame or a numeric matrix");
}

// P(j,k) = |t_j - t_k|^expFix, so that phifunction(t, tau, expFix) == exp(-tau * P)
static arma::mat powDist(const arma::vec& t, double expFix) {
  const arma::uword n = t.n_elem;
  arma::mat P(n, n, arma::fill::zeros);
  const bool linear = (expFix == 1.0);
  for (arma::uword j = 0; j < n; ++j)
    for (arma::uword i = j + 1; i < n; ++i) {
      double d = std::abs(t[i] - t[j]);
      if (!linear) d = std::pow(d, expFix);
      P(i, j) = P(j, i) = d;
    }
    return P;
}

static void checkExpFix(double expFix) {
  if (!std::isfinite(expFix) || expFix <= 0)
    Rcpp::stop("expFix must be a positive finite number");
  if (expFix > 2)
    Rcpp::warning("expFix > 2: correlation matrices may not be positive definite");
}

// tau-independent pieces for one subject, computed once
struct TauPrep {
  std::vector<std::string> lev;    // group labels, order of first appearance
  std::vector<arma::vec>   times;  // time points per group
  std::vector<arma::mat>   Y;      // n_g x p outcomes per group
  std::vector<arma::mat>   P;      // |t_j - t_k|^expFix per group
  std::vector<arma::mat>   W;      // p x p precision matrix per group
  double p = 0.0;
};

// X: col 0 = time, cols 1..p = variables (output of dataToMatrix); g: group label per row
static TauPrep prepareTau(const arma::mat& X, const std::vector<std::string>& g,
                          Rcpp::List wi, double expFix) {
  if (g.size() != X.n_rows)
    Rcpp::stop("group should be the same length as the number of rows of data!");
  if (wi.size() == 0) Rcpp::stop("wi must be a non-empty list");
  if (Rf_isNull(wi.names())) Rcpp::stop("wi must be a named list (names = group labels)");

  const arma::uword p = X.n_cols - 1;
  if (Rcpp::as<arma::mat>(wi[0]).n_rows != p)
    Rcpp::stop("number of variables in data (%d) differs from nrow(wi[[1]]) (%d)",
               (int)p, (int)Rcpp::as<arma::mat>(wi[0]).n_rows);

  Rcpp::CharacterVector wnames = wi.names();
  std::unordered_map<std::string, int> wpos;
  for (int k = 0; k < wnames.size(); ++k)
    wpos.emplace(Rcpp::as<std::string>(wnames[k]), k);   // first match, like [[ ]]

  TauPrep pr;
  pr.p = static_cast<double>(p);
  std::vector<std::vector<arma::uword>> idx;
  splitByGroup(g, pr.lev, idx);

  for (std::size_t k = 0; k < pr.lev.size(); ++k) {
    auto wit = wpos.find(pr.lev[k]);
    if (wit == wpos.end())
      Rcpp::stop("no precision matrix in wi named '%s'", pr.lev[k]);
    arma::mat W = Rcpp::as<arma::mat>(wi[wit->second]);
    if (W.n_rows != p || W.n_cols != p)
      Rcpp::stop("wi[['%s']] must be a %d x %d matrix", pr.lev[k], (int)p, (int)p);

    const arma::uvec rows = arma::conv_to<arma::uvec>::from(idx[k]);
    const arma::mat  Xg   = X.rows(rows);          // materialise first, then .col/.cols
    const arma::vec  t    = Xg.col(0);
    pr.times.push_back(t);
    pr.Y.push_back(Xg.cols(1, p));
    pr.P.push_back(powDist(t, expFix));
    pr.W.push_back(W);
  }
  return pr;
}

// log alpha - alpha*tau + sum_g [ -p/2 log|Phi_g| - 1/2 tr(Y_g' Phi_g^{-1} Y_g W_g) ]
static double logDensityTau(double tau, const TauPrep& pr, double alpha) {
  double loglik = 0.0;
  for (std::size_t k = 0; k < pr.Y.size(); ++k) {
    const arma::mat Phi = arma::exp(-tau * pr.P[k]);   // == phifunction(t, tau, expFix)
    double logdet;
    arma::mat A, L;
    if (arma::chol(L, Phi, "lower")) {
      logdet = 2.0 * arma::accu(arma::log(L.diag()));
      const arma::mat Z = arma::solve(arma::trimatl(L), pr.Y[k]);
      A = Z.t() * Z;
    } else {
      double val, sign;
      arma::log_det(val, sign, Phi);
      logdet = (sign > 0) ? val : std::numeric_limits<double>::quiet_NaN();
      A = pr.Y[k].t() * arma::solve(Phi, pr.Y[k]);    // throws -> R error if singular
    }
    loglik += -0.5 * pr.p * logdet - 0.5 * arma::accu(A % pr.W[k]);
  }
  return std::log(alpha) + loglik - alpha * tau;
}

// Draw n taus from Exp(alpha), keep finite ones, return samples + normalised weights
static void importanceSampleCore(int n, const TauPrep& pr, double alpha,
                                 std::vector<double>& samp, std::vector<double>& weight) {
  const double scale = 1.0 / alpha;              // R::rexp / R::dexp use scale = 1/rate
  std::vector<double> draws(n);
  for (int i = 0; i < n; ++i) draws[i] = R::rexp(scale);

  samp.clear(); weight.clear();
  std::vector<double> logw;
  samp.reserve(n); logw.reserve(n);
  for (int i = 0; i < n; ++i) {
    if (i % 100 == 0) Rcpp::checkUserInterrupt();
    const double l1 = logDensityTau(draws[i], pr, alpha);
    if (!std::isfinite(l1)) continue;
    samp.push_back(draws[i]);
    logw.push_back(l1 - R::dexp(draws[i], scale, 1));
  }
  if (samp.empty()) Rcpp::stop("No valid samples are generated!");

  const double mx = *std::max_element(logw.begin(), logw.end());
  double s = 0.0;
  for (double w : logw) s += std::exp(w - mx);
  const double lse = mx + std::log(s);
  weight.resize(logw.size());
  for (std::size_t k = 0; k < logw.size(); ++k) weight[k] = std::exp(logw[k] - lse);
}




// ---- objective for AA ----------------------------------------------------
struct AAProblem {
  std::vector<arma::mat> Y;   // per subject: n_i x p outcome matrix
  std::vector<arma::mat> P;   // per subject: n_i x n_i powered distances
  arma::mat W;                // p x p precision matrix
  double p = 0.0;
  double lower = 0.0, upper = 0.0;
  bool failed = false;
  std::string msg;
};

// Negative log-likelihood: sum_i [ p/2 log|Phi_i| + 1/2 tr(Phi_i^{-1} Y_i W Y_i') ]
static double aaNegLogLik(double tau, const AAProblem& pr) {
  double obj = 0.0;
  for (std::size_t i = 0; i < pr.Y.size(); ++i) {
    const arma::mat Phi = arma::exp(-tau * pr.P[i]);
    double logdet, tr;
    arma::mat L;
    if (arma::chol(L, Phi, "lower")) {
      logdet = 2.0 * arma::accu(arma::log(L.diag()));
      const arma::mat Z = arma::solve(arma::trimatl(L), pr.Y[i]);
      tr = arma::accu((Z * pr.W) % Z);                       // tr(Z W Z')
    } else {                                                  // not numerically PD
      double val, sign;
      if (!arma::log_det(val, sign, Phi) || sign <= 0)
        return std::numeric_limits<double>::quiet_NaN();     // log(det) < 0 -> NaN, as in R
      arma::mat X;
      if (!arma::solve(X, Phi, pr.Y[i], arma::solve_opts::no_approx))
        return std::numeric_limits<double>::quiet_NaN();
      logdet = val;
      tr = arma::accu((pr.Y[i].t() * X) % pr.W.t());         // tr(Y' Phi^-1 Y W)
    }
    obj += 0.5 * pr.p * logdet + 0.5 * tr;
  }
  return obj;
}

static const double AA_FAILVAL = 1e300;

static double aaSafeEval(double tau, AAProblem* pr) {
  double v;
  try {
    v = aaNegLogLik(tau, *pr);
  } catch (std::exception& e) {
    pr->failed = true; pr->msg = e.what(); return AA_FAILVAL;
  }
  if (!std::isfinite(v)) {
    pr->failed = true; pr->msg = "L-BFGS-B needs finite values of 'fn'"; return AA_FAILVAL;
  }
  return v;
}

// callbacks for lbfgsb()
static double aa_fn(int, double* x, void* ex) {
  AAProblem* pr = static_cast<AAProblem*>(ex);
  if (pr->failed) return AA_FAILVAL;
  return aaSafeEval(x[0], pr);
}

// central difference exactly as optim() does for L-BFGS-B (ndeps = 1e-3, truncated at bounds)
static void aa_gr(int, double* x, double* df, void* ex) {
  AAProblem* pr = static_cast<AAProblem*>(ex);
  if (pr->failed) { df[0] = 0.0; return; }
  const double x0 = x[0];
  double eps = 1e-3, epsused = 1e-3;
  double tmp = x0 + eps;
  if (tmp > pr->upper) { tmp = pr->upper; epsused = tmp - x0; }
  const double v1 = aaSafeEval(tmp, pr);
  tmp = x0 - eps;
  if (tmp < pr->lower) { tmp = pr->lower; eps = x0 - tmp; }
  const double v2 = aaSafeEval(tmp, pr);
  df[0] = (v1 - v2) / (epsused + eps);
  if (pr->failed || !std::isfinite(df[0])) {
    if (!pr->failed) { pr->failed = true; pr->msg = "non-finite finite-difference value"; }
    df[0] = 0.0;
  }
}




//' density function in EM algorithm
//'
//' Computes conditional likelihood for heterogeneous models
//'
//' @param tau the dampening rate
//' @param datai the data set for subject i (col 1 = id, col 2 = time, rest = variables)
//' @param wi named list of precision matrices, names matching the values of groupi
//' @param alpha the rate in exponential distribution
//' @param groupi the data point indices (group labels), one per row of datai
//' @param expFix a scalar specifying the form of the correlation function
//' @noRd
//' @returns a numeric standing for the likelihood for a given subject
// [[Rcpp::export]]
double conDensityTau(double tau, SEXP datai, Rcpp::List wi, double alpha,
                     SEXP groupi, double expFix = 1.0) {
  if (!(tau > 0) || !std::isfinite(tau)) Rcpp::stop("tau must be a positive finite number");
  checkExpFix(expFix);
  const arma::mat X = dataToMatrix(datai);
  const TauPrep pr = prepareTau(X, groupToString(groupi), wi, expFix);
  return logDensityTau(tau, pr, alpha);
}


// [[Rcpp::depends(RcppArmadillo)]]


//' Function for generating the samples from posterior distribution in EM algorithm
//'
//' @param n the number of random samples
//' @param datai the data for subject i
//' @param wi the given precision matrix (named list, names = group labels)
//' @param alpha the exponential distribution with rate alpha
//' @param groupi specify how datai is grouped
//' @param expFix a scalar specifying the form of the correlation function
//' @noRd
//' @returns a data frame for samples and their weights
 // [[Rcpp::export]]
 Rcpp::DataFrame importanceSample(int n, SEXP datai, Rcpp::List wi, double alpha,
                                  SEXP groupi, double expFix = 1.0) {
   if (n <= 0) Rcpp::stop("n must be a positive integer");
   if (!(alpha > 0) || !std::isfinite(alpha)) Rcpp::stop("alpha must be a positive finite number");
   checkExpFix(expFix);

   const arma::mat X = dataToMatrix(datai);
   const TauPrep pr = prepareTau(X, groupToString(groupi), wi, expFix);   // once, not n times

   std::vector<double> samp, weight;
   importanceSampleCore(n, pr, alpha, samp, weight);
   return Rcpp::DataFrame::create(
     Rcpp::Named("sample") = Rcpp::NumericVector(samp.begin(), samp.end()),
     Rcpp::Named("weight") = Rcpp::NumericVector(weight.begin(), weight.end()));
 }


//' Function for computing estimates in EM algorithm
//'
//' @param importancesSample the random samples (output of importanceSample: columns sample, weight)
//' @param datai the data for subject i (col 1 = id, col 2 = time, ...)
//' @param groupi specify how datai is grouped
//' @param expFix a scalar specifying the form of the correlation function
//' @noRd
//' @returns a list with estimateTau and a named list estimatePhi
 // [[Rcpp::export]]
 Rcpp::List importanceEstimates(SEXP importancesSample, SEXP datai, SEXP groupi,
                                double expFix = 1.0) {
   // 1. Weighted estimate of tau
   arma::vec s, w;
   sampleAndWeight(importancesSample, s, w);
   if (s.n_elem == 0) Rcpp::stop("importancesSample contains no samples");
   const double estimateTau = arma::dot(s, w);     // = sum(sample * weight)

   // 2. Split time points by group (order of first appearance)
   const arma::vec tt = timeColumn(datai);
   const std::vector<std::string> g = groupToString(groupi);
   if (g.size() != tt.n_elem)
     Rcpp::stop("groupi should have the same length as the number of rows of datai!");

   std::vector<std::string> lev;
   std::vector<std::vector<arma::uword>> idx;
   splitByGroup(g, lev, idx);

   // 3. Correlation matrix per group
   Rcpp::List estimatePhi(lev.size());
   Rcpp::CharacterVector nms(lev.size());
   for (std::size_t k = 0; k < lev.size(); ++k) {
     const arma::uvec rows = arma::conv_to<arma::uvec>::from(idx[k]);
     const arma::vec  tk   = tt.elem(rows);
     estimatePhi[k] = Rcpp::wrap(phifunction(tk, estimateTau, expFix));
     nms[k] = lev[k];
   }
   estimatePhi.names() = nms;

   return Rcpp::List::create(Rcpp::Named("estimateTau") = estimateTau,
                             Rcpp::Named("estimatePhi") = estimatePhi);
 }




//' Find the tau's for homogeneous model
//' @param B a list of length 1 or 2 (or a single matrix). Each entry is a p by p
//'   precision matrix; with length 2 they represent the pre- and post-treatment network.
//' @param data a data frame (columns: id, time, p variables) or a list of such data frames
//' @param expFix the parameter in the correlation function
//' @param maxit currently unused (optim defaults are used, as in the R version)
//' @param tol currently unused (optim defaults are used, as in the R version)
//' @param lower lower bound(s) for tau; element k is used for data set k (recycled)
//' @param upper upper bound(s) for tau; element k is used for data set k (recycled)
//' @noRd
//' @returns a list with corMatrix (list of named lists of matrices) and tau
 // [[Rcpp::export]]
 Rcpp::List AA(SEXP B, SEXP data, double expFix = 1.0, int maxit = 50, double tol = 1e-4,
               Rcpp::NumericVector lower = Rcpp::NumericVector::create(0.01, 0.1),
               Rcpp::NumericVector upper = Rcpp::NumericVector::create(10.0, 5.0)) {
   (void)maxit; (void)tol;   // kept for compatibility with the R signature

   if (!std::isfinite(expFix) || expFix <= 0) Rcpp::stop("expFix must be a positive finite number");
   if (lower.size() == 0 || upper.size() == 0) Rcpp::stop("lower and upper must not be empty");

   // B: single matrix -> list
   Rcpp::List Bl;
   if (Rf_isMatrix(B))              Bl = Rcpp::List::create(B);
   else if (TYPEOF(B) == VECSXP)    Bl = Rcpp::List(B);
   else Rcpp::stop("B must be a matrix or a list of matrices");

   // data: single data.frame / matrix -> list
   Rcpp::List Dl;
   if (Rf_inherits(data, "data.frame") || Rf_isMatrix(data)) Dl = Rcpp::List::create(data);
   else if (TYPEOF(data) == VECSXP)                          Dl = Rcpp::List(data);
   else Rcpp::stop("data must be a data.frame or a list of data.frames");

   if (Bl.size() != Dl.size()) Rcpp::stop("B should have same length as data!");

   const R_xlen_t K = Bl.size();
   Rcpp::List corMatrix(K);
   Rcpp::NumericVector Tau(K);

   for (R_xlen_t k = 0; k < K; ++k) {
     SEXP dd = Dl[k];
     const arma::mat X = dataToMatrix(dd);            // col 0 = time, cols 1..p = variables
     const arma::uword p = X.n_cols - 1;
     const arma::mat W = Rcpp::as<arma::mat>(Bl[k]);
     if (W.n_rows != p || W.n_cols != p)
       Rcpp::stop("B[[%d]] must be a %d x %d matrix (p = ncol(data) - 2)", (int)(k + 1), (int)p, (int)p);
     if (!X.is_finite())
       Rcpp::stop("data[[%d]]: time and variable columns must not contain NA/NaN/Inf", (int)(k + 1));

     // split rows by subject ID (order of first appearance)
     Rcpp::RObject idcol = firstColumn(dd);
     const std::vector<std::string> g = groupToString(idcol);
     if (g.size() != X.n_rows) Rcpp::stop("internal error: id column length mismatch");
     std::vector<std::string> lev;
     std::vector<std::vector<arma::uword>> idx;
     splitByGroup(g, lev, idx);

     AAProblem pr;
     pr.W = W;
     pr.p = static_cast<double>(p);
     pr.lower = lower[k % lower.size()];
     pr.upper = upper[k % upper.size()];
     if (!(pr.lower > 0) || !std::isfinite(pr.upper) || pr.lower > pr.upper)
       Rcpp::stop("need 0 < lower <= upper < Inf (data set %d)", (int)(k + 1));

     std::vector<arma::vec> times(lev.size());
     pr.Y.reserve(lev.size());
     pr.P.reserve(lev.size());
     for (std::size_t s = 0; s < lev.size(); ++s) {
       const arma::uvec rows = arma::conv_to<arma::uvec>::from(idx[s]);
       const arma::mat  Xs   = X.rows(rows);        // copy the subject's rows into a real matrix
       times[s] = Xs.col(0);                        // time column
       pr.Y.push_back(Xs.cols(1, p));               // outcome columns
       pr.P.push_back(powDist(times[s], expFix));
     }

     // optim(1, likefun, method = "L-BFGS-B", lower, upper) with optim's defaults
     double x = 1.0, lo = pr.lower, up = pr.upper, Fmin = 0.0;
     int nbd = 2, fail = 0, fncount = 0, grcount = 0;
     char msg[60];
     lbfgsb(1, 5, &x, &lo, &up, &nbd, &Fmin, aa_fn, aa_gr, &fail, &pr,
            1e7, 0.0, &fncount, &grcount, 100, msg, 0, 10);

     if (pr.failed)
       Rcpp::stop("optimisation of tau failed for data set %d: %s", (int)(k + 1), pr.msg);
     if (fail != 0)
       Rcpp::warning("L-BFGS-B for data set %d: %s (code %d)", (int)(k + 1), std::string(msg), fail);

     Tau[k] = x;

     // correlation matrices at the estimated tau, named by subject ID
     Rcpp::List A(lev.size());
     Rcpp::CharacterVector nms(lev.size());
     for (std::size_t s = 0; s < lev.size(); ++s) {
       A[s] = Rcpp::wrap(phifunction(times[s], x, expFix));
       nms[s] = lev[s];
     }
     A.names() = nms;
     corMatrix[k] = A;
   }

   return Rcpp::List::create(Rcpp::Named("corMatrix") = corMatrix,
                             Rcpp::Named("tau")       = Tau);
 }



//' Estimate the phimatrix in heterogeneous model
//'
//' @param data longitudinal data set (col 1 = subject id, col 2 = time, rest = variables)
//' @param wi named list of precision matrices (names = group labels)
//' @param alpha exponential distribution with rate alpha
//' @param group specify how data is grouped (one label per row of data)
//' @param l number of random samples in importance sampling
//' @param expFix a scalar specifying the form of the correlation function
//' @noRd
//' @returns a list with Tau (subjects x 1 matrix) and AA (per group, per subject correlation matrices)
// [[Rcpp::export]]
 Rcpp::List AAheter(SEXP data, Rcpp::List wi, double alpha, SEXP group,
                    int l = 5000, double expFix = 1.0) {
   if (l <= 0) Rcpp::stop("l must be a positive integer");
   if (!(alpha > 0) || !std::isfinite(alpha)) Rcpp::stop("alpha must be a positive finite number");
   checkExpFix(expFix);

   const arma::mat X = dataToMatrix(data);                  // col 0 = time, cols 1..p = variables
   Rcpp::RObject idcol = firstColumn(data);
   const std::vector<std::string> ids = groupToString(idcol);   // as.character(data[, 1])
   const std::vector<std::string> g   = groupToString(group);   // as.character(group)
   if (g.size() != X.n_rows)
     Rcpp::stop("group should have the same length as the number of rows of data!");

   std::vector<std::string> subjects, glev;
   std::vector<std::vector<arma::uword>> sidx, gidx;
   splitByGroup(ids, subjects, sidx);                       // unique(data[, 1])
   splitByGroup(g,   glev,     gidx);                       // unique(group)

   // The R version draws simTau (l values) but never uses it. Consume the same
   // draws so that results under set.seed() are identical. Safe to delete.
   for (int i = 0; i < l; ++i) (void)R::rexp(1.0 / alpha);

   const std::size_t nS = subjects.size();
   Rcpp::NumericMatrix Tau(nS, 1);
   Rcpp::List A(nS);                                        // per subject: named list of matrices
   std::vector<std::vector<std::string>> Alev(nS);          // group labels of each A[[s]]

   for (std::size_t s = 0; s < nS; ++s) {
     Rcpp::checkUserInterrupt();
     const arma::uvec rows = arma::conv_to<arma::uvec>::from(sidx[s]);
     const arma::mat  Xs   = X.rows(rows);
     std::vector<std::string> gs;
     gs.reserve(sidx[s].size());
     for (arma::uword r : sidx[s]) gs.push_back(g[r]);

     const TauPrep pr = prepareTau(Xs, gs, wi, expFix);

     // importanceSample()
     std::vector<double> samp, w;
     importanceSampleCore(l, pr, alpha, samp, w);

     // importanceEstimates()
     const double est = arma::dot(arma::conv_to<arma::vec>::from(samp),
                                  arma::conv_to<arma::vec>::from(w));
     Tau(s, 0) = est;

     Rcpp::List phi(pr.lev.size());
     Rcpp::CharacterVector nms(pr.lev.size());
     for (std::size_t k = 0; k < pr.lev.size(); ++k) {
       phi[k] = Rcpp::wrap(arma::mat(arma::exp(-est * pr.P[k])));   // == phifunction(t, est, expFix)
       nms[k] = pr.lev[k];
     }
     phi.names() = nms;
     A[s] = phi;
     Alev[s] = pr.lev;
   }

   Rcpp::CharacterVector subjNames(subjects.begin(), subjects.end());
   Tau.attr("dimnames") = Rcpp::List::create(subjNames, R_NilValue);   // rownames = subjects

   // Rearrange: AA[[group]][[subject]] = A[[subject]][[group]]
   std::unordered_map<std::string, std::size_t> subjPos;
   for (std::size_t s = 0; s < nS; ++s) subjPos.emplace(subjects[s], s);

   Rcpp::List AA(glev.size());
   for (std::size_t k = 0; k < glev.size(); ++k) {
     std::vector<std::size_t> subs;                 // subjects in group k, order of first appearance
     std::unordered_set<std::size_t> seen;
     for (arma::uword r : gidx[k]) {
       const std::size_t si = subjPos[ids[r]];
       if (seen.insert(si).second) subs.push_back(si);
     }

     Rcpp::List Ak(subs.size());
     Rcpp::CharacterVector nk(subs.size());
     for (std::size_t j = 0; j < subs.size(); ++j) {
       const std::size_t si = subs[j];
       nk[j] = subjects[si];
       auto it = std::find(Alev[si].begin(), Alev[si].end(), glev[k]);
       if (it != Alev[si].end()) {                  // otherwise stays NULL, as in R
         Rcpp::List Ai = A[si];
         SEXP m = Ai[it - Alev[si].begin()];
         Ak[j] = m;
       }
     }
     Ak.names() = nk;
     AA[k] = Ak;
   }
   AA.names() = Rcpp::CharacterVector(glev.begin(), glev.end());

   return Rcpp::List::create(Rcpp::Named("Tau") = Tau,
                             Rcpp::Named("AA")  = AA);
 }



// =============================================================================
// simulate_randomTau
// =============================================================================
#include <numeric>          // std::iota
#include <R_ext/Random.h>   // R_unif_index

// Random symmetric 0/1 structure, as in the R code:
//   candidates = strict upper triangle in column-major order
//   (= which(real_stru == 0, arr.ind = TRUE)), choose m without replacement
//   (same algorithm as R's sample()), return S + t(S) + diag(p).
static arma::mat randomStructure(int p, int m, const char* argname) {
  std::vector<arma::uword> rr, cc;
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < j; ++i) { rr.push_back(i); cc.push_back(j); }
    const int N = static_cast<int>(rr.size());
  if (m < 0 || m > N)
    Rcpp::stop("'%s' must be between 0 and p*(p-1)/2 = %d", argname, N);

  std::vector<int> x(N);
  std::iota(x.begin(), x.end(), 0);
  arma::mat S(p, p, arma::fill::zeros);
  int nn = N;
  for (int k = 0; k < m; ++k) {
    const int j   = static_cast<int>(R_unif_index(static_cast<double>(nn)));
    const int pos = x[j];
    x[j] = x[--nn];
    S(rr[pos], cc[pos]) = 1.0;
  }
  return S + S.t() + arma::eye<arma::mat>(p, p);
}

// A with A A' = S, eigen-based like MASS::mvrnorm (negative eigenvalues -> 0)
static arma::mat sqrtFactor(const arma::mat& S) {
  arma::vec ev;
  arma::mat V;
  if (!arma::eig_sym(ev, V, arma::mat(0.5 * (S + S.t()))))
    Rcpp::stop("eigen decomposition failed");
  const double tol = 1e-6 * std::abs(ev.max());
  if (arma::any(ev < -tol)) Rcpp::stop("'Sigma' is not positive definite");
  arma::rowvec s(ev.n_elem);
  for (arma::uword k = 0; k < ev.n_elem; ++k) s[k] = ev[k] > 0 ? std::sqrt(ev[k]) : 0.0;
  return V.each_row() % s;
}

// ---- networks, precision and covariance matrices (shared by the simulators) ----
struct SimNetworks {
  arma::mat stru1, stru2;     // adjacency matrices (pre, post)
  arma::mat sigma1, sigma2;   // covariance matrices (pre, post)
};

// Same steps and RNG order as the R code:
// m1 edges -> m2 differences -> theta ~ N(0, 2^2) -> MakePositiveDefinite -> solve
static SimNetworks simulateNetworks(int p, int m1, int m2) {
  SimNetworks out;
  out.stru1 = randomStructure(p, m1, "m1");
  const arma::mat disturbance = randomStructure(p, m2, "m2");
  out.stru2 = out.stru1 + disturbance;
  out.stru2.transform([](double v) { return std::fmod(v, 2.0); });   // %% 2

  Rcpp::NumericVector th = Rcpp::rnorm(p * p, 0.0, 2.0);
  arma::mat theta(th.begin(), p, p);
  for (int j = 0; j < p; ++j)
    for (int i = j; i < p; ++i) theta(i, j) = 0.0;      // lower triangle incl. diag = 0
  theta = theta + theta.t() + arma::eye<arma::mat>(p, p);

  arma::mat theta1 = theta % out.stru1;
  arma::mat theta2 = theta % out.stru2;

  Rcpp::Environment fakeNS = Rcpp::Environment::namespace_env("fake");
  Rcpp::Function makePD = fakeNS["MakePositiveDefinite"];
  Rcpp::List pd1 = makePD(Rcpp::Named("omega") = theta1,
                          Rcpp::Named("pd_strategy") = "diagonally_dominant",
                          Rcpp::Named("scale") = true);
  Rcpp::List pd2 = makePD(Rcpp::Named("omega") = theta2,
                          Rcpp::Named("pd_strategy") = "diagonally_dominant",
                          Rcpp::Named("scale") = true);
  theta1 = Rcpp::as<arma::mat>(pd1["omega"]);
  theta2 = Rcpp::as<arma::mat>(pd2["omega"]);

  if (!arma::inv(out.sigma1, theta1)) Rcpp::stop("theta1 is singular");
  if (!arma::inv(out.sigma2, theta2)) Rcpp::stop("theta2 is singular");
  return out;
}

// Adjacency matrix with dimnames V1..Vp (= feature names of the simulated data)
static Rcpp::NumericMatrix namedStructure(const arma::mat& M) {
  const int p = M.n_rows;
  Rcpp::NumericMatrix out(Rcpp::wrap(M));
  Rcpp::CharacterVector nm(p);
  for (int k = 0; k < p; ++k) nm[k] = "V" + std::to_string(k + 1);
  out.attr("dimnames") = Rcpp::List::create(nm, nm);
  return out;
}

// One draw from N(0, kronecker(Phi, Sigma)) returned as m x p (row = time point).
// Y = A Z B' with A A' = Sigma, B B' = Phi gives cov(vec(Y)) = Phi %x% Sigma
// without forming the (m*p) x (m*p) matrix.
static arma::mat simSubject(const arma::mat& A, const arma::mat& Phi) {
  const arma::uword p = A.n_rows, m = Phi.n_rows;
  const arma::mat B = sqrtFactor(Phi);
  Rcpp::NumericVector z = Rcpp::rnorm(p * m);
  const arma::mat Z(z.begin(), p, m);          // copies
  return (A * Z * B.t()).t();
}

// data.frame(subject, time, V1, ..., Vp), same layout as the R version
static Rcpp::List buildSimDataFrame(const std::vector<arma::vec>& tp,
                                    const std::vector<arma::mat>& Y, int p) {
  const int n = static_cast<int>(tp.size());
  int N = 0;
  for (int i = 0; i < n; ++i) N += tp[i].n_elem;

  Rcpp::CharacterVector subject(N);
  Rcpp::NumericVector   time(N);
  arma::mat V(N, p);
  int r = 0;
  for (int i = 0; i < n; ++i) {
    const std::string lab = "subject" + std::to_string(i + 1);
    const int mi = tp[i].n_elem;
    for (int j = 0; j < mi; ++j) {
      subject[r + j] = lab;
      time[r + j]    = tp[i][j];
    }
    V.rows(r, r + mi - 1) = Y[i];
    r += mi;
  }

  Rcpp::List out(p + 2);
  Rcpp::CharacterVector nm(p + 2);
  out[0] = subject; nm[0] = "subject";
  out[1] = time;    nm[1] = "time";
  for (int k = 0; k < p; ++k) {
    out[k + 2] = Rcpp::NumericVector(V.colptr(k), V.colptr(k) + N);
    nm[k + 2]  = "V" + std::to_string(k + 1);
  }
  out.attr("names")     = nm;
  out.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -N);
  out.attr("class")     = "data.frame";
  return out;
}

//' Simulate data with random tau
//' @param n the number of subjects
//' @param p the dimension of the normal distribution
//' @param m1 the number of edges
//' @param tt the average length of data for each subject
//' @param m2 the number of differences between two networks
//' @param alpha the parameter in the exponential distribution
//' @param group 1 = pre only, otherwise pre and post
//' @noRd
//' @returns a list with data, network, tau and alpha
 // [[Rcpp::export]]
 Rcpp::List simulate_randomTau(int n, int p, int m1, int tt, int m2,
                               double alpha, int group) {
   if (n < 1)  Rcpp::stop("n must be >= 1");
   if (p < 1)  Rcpp::stop("p must be >= 1");
   if (tt < 1) Rcpp::stop("tt must be >= 1");
   if (!(alpha > 0) || !std::isfinite(alpha)) Rcpp::stop("alpha must be a positive finite number");
   const bool twoGroups = (group != 1);          // R: if (group == 1) ... else ...

   // ---- 1. true tau, time points, correlation matrices ----------------------
   Rcpp::NumericVector trueTau = Rcpp::rexp(n, alpha);
   std::vector<arma::vec> timepoint1(n), timepoint2(n);
   std::vector<arma::mat> cc1(n), cc2(n);

   for (int i = 0; i < n; ++i) {
     const int m3  = 1 + static_cast<int>(R_unif_index(static_cast<double>(tt)));  // sample(1:tt, 1)
     const int len = twoGroups ? 2 * m3 : m3;
     Rcpp::NumericVector e = Rcpp::rexp(len, 1.0);                                 // stats::rexp(len)
     const arma::vec cs = arma::cumsum(arma::vec(e.begin(), len));

     timepoint1[i] = cs.head(m3);
     cc1[i] = phifunction(timepoint1[i], trueTau[i], 1.0);
     if (twoGroups) {
       timepoint2[i] = cs.subvec(m3, 2 * m3 - 1);
       cc2[i] = phifunction(timepoint2[i], trueTau[i], 1.0);
     }
   }

   // ---- 2. network structures -----------------------------------------------
   const arma::mat real_stru1  = randomStructure(p, m1, "m1");
   const arma::mat disturbance = randomStructure(p, m2, "m2");
   arma::mat real_stru2 = real_stru1 + disturbance;
   real_stru2.transform([](double v) { return std::fmod(v, 2.0); });     // %% 2

   // ---- 3. precision and covariance matrices ---------------------------------
   Rcpp::NumericVector th = Rcpp::rnorm(p * p, 0.0, 2.0);
   arma::mat theta(th.begin(), p, p);
   for (int j = 0; j < p; ++j)
     for (int i = j; i < p; ++i) theta(i, j) = 0.0;     // lower triangle incl. diag = 0
   theta = theta + theta.t() + arma::eye<arma::mat>(p, p);

   arma::mat theta1 = theta % real_stru1;
   arma::mat theta2 = theta % real_stru2;

   Rcpp::Environment fakeNS = Rcpp::Environment::namespace_env("fake");
   Rcpp::Function makePD = fakeNS["MakePositiveDefinite"];
   Rcpp::List pd1 = makePD(Rcpp::Named("omega") = theta1,
                           Rcpp::Named("pd_strategy") = "diagonally_dominant",
                           Rcpp::Named("scale") = true);
   Rcpp::List pd2 = makePD(Rcpp::Named("omega") = theta2,
                           Rcpp::Named("pd_strategy") = "diagonally_dominant",
                           Rcpp::Named("scale") = true);
   theta1 = Rcpp::as<arma::mat>(pd1["omega"]);
   theta2 = Rcpp::as<arma::mat>(pd2["omega"]);

   arma::mat sigma1, sigma2;
   if (!arma::inv(sigma1, theta1)) Rcpp::stop("theta1 is singular");
   if (!arma::inv(sigma2, theta2)) Rcpp::stop("theta2 is singular");

   // ---- 4. simulate data: pre for all subjects, then post ---------------------
   const arma::mat A1 = sqrtFactor(sigma1);
   std::vector<arma::mat> Y1(n);
   for (int i = 0; i < n; ++i) Y1[i] = simSubject(A1, cc1[i]);
   Rcpp::List a1 = buildSimDataFrame(timepoint1, Y1, p);

   if (twoGroups) {
     const arma::mat A2 = sqrtFactor(sigma2);
     std::vector<arma::mat> Y2(n);
     for (int i = 0; i < n; ++i) Y2[i] = simSubject(A2, cc2[i]);
     Rcpp::List a2 = buildSimDataFrame(timepoint2, Y2, p);

     return Rcpp::List::create(
       Rcpp::Named("data")    = Rcpp::List::create(Rcpp::Named("pre") = a1,
                   Rcpp::Named("post") = a2),
                   Rcpp::Named("network") = Rcpp::List::create(Rcpp::Named("pre") = real_stru1,
                               Rcpp::Named("post") = real_stru2),
                               Rcpp::Named("tau")     = trueTau,
                               Rcpp::Named("alpha")   = alpha);
   }

   return Rcpp::List::create(
     Rcpp::Named("data")    = Rcpp::List::create(Rcpp::Named("pre") = a1),
     Rcpp::Named("network") = Rcpp::List::create(Rcpp::Named("pre") = real_stru1),
     Rcpp::Named("tau")     = trueTau,
     Rcpp::Named("alpha")   = alpha);
 }



//' Simulate longitudinal data with fixed tau
//'
//' @param n the number of subjects
//' @param p the dimension of the normal distribution
//' @param m1 the number of edges
//' @param tau the dampening rate: length 1 (pre only) or 2 (pre and post)
//' @param tt the maximum number of time points per subject
//' @param m2 the number of differences between the two networks
//' @noRd
//' @returns a list with data, network and tau
 // [[Rcpp::export]]
 Rcpp::List simulate_long(int n, int p, int m1, Rcpp::NumericVector tau,
                          int tt = 5, int m2 = 0) {
   if (n < 1)  Rcpp::stop("n must be >= 1");
   if (p < 1)  Rcpp::stop("p must be >= 1");
   if (tt < 1) Rcpp::stop("tt must be >= 1");
   if (tau.size() != 1 && tau.size() != 2) Rcpp::stop("tau must have length 1 or 2");
   for (R_xlen_t k = 0; k < tau.size(); ++k)
     if (!(tau[k] > 0) || !std::isfinite(tau[k]))
       Rcpp::stop("tau must contain positive finite numbers");
     const bool twoGroups = (tau.size() == 2);

     // ---- 1. time points and correlation matrices -----------------------------
     std::vector<arma::vec> timepoint1(n), timepoint2(n);
     std::vector<arma::mat> cc1(n), cc2(n);
     for (int i = 0; i < n; ++i) {
       const int m3  = 1 + static_cast<int>(R_unif_index(static_cast<double>(tt)));  // sample(1:tt, 1)
       const int len = twoGroups ? 2 * m3 : m3;
       Rcpp::NumericVector e = Rcpp::rexp(len, 1.0);                                 // stats::rexp(len)
       const arma::vec cs = arma::cumsum(arma::vec(e.begin(), len));

       timepoint1[i] = cs.head(m3);
       if (twoGroups) timepoint2[i] = cs.subvec(m3, 2 * m3 - 1);
     }
     for (int i = 0; i < n; ++i) {
       cc1[i] = phifunction(timepoint1[i], tau[0], 1.0);
       if (twoGroups) cc2[i] = phifunction(timepoint2[i], tau[1], 1.0);
     }

     // ---- 2. networks and covariance matrices ---------------------------------
     const SimNetworks net = simulateNetworks(p, m1, m2);

     // ---- 3. data: pre for all subjects, then post -----------------------------
     const arma::mat A1 = sqrtFactor(net.sigma1);
     std::vector<arma::mat> Y1(n);
     for (int i = 0; i < n; ++i) Y1[i] = simSubject(A1, cc1[i]);
     Rcpp::List a1 = buildSimDataFrame(timepoint1, Y1, p);

     Rcpp::NumericMatrix real_stru1 = namedStructure(net.stru1);

     if (twoGroups) {
       const arma::mat A2 = sqrtFactor(net.sigma2);
       std::vector<arma::mat> Y2(n);
       for (int i = 0; i < n; ++i) Y2[i] = simSubject(A2, cc2[i]);
       Rcpp::List a2 = buildSimDataFrame(timepoint2, Y2, p);

       return Rcpp::List::create(
         Rcpp::Named("data")    = Rcpp::List::create(Rcpp::Named("pre") = a1,
                     Rcpp::Named("post") = a2),
                     Rcpp::Named("network") = Rcpp::List::create(Rcpp::Named("pre") = real_stru1,
                                 Rcpp::Named("post") = namedStructure(net.stru2)),
                                 Rcpp::Named("tau")     = tau);
     }

     // one tau: not nested, as in the R version
     return Rcpp::List::create(
       Rcpp::Named("data")    = a1,
       Rcpp::Named("network") = real_stru1,
       Rcpp::Named("tau")     = tau);
 }


// ---- helpers for BB --------------------------------------------------------

static inline double softThr(double x, double t) {
  return x > t ? x - t : (x < -t ? x + t : 0.0);
}

// log|M| for a symmetric matrix; NaN if det <= 0 (like log(det()) in R)
static double logDetSym(const arma::mat& M) {
  arma::mat L;
  if (arma::chol(L, M, "lower")) return 2.0 * arma::accu(arma::log(L.diag()));
  double val, sign;
  if (!arma::log_det(val, sign, M) || sign <= 0)
    return std::numeric_limits<double>::quiet_NaN();
  return val;
}

// colnames(x)[-c(1, 2)] for a data.frame or matrix; NULL if unavailable
static Rcpp::RObject featureNamesOf(SEXP x) {
  Rcpp::RObject nm;
  if (Rf_inherits(x, "data.frame")) {
    nm = Rf_getAttrib(x, R_NamesSymbol);
  } else if (Rf_isMatrix(x)) {
    Rcpp::RObject dn = Rf_getAttrib(x, R_DimNamesSymbol);
    if (!dn.isNULL()) nm = VECTOR_ELT(dn, 1);
  }
  if (nm.isNULL()) return R_NilValue;
  Rcpp::CharacterVector all(nm);
  if (all.size() < 3) return R_NilValue;
  return Rcpp::CharacterVector(all.begin() + 2, all.end());
}

// tau-independent pieces of one stage
struct BBGroup {
  arma::mat   Sraw;        // sum_j X_j' A_j^{-1} X_j   (amatrix[[i]] in R)
  double      extra = 0.0; // sum_j p * log|A_j|
  double      n = 0.0;     // nrow(dd)
  double      nn = 0.0;    // number of subjects
  arma::uword p = 0;
};

static BBGroup bbPrepare(SEXP dd, SEXP AiS, int stage) {
  if (TYPEOF(AiS) != VECSXP) Rcpp::stop("A[[%d]] must be a list of matrices", stage);
  Rcpp::List Ai(AiS);

  const arma::mat X = dataToMatrix(dd);              // col 0 = time, cols 1..p = variables
  if (!X.is_finite())
    Rcpp::stop("data[[%d]]: time and variable columns must not contain NA/NaN/Inf", stage);

  BBGroup G;
  G.p = X.n_cols - 1;
  G.n = static_cast<double>(X.n_rows);

  Rcpp::RObject idcol = firstColumn(dd);
  const std::vector<std::string> ids = groupToString(idcol);
  std::vector<std::string> subj;
  std::vector<std::vector<arma::uword>> idx;
  splitByGroup(ids, subj, idx);                      // unique(dd[, 1]) order
  G.nn = static_cast<double>(subj.size());

  if (static_cast<std::size_t>(Ai.size()) != subj.size())
    Rcpp::stop("Data do not match! A[[%d]] has %d matrices but data[[%d]] has %d subjects",
               stage, (int)Ai.size(), stage, (int)subj.size());

  // Match matrices to subjects by name if possible, otherwise by position
  std::vector<R_xlen_t> pos(subj.size());
  bool byName = false;
  if (!Rf_isNull(Ai.names())) {
    Rcpp::CharacterVector an = Ai.names();
    std::unordered_map<std::string, R_xlen_t> ap;
    for (R_xlen_t k = 0; k < an.size(); ++k) ap.emplace(Rcpp::as<std::string>(an[k]), k);
    byName = true;
    for (std::size_t s = 0; s < subj.size(); ++s) {
      auto it = ap.find(subj[s]);
      if (it == ap.end()) { byName = false; break; }
      pos[s] = it->second;
    }
  }
  if (!byName) for (std::size_t s = 0; s < subj.size(); ++s) pos[s] = (R_xlen_t)s;

  G.Sraw.zeros(G.p, G.p);
  const double pd = static_cast<double>(G.p);
  for (std::size_t s = 0; s < subj.size(); ++s) {
    SEXP aij = Ai[pos[s]];
    if (Rf_isNull(aij))
      Rcpp::stop("A[[%d]] has no matrix for subject '%s'", stage, subj[s]);
    const arma::mat Aij  = Rcpp::as<arma::mat>(aij);
    const arma::uvec rows = arma::conv_to<arma::uvec>::from(idx[s]);
    const arma::mat Xr   = X.rows(rows);
    const arma::mat Y    = Xr.cols(1, G.p);          // n_ij x p
    if (Aij.n_rows != Y.n_rows || Aij.n_cols != Y.n_rows)
      Rcpp::stop("The format of A does not match the format of data!! (stage %d, subject '%s')",
                 stage, subj[s]);

    arma::mat L;
    double logdet;
    if (arma::chol(L, Aij, "lower")) {
      logdet = 2.0 * arma::accu(arma::log(L.diag()));
      const arma::mat Z = arma::solve(arma::trimatl(L), Y);
      G.Sraw += Z.t() * Z;                           // Y' A^{-1} Y
    } else {
      arma::mat AY;
      if (!arma::solve(AY, Aij, Y, arma::solve_opts::no_approx))
        Rcpp::stop("correlation matrix of subject '%s' (stage %d) is singular", subj[s], stage);
      G.Sraw += Y.t() * AY;
      logdet = logDetSym(Aij);
    }
    G.extra += pd * logdet;
  }
  G.Sraw = 0.5 * (G.Sraw + G.Sraw.t());
  return G;
}

// Theta-step: argmin -log|T| + tr(S T) + rho/2 ||T - (Z - U)||_F^2  (closed form)
static arma::mat admmThetaUpdate(const arma::mat& S, const arma::mat& Z,
                                 const arma::mat& U, double rho) {
  arma::mat M = rho * (Z - U) - S;
  M = 0.5 * (M + M.t());
  arma::vec d;
  arma::mat V;
  if (!arma::eig_sym(d, V, M)) Rcpp::stop("eigen decomposition failed in ADMM");
  const arma::vec th = (d + arma::sqrt(arma::square(d) + 4.0 * rho)) / (2.0 * rho);
  return (V.each_row() % th.t()) * V.t();
}

struct ADMMControl {
  double rho;
  int    maxit;
  double epsAbs, epsRel;
};

// ADMM for the graphical lasso (K = 1) and the fused graphical lasso (K = 2).
// K = 1: penDiag decides whether the diagonal is penalized (glasso default: yes).
// K = 2: diagonal unpenalized, as in the CVXR masks.
static std::vector<arma::mat> admmGraphLasso(const std::vector<arma::mat>& S,
                                             double lam1, double lam2, bool penDiag,
                                             const ADMMControl& ctl,
                                             int& iterOut, bool& convOut) {
  const std::size_t K = S.size();
  const arma::uword p = S[0].n_rows;
  double rho = ctl.rho;

  std::vector<arma::mat> Th(K), Z(K), U(K), Zold(K);
  for (std::size_t k = 0; k < K; ++k) {
    arma::vec d = S[k].diag() + lam1;
    d.transform([](double v) { return v > 1e-8 ? 1.0 / v : 1.0; });
    Z[k] = arma::diagmat(d);                         // warm start
    U[k].zeros(p, p);
  }

  convOut = false;
  iterOut = 0;
  const double sqrtKp = std::sqrt(static_cast<double>(K)) * static_cast<double>(p);

  for (int it = 1; it <= ctl.maxit; ++it) {
    if (it % 50 == 0) Rcpp::checkUserInterrupt();

    // 1. Theta-step
    for (std::size_t k = 0; k < K; ++k) Th[k] = admmThetaUpdate(S[k], Z[k], U[k], rho);

    // 2. Z-step (elementwise, exact)
    for (std::size_t k = 0; k < K; ++k) Zold[k] = Z[k];
    const double t1 = lam1 / rho, t2 = lam2 / rho;
    if (K == 1) {
      const arma::mat Aa = Th[0] + U[0];
      for (arma::uword j = 0; j < p; ++j)
        for (arma::uword i = 0; i <= j; ++i) {
          const double v = (i == j && !penDiag) ? Aa(i, j) : softThr(Aa(i, j), t1);
          Z[0](i, j) = Z[0](j, i) = v;
        }
    } else {
      const arma::mat A1 = Th[0] + U[0], A2 = Th[1] + U[1];
      for (arma::uword j = 0; j < p; ++j)
        for (arma::uword i = 0; i <= j; ++i) {
          const double a1 = A1(i, j), a2 = A2(i, j);
          double z1, z2;
          if (i == j) {                              // diagonal: no penalty
            z1 = a1; z2 = a2;
          } else {                                   // fused step, then soft-threshold
            if (a1 > a2 + 2.0 * t2)      { z1 = a1 - t2; z2 = a2 + t2; }
            else if (a2 > a1 + 2.0 * t2) { z1 = a1 + t2; z2 = a2 - t2; }
            else                         { z1 = z2 = 0.5 * (a1 + a2); }
            z1 = softThr(z1, t1);
            z2 = softThr(z2, t1);
          }
          Z[0](i, j) = Z[0](j, i) = z1;
          Z[1](i, j) = Z[1](j, i) = z2;
        }
    }

    // 3. dual update and residuals
    double r2 = 0.0, s2 = 0.0, nTh = 0.0, nZ = 0.0, nU = 0.0;
    for (std::size_t k = 0; k < K; ++k) {
      const arma::mat D = Th[k] - Z[k];
      U[k] += D;
      r2  += arma::accu(arma::square(D));
      s2  += arma::accu(arma::square(Z[k] - Zold[k]));
      nTh += arma::accu(arma::square(Th[k]));
      nZ  += arma::accu(arma::square(Z[k]));
      nU  += arma::accu(arma::square(U[k]));
    }
    const double r = std::sqrt(r2), s = rho * std::sqrt(s2);
    const double epsPri  = sqrtKp * ctl.epsAbs + ctl.epsRel * std::max(std::sqrt(nTh), std::sqrt(nZ));
    const double epsDual = sqrtKp * ctl.epsAbs + ctl.epsRel * rho * std::sqrt(nU);
    iterOut = it;
    if (r <= epsPri && s <= epsDual) { convOut = true; break; }

    // 4. adaptive rho (residual balancing); scaled dual must be rescaled
    if (r > 10.0 * s)      { rho *= 2.0; for (auto& u : U) u /= 2.0; }
    else if (s > 10.0 * r) { rho /= 2.0; for (auto& u : U) u *= 2.0; }
  }

  // Return the sparse iterate Z if it is positive definite, otherwise Theta
  std::vector<arma::mat> out(K);
  for (std::size_t k = 0; k < K; ++k) {
    arma::mat L;
    out[k] = arma::chol(L, Z[k]) ? Z[k] : Th[k];
    out[k] = 0.5 * (out[k] + out[k].t());
  }
  return out;
}



//' Estimate the precision matrices (ADMM)
//'
//' @param A a list of length 1 or 2 (stages). Each entry is a list of the phi
//'   matrices of all subjects in that stage.
//' @param data a list of data frames (col 1 = id, col 2 = time, rest = variables)
//' @param lambda tuning parameter(s): lambda(1) sparsity, lambda(2) fusion (m = 2)
//' @param random logical, type of the model
//' @param tau needed if random = TRUE (enters the likelihood via log(mean(tau)))
//' @param rho initial ADMM step size (adapted automatically)
//' @param maxit maximum number of ADMM iterations
//' @param eps_abs absolute ADMM tolerance
//' @param eps_rel relative ADMM tolerance
//' @noRd
//' @returns a list with wi (precision matrices), w (covariance matrices) and ll
 // [[Rcpp::export]]
 Rcpp::List BB(SEXP A, SEXP data, Rcpp::NumericVector lambda, bool random = false,
               Rcpp::Nullable<Rcpp::NumericVector> tau = R_NilValue,
               double rho = 1.0, int maxit = 10000,
               double eps_abs = 1e-6, double eps_rel = 1e-5) {
   if (TYPEOF(A) != VECSXP || TYPEOF(data) != VECSXP) Rcpp::stop("A and data must be lists!");
   if (Rf_inherits(data, "data.frame"))
     Rcpp::stop("data must be a list of data frames, not a single data frame");
   Rcpp::List Al(A), Dl(data);
   const R_xlen_t m = Al.size();
   if (m != Dl.size()) Rcpp::stop(" List A should have same length as list  data!");
   if (m < 1 || m > 2) Rcpp::stop("A and data must have length 1 or 2");

   const R_xlen_t nlam = (m == 2) ? 2 : 1;
   if (lambda.size() < nlam) Rcpp::stop("lambda must have length %d", (int)nlam);
   for (R_xlen_t k = 0; k < nlam; ++k)
     if (!std::isfinite(lambda[k]) || lambda[k] < 0)
       Rcpp::stop("lambda must contain non-negative finite numbers");
     if (!(rho > 0) || !std::isfinite(rho)) Rcpp::stop("rho must be a positive finite number");
     if (maxit < 1) Rcpp::stop("maxit must be >= 1");

     double logMeanTau = 0.0;
     if (random) {
       if (tau.isNull()) Rcpp::stop("tau must be supplied when random = TRUE");
       Rcpp::NumericVector tv(tau.get());
       if (tv.size() == 0) Rcpp::stop("tau must not be empty");
       const double mt = Rcpp::mean(tv);
       if (!(mt > 0) || !std::isfinite(mt)) Rcpp::stop("mean(tau) must be a positive finite number");
       logMeanTau = std::log(mt);
     }

     Rcpp::RObject glev = Rf_getAttrib(A, R_NamesSymbol);
     Rcpp::RObject featureNames = featureNamesOf(Dl[0]);

     // ---- tau-independent pieces per stage --------------------------------------
     std::vector<BBGroup> G;
     G.reserve(m);
     for (R_xlen_t i = 0; i < m; ++i) G.push_back(bbPrepare(Dl[i], Al[i], (int)(i + 1)));
     const arma::uword p = G[0].p;
     if (m == 2 && G[1].p != p) Rcpp::stop("all data sets must have the same number of variables");
     if (!featureNames.isNULL() && Rf_xlength(featureNames) != (R_xlen_t)p) featureNames = R_NilValue;

     std::vector<arma::mat> S(m);
     for (R_xlen_t i = 0; i < m; ++i) S[i] = G[i].Sraw / G[i].n;

     // ---- ADMM ------------------------------------------------------------------
     const ADMMControl ctl{rho, maxit, eps_abs, eps_rel};
     int iters = 0;
     bool conv = false;
     const std::vector<arma::mat> est =
       admmGraphLasso(S, lambda[0], (m == 2) ? lambda[1] : 0.0, /*penDiag=*/ m == 1,
                      ctl, iters, conv);
     if (!conv)
       Rcpp::warning("ADMM did not converge in %d iterations; consider increasing maxit", maxit);

     // ---- likelihood and output -------------------------------------------------
     double likelihood = 0.0;
     Rcpp::List wi(m), w(m);
     for (R_xlen_t i = 0; i < m; ++i) {
       const arma::mat& Th = est[i];
       likelihood += -G[i].extra / G[i].n + logDetSym(Th) - arma::accu(Th % G[i].Sraw) / G[i].n;
       if (random) likelihood -= 2.0 * G[i].nn * logMeanTau / G[i].n;

       arma::mat Wi;
       if (!arma::inv_sympd(Wi, Th) && !arma::inv(Wi, Th))
         Rcpp::stop("estimated precision matrix %d is singular", (int)(i + 1));

       Rcpp::NumericMatrix thR(Rcpp::wrap(Th)), wR(Rcpp::wrap(Wi));
       if (!featureNames.isNULL()) {
         thR.attr("dimnames") = Rcpp::List::create(featureNames, featureNames);
         wR.attr("dimnames")  = Rcpp::List::create(featureNames, featureNames);
       }
       wi[i] = thR;
       w[i]  = wR;
     }
     if (!glev.isNULL()) wi.attr("names") = glev;

     return Rcpp::List::create(Rcpp::Named("wi") = wi,
                               Rcpp::Named("w")  = w,
                               Rcpp::Named("ll") = likelihood / 2.0);
 }
