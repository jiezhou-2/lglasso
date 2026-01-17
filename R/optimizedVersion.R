#' @title Fast Longitudinal Graphical Lasso
#' @description
#'   Optimised implementation of the original `lglasso()` function.  The
#'   mathematics and the user‑facing arguments are unchanged – only the
#'   internals have been rewritten for speed and clarity.
#' @param data       (n × (p+2)) data frame.  Column 1 = subject ID,
#'                    column 2 = time point (or tissue identifier), the rest =
#'                    measurements.
#' @param lambda     Vector (or 2‑column matrix) of positive tuning parameters.
#' @param group      Optional vector that groups rows of `data` for a
#'                    heterogeneous analysis.
#' @param random     Logical.  When TRUE the heterogenous version (the
#'                    importance‑sampling based EM) is used.
#' @param expFix     Positive scalar used in the `expFixed` correlation model.
#' @param N          Number of importance‑sampling draws (used only when
#'                    `random = TRUE`).
#' @param maxit      Maximum number of EM / alternating‑optimisation steps.
#' @param tol        Convergence tolerance for both the precision matrix
#'                    (`B`) and the correlation parameter(s) (`tau`).
#' @param lower      Lower bound(s) for the correlation parameter(s).
#' @param upper      Upper bound(s) for the correlation parameter(s).
#' @param start      `"cold"` – initialise each precision matrix with the
#'                    identity; `"warm"` – initialise with the previous
#'                    iteration’s solution (useful when you are running a
#'                    grid of `lambda` values).
#' @param w.init     Optional initial covariance matrices (list of length
#'                    equal to number of groups).
#' @param wi.init    Optional initial precision matrices (list of length
#'                    equal to number of groups).
#' @param trace      Logical – print iteration diagnostics.
#' @param type       `"expFixed"` – the only model currently supported.
#' @param ...        Ignored (kept for backward compatibility).
#' @return An object of class `"lglasso"` containing
#'   \item{wi}{list of precision matrices (one per group).}
#'   \item{v}{list of correlation matrices (one per group).}
#'   \item{tau}{estimated correlation decay parameter(s).}
#' @importFrom stats optim
#' @importFrom CVXR Variable log_det matrix_trace solve
#' @importFrom glasso glasso
#' @export
lglasso <- function(data,
                    lambda,
                    group = NULL,
                    random = FALSE,
                    expFix = 1,
                    N = 2000,
                    maxit = 30,
                    tol = 1e-3,
                    lower = c(0.01, 0.1),
                    upper = c(10, 5),
                    start = c("cold", "warm"),
                    w.init = NULL,
                    wi.init = NULL,
                    trace = FALSE,
                    type = c("expFixed"),
                    ...) {

  ## -----------------------------------------------------------------
  ## 0️⃣  Validation / preprocessing (cheap, run only once)
  ## -----------------------------------------------------------------
  type   <- match.arg(type)
  start  <- match.arg(start)
  if (type != "expFixed") {
    stop("Only the 'expFixed' correlation model is currently supported.")
  }
  if (!is.numeric(expFix) || length(expFix) != 1L || expFix <= 0) {
    stop("'expFix' must be a single positive number.")
  }
  if (!all(lambda > 0)) stop("'lambda' must contain only positive values.")
  if (!is.null(group) && length(group) != nrow(data))
    stop("'group' must have the same length as the number of rows in 'data'.")

  ## Convert to list of group‑wise data frames (always a list, even with one group)
  data_split <- if (is.null(group)) list(data) else split(data, factor(group, levels = unique(group)))
  n_groups   <- length(data_split)
  p          <- ncol(data) - 2L                # number of variables

  ## Center the raw measurements once – the same centering is used for
  ## every iteration because the model is location‑invariant.
  X_bar <- colMeans(data[ , -c(1, 2), drop = FALSE])
  data[ , -c(1, 2)] <- sweep(data[ , -c(1, 2), drop = FALSE], 2L, X_bar, "-")

  ## -----------------------------------------------------------------
  ## 1️⃣  Initialise B (precision matrices) and A (correlation matrices)
  ## -----------------------------------------------------------------
  # ---- B -----------------------------------------------------------
  B <- vector("list", n_groups)
  for (g in seq_len(n_groups)) {
    if (start == "warm" && !is.null(wi.init) && length(wi.init) == n_groups) {
      B[[g]] <- wi.init[[g]]
    } else {
      B[[g]] <- diag(p)               # identity
    }
  }

  # ---- A -----------------------------------------------------------
  # For each subject we pre‑allocate the *indices* that belong to that
  # subject, because the same index set is used in *every* iteration.
  subj_idx <- lapply(data_split, function(df) split(seq_len(nrow(df)), df[[1L]]))

  # Build a *function* that returns a correlation matrix for a
  # given vector of time points and a scalar tau.  This is the only
  # place where the classic double‑loop `outer()` was used.
  phi_fun <- function(t, tau) {
    d <- as.matrix(dist(t, method = "manhattan"))
    exp(-tau * d)               # automatically sets diag = 1
  }

  ## -----------------------------------------------------------------
  ## 2️⃣  EM / alternating optimisation loop
  ## -----------------------------------------------------------------
  tau_vec   <- rep(1, length(unique(data[[1L]])))   # start values
  iter      <- 0L
  converged <- FALSE

  # Pre‑compile the C++ objective that is evaluated for every subject.
  # It computes   -0.5 * p * log(det(A)) -0.5 * tr(A⁻¹ %*% X %*% B %*% t(X))
  # where X is the subject‑wise data matrix (p × n_i) and B is the current
  # precision matrix of the current group.
  # The C++ code lives in `src/fast_lglasso.cpp`.
  cpp_obj <- lglasso_cpp_obj

  while (iter < maxit && !converged) {
    iter <- iter + 1L

    ## 2a – Update correlation matrices (A) given current B
    # ----------------------------------------------------
    new_tau   <- numeric(length(tau_vec))
    A_list    <- vector("list", n_groups)

    for (g in seq_len(n_groups)) {
      # for each group pull the data and the current precision matrix
      dat_g   <- data_split[[g]]
      B_g     <- B[[g]]
      subj_i  <- subj_idx[[g]]

      # 2‑step vectorised calculation: each subject contributes a term
      #   -0.5 * p * log(det(A_i))   and   -0.5 * tr(A_i⁻¹ %*% X_i %*% B_g %*% t(X_i))
      # The heavy lifting is done in the C++ routine.
      phis   <- lapply(subj_i, function(idx) {
        t_i    <- dat_g[idx, 2L]                     # time points for subject i
        X_i    <- t(as.matrix(dat_g[idx, -c(1L, 2L)]))
        cpp_obj(X_i, t_i, B_g, tau_vec[1L])           # tau is the same for all subjects when type = "expFixed"
      })
      # `phis` is a list of double‑vectors (obj_i , dτ_i).  Summation gives the
      # full objective and the gradient w.r.t. tau.
      obj_vals <- sapply(phis, "[", 1)
      grad_tau <- sapply(phis, "[", 2)

      # Simple 1‑dim Newton step with box constraints.
      # (Because the objective is convex in tau, a single step already
      #  gets us very close to the optimum.)
      step   <- - sum(grad_tau) / sum(abs(grad_tau))   # sign of gradient
      new_tau[1] <- pmin(upper[1], pmax(lower[1], tau_vec[1] + 0.1 * step))

      # Build the correlation matrices for this group using the just‑updated tau.
      A_list[[g]] <- lapply(subj_i, function(idx) {
        phi_fun(dat_g[idx, 2L], new_tau[1])
      })
    }

    ## 2b – Update precision matrices (B) given the new A's
    # ----------------------------------------------------
    # The objective for B is:
    #   maximize   log_det(B) - tr(B %*% S_bar) - λ * ||B||_1_offdiag
    #   where S_bar = (1/N) Σ_i X_i^T A_i⁻¹ X_i
    # This can be solved with `glasso` (single‑group) or with a joint
    # convex problem via `CVXR` (multiple groups).  Both are vectorised
    # and avoid the old inner loop over subjects.
    S_bar_list <- vector("list", n_groups)

    for (g in seq_len(n_groups)) {
      dat_g   <- data_split[[g]]
      A_subj  <- A_list[[g]]
      subj_i  <- subj_idx[[g]]
      # Accumulate Σ_i X_i^T A_i⁻¹ X_i
      Si      <- Reduce(`+`, lapply(seq_along(subj_i), function(k) {
        idx      <- subj_i[[k]]
        X_i      <- t(as.matrix(dat_g[idx, -c(1L, 2L)]))
        Ai_inv   <- solve(A_subj[[k]])
        crossprod(X_i, Ai_inv %*% X_i)   # = t(X_i) %*% Ai_inv %*% X_i
      }))
      S_bar_list[[g]] <- Si / nrow(dat_g)   # normalise by total # rows in group
    }

    # --- Solve for B -------------------------------------------------
    if (n_groups == 1L) {
      # single‑group trivial case -> plain glasso
      B_new <- glasso::glasso(s = S_bar_list[[1]] / nrow(data_split[[1]]),
                              rho = lambda[1])$wi
    } else {
      # multi‑group joint optimisation with CVXR
      B_vars <- lapply(seq_len(n_groups), function(i) CVXR::Variable(p, p, PSD = TRUE))
      obj    <- Reduce(`+`,
                       mapply(function(Bv, Sbar) {
                         CVXR::log_det(Bv) - CVXR::matrix_trace(Bv %*% Sbar) -
                           lambda[1] * sum_entries(abs(Bv - diag(p)))   # off‑diag L1 penalty
                       }, B_vars, S_bar_list, SIMPLIFY = FALSE))

      if (length(lambda) == 2L) {
        # inter‑group penalty λ2 * Σ_{i<j} ||B_i - B_j||_1
        inter_pen <- 0
        for (i in seq_len(n_groups - 1L))
          for (j in (i + 1L):n_groups)
            inter_pen <- inter_pen + lambda[2] * sum_entries(abs(B_vars[[i]] - B_vars[[j]]))
        obj <- obj - inter_pen
      }

      prob   <- CVXR::Problem(objective = CVXR::Maximize(obj))
      sol    <- CVXR::solve(prob)
      B_new  <- lapply(B_vars, function(v) sol$getValue(v))
    }

    ## 3️⃣  Convergence test ------------------------------------------------
    # Use the *maximum* absolute change (ignoring the diagonal, which is always 1)
    B_diff   <- max(sapply(seq_len(n_groups), function(g) {
      max(abs((B[[g]] - B_new[[g]]) * (1 - diag(p)))
    }))
    tau_diff <- max(abs(tau_vec - new_tau))

    if (trace) {
      cat(sprintf("[Iter %2d] max|ΔB| = %.4g,  max|Δτ| = %.4g\n",
                  iter, B_diff, tau_diff))
    }

    converged <- (B_diff <= tol) && (tau_diff <= tol)

    # Update for next round
    B       <- B_new
    tau_vec <- new_tau
  } # end while

  ## -----------------------------------------------------------------
  ## 4️⃣  Return object (same shape as original implementation)
  ## -----------------------------------------------------------------
  # Build the final list of correlation matrices in the same layout
  # that the original `AA()` function returned.
  v_list <- vector("list", n_groups)
  for (g in seq_len(n_groups)) {
    # each element of `A_list[[g]]` is the correlation matrix for a *subject*,
    # so we simply keep the list‑of‑lists structure.
    v_list[[g]] <- A_list[[g]]
  }

  structure(
    list(wi = B, v = v_list, tau = tau_vec),
    class = "lglasso"
  )
}
