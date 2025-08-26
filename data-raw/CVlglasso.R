CV=function (X = NULL, S = NULL, lam = 10^seq(-2, 2, 0.2), diagonal = FALSE,
          path = FALSE, tol = 1e-04, maxit = 10000, adjmaxit = NULL,
          K = 5, crit.cv = c("loglik", "AIC", "BIC"), start = c("warm",
                                                                "cold"),
          cores = 1, trace = c("progress", "print", "none"),
          ...)
{
  crit.cv = match.arg(crit.cv)
  start = match.arg(start)
  trace = match.arg(trace)
  lam = sort(lam)
  Path = NULL
  initmaxit = maxit
  S.train = S.valid = S
  CV_errors = array(0, c(length(lam), K))
  if (trace == "progress") {
    progress = txtProgressBar(max = K * length(lam), style = 3)
  }
  if (K == 1) {
    n = nrow(S)
    if (path) {
      Path = array(0, c(ncol(S), ncol(S), length(lam)))
    }
  }
  else {
    n = nrow(X)
    ind = sample(n)
  }
  for (k in 1:K) {
    if (K > 1) {
      leave.out = ind[(1 + floor((k - 1) * n/K)):floor(k *
                                                         n/K)]
      X.train = X[-leave.out, , drop = FALSE]
      X_bar = apply(X.train, 2, mean)
      X.train = scale(X.train, center = X_bar, scale = FALSE)
      X.valid = X[leave.out, , drop = FALSE]
      X.valid = scale(X.valid, center = X_bar, scale = FALSE)
      n = nrow(X.valid)
      S.train = crossprod(X.train)/(dim(X.train)[1])
      S.valid = crossprod(X.valid)/(dim(X.valid)[1])
    }
    maxit = initmaxit
    init = S.train
    initOmega = diag(ncol(S.train))
    if (!diagonal) {
      Sminus = S.train
      diag(Sminus) = 0
      alpha = min(c(lam[1]/max(abs(Sminus)), 1))
      init = (1 - alpha) * S.train
      diag(init) = diag(S.train)
    }
    for (i in 1:length(lam)) {
      lam_ = lam[i]
      if (diagonal) {
        diag(init) = diag(S.train) + lam_
      }
      GLASSO = glasso(s = S.train, rho = lam_, thr = tol,
                      maxit = maxit, penalize.diagonal = diagonal,
                      start = "warm", w.init = init, wi.init = initOmega,
                      trace = FALSE, ...)
      if (start == "warm") {
        init = GLASSO$w
        initOmega = GLASSO$wi
        maxit = adjmaxit
      }
      CV_errors[i, k] = (nrow(X)/2) * (sum(GLASSO$wi *
                                             S.valid) - determinant(GLASSO$wi, logarithm = TRUE)$modulus[1])
      if (crit.cv == "AIC") {
        CV_errors[i, k] = CV_errors[i, k] + sum(GLASSO$wi !=
                                                  0)
      }
      if (crit.cv == "BIC") {
        CV_errors[i, k] = CV_errors[i, k] + sum(GLASSO$wi !=
                                                  0) * log(nrow(X))/2
      }
      if (path) {
        Path[, , i] = GLASSO$wi
      }
      if (trace == "progress") {
        setTxtProgressBar(progress, i + (k - 1) * length(lam))
      }
      else if (trace == "print") {
        cat("\nFinished lam = ", paste(lam_, sep = ""))
      }
    }
    if (trace == "print") {
      cat("\nFinished fold ", paste(k, sep = ""))
    }
  }
  AVG = apply(CV_errors, 1, mean)
  best_lam = lam[which.min(AVG)]
  error = min(AVG)
  return(list(lam = best_lam, path = Path, min.error = error,
              avg.error = AVG, cv.error = CV_errors))
}




CVglasso=function (X = NULL, S = NULL, nlam = 10, lam.min.ratio = 0.01,
                   lam = NULL, diagonal = FALSE, path = FALSE, tol = 1e-04,
                   maxit = 10000, adjmaxit = NULL, K = 5, crit.cv = c("loglik",
                                                                      "AIC", "BIC"),
                   start = c("warm", "cold"), cores = 1,
                   trace = c("progress", "print", "none"), ...)
{
  if (is.null(X) && is.null(S)) {
    stop("Must provide entry for X or S!")
  }
  if (!all(lam > 0)) {
    stop("lam must be positive!")
  }
  if (!(all(c(tol, maxit, adjmaxit, K, cores) > 0))) {
    stop("Entry must be positive!")
  }
  if (!(all(sapply(c(tol, maxit, adjmaxit, K, cores, nlam,
                     lam.min.ratio), length) <= 1))) {
    stop("Entry must be single value!")
  }
  if (all(c(maxit, adjmaxit, K, cores)%%1 != 0)) {
    stop("Entry must be an integer!")
  }
  if (cores < 1) {
    stop("Number of cores must be positive!")
  }
  if (cores > 1 && path) {
    cat("\nParallelization not possible when producing solution path. Setting cores = 1...")
    cores = 1
  }
  K = ifelse(path, 1, K)
  if (cores > K) {
    cat("\nNumber of cores exceeds K... setting cores = K")
    cores = K
  }
  if (is.null(adjmaxit)) {
    adjmaxit = maxit
  }
  crit.cv = match.arg(crit.cv)
  start = match.arg(start)
  trace = match.arg(trace)
  call = match.call()
  MIN.error = AVG.error = CV.error = NULL
  n = ifelse(is.null(X), nrow(S), nrow(X))
  if (is.null(S)) {
    S = (nrow(X) - 1)/nrow(X) * cov(X)
  }
  Sminus = S
  diag(Sminus) = 0
  if (is.null(lam)) {
    if (!((lam.min.ratio <= 1) && (lam.min.ratio > 0))) {
      cat("\nlam.min.ratio must be in (0, 1]... setting to 1e-2!")
      lam.min.ratio = 0.01
    }
    if (!((nlam > 0) && (nlam%%1 == 0))) {
      cat("\nnlam must be a positive integer... setting to 10!")
      nlam = 10
    }
    lam.max = max(abs(Sminus))
    lam.min = lam.min.ratio * lam.max
    lam = 10^seq(log10(lam.min), log10(lam.max), length = nlam)
  }
  else {
    lam = sort(lam)
  }
  if ((length(lam) > 1) & (!is.null(X) || path)) {
    if (cores > 1) {
      GLASSO = CVP(X = X, lam = lam, diagonal = diagonal,
                   tol = tol, maxit = maxit, adjmaxit = adjmaxit,
                   K = K, crit.cv = crit.cv, start = start, cores = cores,
                   trace = trace, ...)
      MIN.error = GLASSO$min.error
      AVG.error = GLASSO$avg.error
      CV.error = GLASSO$cv.error
    }
    else {
      if (is.null(X)) {
        X = matrix(0)
      }
      GLASSO = CV(X = X, S = S, lam = lam, diagonal = diagonal,
                  path = path, tol = tol, maxit = maxit, adjmaxit = adjmaxit,
                  K = K, crit.cv = crit.cv, start = start, trace = trace,
                  ...)
      MIN.error = GLASSO$min.error
      AVG.error = GLASSO$avg.error
      CV.error = GLASSO$cv.error
      Path = GLASSO$path
    }
    if ((GLASSO$lam == lam[1]) && (length(lam) != 1) && !path) {
      cat("\nOptimal tuning parameter on boundary... consider providing a smaller lam value or decreasing lam.min.ratio!")
    }
    if (diagonal) {
      init = S + GLASSO$lam
    }
    else {
      alpha = min(c(GLASSO$lam/max(abs(Sminus)), 1))
      init = (1 - alpha) * S
      diag(init) = diag(S)
    }
    lam_ = GLASSO$lam
    GLASSO = glasso(s = S, rho = lam_, thr = tol, maxit = maxit,
                    penalize.diagonal = diagonal, start = "warm", w.init = init,
                    wi.init = diag(ncol(S)), trace = FALSE, ...)
    GLASSO$lam = lam_
  }
  else {
    if (length(lam) > 1) {
      stop("Must set specify X, set path = TRUE, or provide single value for lam.")
    }
    if (diagonal) {
      init = S + lam
    }
    else {
      alpha = min(c(lam/max(abs(Sminus)), 1))
      init = (1 - alpha) * S
      diag(init) = diag(S)
    }
    GLASSO = glasso(s = S, rho = lam, thr = tol, maxit = maxit,
                    penalize.diagonal = diagonal, start = "warm", w.init = init,
                    wi.init = diag(ncol(S)), trace = FALSE, ...)
    GLASSO$lam = lam
  }
  if (diagonal) {
    C = 1
  }
  else {
    C = 1 - diag(ncol(S))
  }
  loglik = (-n/2) * (sum(GLASSO$wi * S) - determinant(GLASSO$wi,
                                                      logarithm = TRUE)$modulus[1] + GLASSO$lam * sum(abs(C *
                                                                                                            GLASSO$wi)))
  tuning = matrix(c(log10(GLASSO$lam), GLASSO$lam), ncol = 2)
  colnames(tuning) = c("log10(lam)", "lam")
  if (!path) {
    Path = NULL
  }
  returns = list(Call = call, Iterations = GLASSO$niter, Tuning = tuning,
                 Lambdas = lam, maxit = maxit, Omega = GLASSO$wi, Sigma = GLASSO$w,
                 Path = Path, Loglik = loglik, MIN.error = MIN.error,
                 AVG.error = AVG.error, CV.error = CV.error)
  class(returns) = "CVglasso"
  return(returns)
}



CVP=function (X = NULL, lam = 10^seq(-2, 2, 0.2), diagonal = FALSE,
              tol = 1e-04, maxit = 10000, adjmaxit = NULL, K = 5, crit.cv = c("loglik",
                                                                              "AIC", "BIC"),
              start = c("warm", "cold"), cores = 1,
              trace = c("progress", "print", "none"), ...)
{
  crit.cv = match.arg(crit.cv)
  start = match.arg(start)
  trace = match.arg(trace)
  lam = sort(lam)
  num_cores = detectCores()
  if (cores > num_cores) {
    cat("\nOnly detected", paste(num_cores, "cores...", sep = " "))
  }
  if (cores > K) {
    cat("\nNumber of cores exceeds K... setting cores = K")
    cores = K
  }
  cluster = makeCluster(cores)
  registerDoParallel(cluster)
  n = nrow(X)
  ind = sample(n)
  k = NULL
  CV = foreach(k = 1:K, .packages = "CVglasso", .combine = "cbind",
               .inorder = FALSE) %dopar% {
                 if (trace == "progress") {
                   progress = txtProgressBar(max = length(lam), style = 3)
                 }
                 leave.out = ind[(1 + floor((k - 1) * n/K)):floor(k *
                                                                    n/K)]
                 X.train = X[-leave.out, , drop = FALSE]
                 X_bar = apply(X.train, 2, mean)
                 X.train = scale(X.train, center = X_bar, scale = FALSE)
                 X.valid = X[leave.out, , drop = FALSE]
                 X.valid = scale(X.valid, center = X_bar, scale = FALSE)
                 S.train = crossprod(X.train)/(dim(X.train)[1])
                 S.valid = crossprod(X.valid)/(dim(X.valid)[1])
                 init = S.train
                 initOmega = diag(ncol(S.train))
                 CV_error = array(0, length(lam))
                 if (!diagonal) {
                   Sminus = S.train
                   diag(Sminus) = 0
                   alpha = min(c(lam[1]/max(abs(Sminus)), 1))
                   init = (1 - alpha) * S.train
                   diag(init) = diag(S.train)
                 }
                 for (i in 1:length(lam)) {
                   lam_ = lam[i]
                   if (diagonal) {
                     diag(init) = diag(S.train) + lam_
                   }
                   GLASSO = glasso(s = S.train, rho = lam_, thr = tol,
                                   maxit = maxit, penalize.diagonal = diagonal,
                                   start = "warm", w.init = init, wi.init = initOmega,
                                   trace = FALSE, ...)
                   if (start == "warm") {
                     init = GLASSO$w
                     initOmega = GLASSO$wi
                     maxit = adjmaxit
                   }
                   CV_error[i] = (nrow(X)/2) * (sum(GLASSO$wi * S.valid) -
                                                  determinant(GLASSO$wi, logarithm = TRUE)$modulus[1])
                   if (crit.cv == "AIC") {
                     CV_error[i] = CV_error[i] + sum(GLASSO$wi !=
                                                       0)
                   }
                   if (crit.cv == "BIC") {
                     CV_error[i] = CV_error[i] + sum(GLASSO$wi !=
                                                       0) * log(nrow(X))/2
                   }
                   if (trace == "progress") {
                     setTxtProgressBar(progress, i)
                   }
                   else if (trace == "print") {
                     cat("\nFinished lam = ", paste(lam_, sep = ""))
                   }
                 }
                 return(CV_error)
               }
  AVG = apply(CV, 1, mean)
  best_lam = lam[which.min(AVG)]
  error = min(AVG)
  stopCluster(cluster)
  return(list(lam = best_lam, min.error = error, avg.error = AVG,
              cv.error = CV))
}
