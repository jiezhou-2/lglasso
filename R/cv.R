cv.lglasso=function (data = NULL,  lam = 10^seq(-2, 2, 0.2), type="general", group=NULL, expFix=1, diagonal = TRUE,
             tol = 1e-04, maxit = 10000, K = 5, start = c("warm","cold"), trace = c("progress", "print", "none"),
             ...)
{
  X=data[,-c(1,2)]
  crit.cv = match.arg(crit.cv)
  start = match.arg(start)
  trace = match.arg(trace)
  lam = sort(lam)
  Path = NULL
 # initmaxit = maxit
  CV_errors = array(0, c(length(lam), K))
  if (trace == "progress") {
    progress = txtProgressBar(max = K * length(lam), style = 3)
  }
  if (K == 1) {
      Path = array(0, c(ncol(X), ncol(X), length(lam)))
      S.train=cov(X)
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
   # maxit = initmaxit
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
      GLASSO = lglasso(data = data, rho = lam_, thr = tol,
                      maxit = maxit, penalize.diagonal = diagonal,
                      start = "warm", w.init = init, wi.init = initOmega,
                      trace = FALSE, ...)
      if (start == "warm") {
        init = GLASSO$w
        initOmega = GLASSO$wi
      }
      CV_errors[i, k] = (nrow(X)/2) * (sum(GLASSO$wi *
                                             S.valid) - determinant(GLASSO$wi, logarithm = TRUE)$modulus[1])

      if (K==1) {
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


