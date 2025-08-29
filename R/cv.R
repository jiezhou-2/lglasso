cv.lglasso=function (data = NULL,  nlam=10,lam.min.ratio=0.01, lam = NULL, type="general", group=NULL, expFix=1, diagonal = TRUE,
             tol = 1e-04, maxit = 10000, K = 1, start = c("warm","cold"), trace = c("progress", "print", "none"),
             ...)
{

  if (is.null(data)) {
    stop("Must provide entry for data!")
  }
  if (!all(lam > 0)) {
    stop("lam must be positive!")
  }

  if (!(all(c(tol, maxit, K) > 0))) {
    stop("Entry must be positive!")
  }


  if (all(c(maxit, K)%%1 != 0)) {
    stop("Entry must be an integer!")
  }

  X=data[,-c(1,2)]
  S = (nrow(X) - 1)/nrow(X) * cov(X)
  crit.cv = match.arg(crit.cv)
  start = match.arg(start)
  trace = match.arg(trace)
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




    if (K > 1) {
      leave.out = ind[(1 + floor((k - 1) * n/K)):floor(k *
                                                         n/K)]
      X.train = X[-leave.out, , drop = FALSE]
      X_bar = apply(X.train, 2, mean)
      X.train = scale(X.train, center = X_bar, scale = FALSE)
      X.valid = X[leave.out, , drop = FALSE]
      X.valid = scale(X.valid, center = X_bar, scale = FALSE)
      S.train = crossprod(X.train)/(dim(X.train)[1])
      S.valid = crossprod(X.valid)/(dim(X.valid)[1])
    }

    # init = S.train
    # initOmega = diag(ncol(S.train))

lapply(lamda, lglasso)


    for (i in 1:length(lam)) {
      #lam_ = lam[i]

       # diag(init) = diag(S.train) + lam_

      LGLASSO = lglasso(data = data, lambda =   lam_, tol = tol,
                      maxit = maxit, penalize.diagonal = TRUE,
                      start = "warm", w.init = init, wi.init = initOmega,
                      trace = FALSE, ...)
      if (start == "warm") {
        init = LGLASSO$w
        initOmega = LGLASSO$wi
      }
      CV_errors[i, k] = (nrow(X)/2) * (sum(GLASSO$wi *
                                             S.valid) - determinant(GLASSO$wi, logarithm = TRUE)$modulus[1])

      if (K==1) {
        Path[, , i] = LGLASSO$wi
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

  AVG = apply(CV_errors, 1, mean)
  best_lam = lam[which.min(AVG)]
  error = min(AVG)
  return(list(lam = best_lam, path = Path, min.error = error,
              avg.error = AVG, cv.error = CV_errors))
}


