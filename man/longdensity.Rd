\name{Complete likelihood function for ith subject}
\alias{logdensity}
\title{Complete likelihood function}
\description{
 Complete likelihood function used in EM algorithm of heterogeneous marginal graphical lasso model
}
\usage{
logdensity(idata,omega,tau,alpha,ty)
}
\arguments{
  \item{idata}{ Data matrix for the subject i in which the first column is id for subject, the second column is
      the time points of observation.  Columns 2 to (p+2) is the observations for p variables.}
 \item{omega}{Precision matrix.}
  \item{tau}{Correlation parameter.}
  \item{alpha}{Parameter in exponential distribution}
  \item{ty}{Type of correlation function.}
}
\value{
 Value of complete likelihood function at given value of omega, tau and alpha
}


\author{Jie Zhou}

