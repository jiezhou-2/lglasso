\name{Maximum Likelihood Estiamtion (mle)}
\alias{mle_net}
\alias{mle_alpha}
\alias{mle_tau}
\title{Set of functions for computing the mle.}
\description{
 These functions are used to compute the likelihood function  for stationary Gaussian graphical model.
}
\usage{
mle_net(data,priori)
mle_alpha(data,alpha0,omega, ty, tole, lower,upper)
mle_tau(data, omega, ty,lower,upper)
}
\arguments{
  \item{data}{Data matrix  in which the first column is subject id, the second column is
 the time points of observation.  Columns 2 to (p+2) is the observations for p variables.}
  \item{priori}{Given structure of precision matrix.}
  \item{alpha0}{Initial value for the parameter in exponential distribution}
  \item{omega}{Precision matrix}
  \item{tole}{Error tolerance for determination of convergence of EM algorithm}
  \item{lower}{Lower bound for prediction of correlation parameter tau}
  \item{upper}{Upper bound for prediction of correlation parameter tau}
  \item{ty}{Type of correlation function}
}
\value{
 Mle of precision matrix, parameter in exponential distribution and correlation parameters for mle_net, mle_alpha and. mle_tau respectively.
}
\author{Jie Zhou}

