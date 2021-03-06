\name{Maximum Likelihood Estimate of Precision Matrix and Correlation Parameters for Given Network}
\alias{mle}
\title{Maximum Likelihood Estimate of Precision Matrix and Correlation Parameters for Given Network}
\description{
Maximum Likelihood Estimate of Precision Matrix and Correlation Parameters for Given Network
}
\usage{
mle(data,network,heter,ty,tole,lower,upper)
}
\arguments{
  \item{data}{Data matrix  in which the first column is subject id, the second column is
 the time points of observation.  Columns 2 to (p+2) is the observations for p variables.}
  \item{heter}{Binary variable \code{TRUE} or \code{FALSE}, indicating heterogeneous model or homogeneous model is fitted. In heterogeneous model, subjects are allowed to have his/her own temporal correlation parameter \code{tau_i}; while in homogeneous model, all the subjects are assumed to  share the same temporal correlation parameter,i.e., \code{tau_1=tau_2=...tau_m}.}
  \item{network}{The structure of precision matrix}
  \item{tole}{Error tolerance for determination of convergence of EM algorithm}
  \item{lower}{Lower bound for prediction of correlation parameter tau}
  \item{upper}{Upper bound for prediction of correlation parameter tau}
  \item{ty}{Type of correlation function}
}
\value{
A list which include the maximum likelihood estimate of precision matrix, correlation parameter \code{tau}. If \code{heter=TRUE},
the output also include the estimate of alpha where \code{tau~exp(alpha)}
}

\author{Jie Zhou}

