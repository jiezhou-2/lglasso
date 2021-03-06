\name{Log likelihood function for ith subject.}
\alias{lli_homo}
\title{Log likelihood function for ith Subject}
\description{
Value of likelihood function for subject i at given omega and tau.
}
\usage{
lli_homo(idata,omega,tau,ty)
}
\arguments{
  \item{idata}{ Data matrix for the subject i in which the first column is id for subject, the second column is
 the time points of observation.  Columns 2 to (p+2) is the observations for p variables.}
 \item{omega}{Precision matrix.}
  \item{tau}{Correlation parameter.}
  \item{ty}{Type of correlation function.}
}
\value{
Value of likelihood function for subject i at given omega and tau
  
}

\author{Jie Zhou}

