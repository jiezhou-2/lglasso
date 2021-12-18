\name{Log likelihood function}
\alias{ll_homo}
\title{Log likelihood function}
\description{
Value of likelihood function at given parameter. 
}
\usage{
ll_homo(data,omega,tau,ty)
}
\arguments{
  \item{data}{Data matrix  in which the first column is subject id, the second column is
the time points of observation.  Columns 2 to (p+2) is the observations for p variables.}
  \item{omega}{Precision matrix.}
  \item{tau}{tau  Correlation parameter}
  \item{ty}{Type of correlation function.}
}
\value{
 Value of likelihood function at given omega and tau  
}
\author{Jie Zhou}

