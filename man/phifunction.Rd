\name{Correlation function}
\alias{phifunction}
\title{Construct the Correlation Matrix}
\description{
 This function is used to construct the likelihood function  for stationary Gaussian graphical model.
}
\usage{
phifunction(t,tau,ty)
}
\arguments{
  \item{t}{A numeric vector representing the time points of observations.}
  \item{tau}{Correlation parameter.}
  \item{ty}{Type of correlation function.}
}
\value{
A square matrix \code{M} with dimension equal to the lengthe of \code{t} 
  
}

\author{Jie Zhou}

