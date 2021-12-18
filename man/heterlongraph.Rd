\name{Heterlongraph function}
\alias{heterlongraph}
\title{Estimating the heterogeneous SGGM }
\description{
 This function is used to estimate the heterogeneous stationary Gaussian graphical model.
}
\usage{
heterlongraph(data,rho,ty,tole,lower,upper)
}
\arguments{
\item{data}{A numeric vector representing the time points of observations.}
\item{rho}{Correlation parameter.}
\item{ty}{Type of correlation function.}
\item{tole}{tolerance error for the iteration}
\item{lower}{lower bound for correlation parameters}
\item{upper}{upper bound for correlation parameters}
}

\value{
A list consisting of estiamted network and correlation parameters. 
}

\author{Jie Zhou}

