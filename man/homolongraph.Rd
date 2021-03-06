\name{Homolongraph function}
\alias{homolongraph}
\title{Estimating the homogeneous SGGM }
\description{
 This function is used to estimate the homogeneous stationary Gaussian graphical model.
}
\usage{
homolongraph(data, rho,ty, tole,lower,upper)
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

