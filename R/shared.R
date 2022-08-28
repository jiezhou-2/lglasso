#' Construct the temporal component fo correlation function
#'
#' @param t Time points of observations
#' @param tau correlation parameter
#' @param type The type of correlation function, which typically take either 0,1 or 2.
#' @author Jie Zhou
#' @return A square matrix with dimension equal to the length of vector t
phifunction=function(t,tau,type=1){
  n=length(t)
  if (n==1){
    M=as.matrix(1)
  }else{
    M=matrix(nrow = n, ncol = n)
    for (i in 1:n) {
      for (j in i:n){
        M[i,j]=exp(-tau*(abs(t[i]-t[j]))^type)
        M[j,i]=M[i,j]
      }
    }
    diag(M)=1
  }
  return(M)
}


#' Quasi covariance matrix for subject i
#' @param idata Data matrix for the subject i in which the first column is subject (cluster) id, the second column stands for
#' the time points () of observation.  Columns 2 to (p+2) is the observations for p variables respectively.
#' @param itau Correlation parameter
#' @param type  Type of correlation function, which typically take either  0, 1 or 2.
#' @author Jie Zhou
#' @return Empirical quasi covariance matrix
iss=function(idata,itau,type){
  t=idata[,2]
  p=ncol(idata)-2
  inversephi= solve(phifunction(t=t,tau = itau,type = type))
  si=matrix(0,nrow = p,ncol = p)
  yy=as.matrix(idata[,-c(1,2)])
  yy=scale(yy, scale = F)
  for (j in 1:length(t)) {
    for (k in 1:length(t)) {
      si=si+inversephi[j,k]*(yy[j,])%*%t(yy[k,])

    }
  }
  return(si)
}



##For a given network matrix, compute
##MLE of precision matrix
#' Title
#'
#' @param data A Longitudinal data set
#'
#'
#' @param priori Given structure of precision matrix
#' @author Jie Zhou
#' @return The maximum likelihood estimation

mle_net=function(data,priori){
  data=data[,-c(1,2)]
  priori=priori+t(priori)
  priori=ifelse(priori==0,0,1)
  diag(priori)=0
  p=dim(data)[2]
  n=dim(data)[1]
  precision=matrix(0,nrow = p,ncol = p)
  ##for the first row
  for (j in 1:p) {
    if (sum(priori[j,])>=nrow(data)) {
      stop("The number of unknown parameters exceeds the sample size!")
    }
    data=scale(data,scale = F)
    y=data[,j]
    index=which(priori[j,]==1)
    if (length(index)==0) {
      precision[j,]=0
      precision[j,j]=1/stats::var(y)
    }else{
      x=data[,index]
      result=stats::lm(y~0+x)
      alpha=result$coefficients
      sigma=t(result$residuals)%*%(result$residuals)/result$df.residual
      precision[j,index]=-alpha/sigma[1]
      precision[j,j]=1/sigma
    }
    precision=(precision+t(precision))/2
  }
  return(precision)
}




#' @title Graphical Lasso for Longitudinal Data
#' @description This function implements the L_1 penalized maximum likelihood estimation for precision matrix (network)  based on correlated data, e.g., irregularly spaced longitudinal
#'  data. It can be regarded as an extension of the package \code{glasso} (Friedman,Hastie and Tibshirani, 2008) which aims
#'  to find the sparse estimate of the network from independent continuous data.
#' @param data Data matrix  in which the first column is subject id, the second column is
#'  time points of observations for temporal data or site id for spatial data.  Columns \code{3} to \code{(p+2)} is the observations for \code{p} variables.
#' @param rho Tuning parameter used in \code{L_1} penalty
#' @param heter Binary variable \code{TRUE} or \code{FALSE}, indicating heterogeneous model or homogeneous model is fitted. In heterogeneous model,
#' subjects are allowed to have his/her own temporal correlation parameter \code{tau_i}; while in homogeneous model, all the subjects are assumed to
#'  share the same temporal correlation parameter,i.e., \code{tau_1=tau_2=...tau_m}.
#' @param type A positive number which specify the correlation function. The general form of correlation function  is given by \code{ exp(tau|t_i-t_j|^type)}.
#' in which \code{type=0} can be used for spatial correlation while \code{type>0} are used for temporal correlation. For latter, the default value is set to be \code{type=1}.
#' @param tole Threshold for convergence. Default value is \code{1e-2}. Iterations stop when maximum
#' absolute difference between consecutive estimates of parameter change is less than \code{tole}.
#' @param lower  Lower bound for predicts of correlation parameter \code{tau}.
#' Default value is \code{1e-2}. The estimate of \code{tau}(\code{alpha}) will be searched in the
#' interval \code{[lower,upper]}, where parameter \code{upper} is explained in the following.
#' @param upper Upper bound for predicts of correlation parameter \code{tau}.
#' @author Jie Zhou
#' @references  Jie Zhou, Jiang Gui, Weston D.Viles, Anne G.Hoen Identifying Microbial Interaction Networks Based on Irregularly Spaced Longitudinal 16S rRNA sequence data. bioRxiv 2021.11.26.470159; doi: https://doi.org/10.1101/2021.11.26.470159
#' @references Friedman J, Tibshirani TH and R. Glasso: Graphical Lasso: Estimation of Gaussian Graphical Models.; 2019. Accessed November 28, 2021. https://CRAN.R-project.org/package=glasso
#' @references Friedman J, Hastie T, Tibshirani TH, Sparse inverse covariance estimation with the graphical lasso, Biostatistics, Volume 9, Issue 3, July 2008, Pages 432â441, https://doi.org/10.1093/biostatistics/kxm045
#' @return  If \code{heter=TRUE}, then a list with three components is returned which are  respectively
#' the estimate of parameter \code{alpha} in exponent distribution, correlation parameter \code{tau} and precision matrix \code{omega}. If \code{heter=FALSE},
#' then a list with two components is returned which are respectively the estimate of correlation parameter \code{tau} and precision matrix \code{omega}.
#' @export

lglasso=function(data,x=NULL, rho,heter=TRUE,type=1, tole=0.01,lower=0.01,upper=10){
  if (heter==TRUE){
    if (is.null(x)){
    aa=heterlongraph(data=data,rho=rho,type=type,tole=tole,lower=lower,upper=upper)
  }else{
    aa=covaheterlongraph(data=data,covariates=x,rho=rho,type=type,tole=tole,lower=lower,upper=upper)
  }
    }else{
      aa=homolongraph(data=data,rho=rho,type=type,tole=tole,lower=lower,upper=upper)
  }

  return(aa)

}




#' Maximum Likelihood Estimate of Precision Matrix and Correlation Parameters for Given Network
#' @param data Data matrix  in which the first column is subject id, the second column is
#'  time points of observations for temporal data or site id for spatial data.
#'   Columns \code{3} to \code{(p+2)} is the observations for \code{p} variables.
#' @param network The network selected by function lglasso
#' @param heter Binary variable \code{TRUE} or \code{FALSE}, indicating heterogeneous model or homogeneous model is fitted. In heterogeneous model,
#' subjects are allowed to have his/her own temporal correlation parameter \code{tau_i}; while in homogeneous model, all the subjects are assumed to
#'  share the same temporal correlation parameter,i.e., \code{tau_1=tau_2=...tau_m}.
#' @param type  A positive number which specify the correlation function. The general form of correlation function  is given by \code{ exp(tau|t_i-t_j|^type)}.
#' in which \code{type=0} can be used for spatial correlation while \code{type>0} are used for temporal correlation. For latter, the default value is set to be \code{type=1}.
#' @param tole Threshold for convergence. Default value is \code{1e-2}. Iterations stop when maximum
#' absolute difference between consecutive estimates of parameter change is less than \code{tole}.
#' @param lower Lower bound for predicts of correlation parameter \code{tau}.
#' Default value is \code{1e-2}. The estimate of \code{tau}(\code{alpha}) will be searched in the
#' interval \code{[lower,upper]}, where parameter \code{upper} is explained in the following.
#' @param upper Upper bound for predicts of correlation parameter \code{tau}.
#'
#' @return A list which include the maximum likelihood estimate of precision matrix, correlation parameter \code{tau}. If \code{heter=TRUE},
#' the output also include the estimate of alpha where \code{tau~exp(alpha)}
#' @author Jie Zhou
#' @export
mle=function(data,x=NULL, network,heter=TRUE,type=1,tole=0.01,lower=0.01,upper=10){
  mlenetwork=mle_net(data = data,priori = network)
  if (heter==TRUE){
    if (is.null(x)){
    mle=mle_alpha(data=data,alpha0=1,omega=mlenetwork, type=type, tole=tole, lower=lower,upper=upper)
    }else{
    mle=covamle_alpha(data=data,covariates=x,alpha0=1,omega=mlenetwork, type=type, tole=tole, lower=lower,upper=upper)
    }
    }else{
    mle=mle_tau(data=data, omega=mlenetwork, type=type,lower=lower,upper=upper)
  }
  if (heter==TRUE){
    result=list(network=mlenetwork,alpha=mle$alpha,tau=mle$tau)
  }
  if (heter==FALSE){
    result=list(network=mlenetwork,tau=mle)
  }
  return(result)
}

