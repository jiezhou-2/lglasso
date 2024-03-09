#' Construct the temporal component fo correlation function
#'
#' @param t Time points of observations
#' @param tau correlation parameter
#' @param type The type of correlation function, which typically take either 0,1 or 2.
#' @author Jie Zhou
#' @export
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
#' @export
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
#' @export
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

  for (j in 1:(p-1)) {
    for (i in (j+1):p){
      precision[i,j]=precision[j,i]
    }
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

#' BIC function for i.i.d normal data
#'
#' @param data The data of interest
#' @param G The estimated network
#' @export
#' @return The value of bic
bicfunction=function(data,G,Tem){
  if (det(G)<=0){
    warning("G is not positive definite matrix!")
    return(bic=NA)
  }else{
    n=nrow(data)
    k=length(which(G!=0))/2
    yy=scale(as.matrix(data[,-c(1,2)]))
    ll=-0.5*sum(apply(yy,1,function(x)  t(x) %*% G %*% x ))+0.5*n*log(det(G))
    bic=-2*ll+k*log(n)+k*log(ncol(data)-2)/Tem
  }
  return(bic=bic)
}

#bbic=0.5*n*log(det(G))
#for (i in 1:nrow(yy)) {
# bbic=bbic-0.5* t(yy[i,])%*% G %*%yy[i,]
#}





#' full log likelihood used in EBIC computation
#' @param idata Data matrix for the subject i in which the first column is id for subject, the second column is
#' the time points of observation.  Columns 2 to (p+2) is the observations for p variables.
#' @param omega Precision matrix
#' @param tau Correlation parameter
#' @param type Type of correlation function, which can take either  "abs" or "qua".
#' @author Jie Zhou
#' @export
#' @return Value of likelihood function for subject i at given omega and tau

lli_homo=function(idata,omega,tau,type){

  if (det(omega)<=0) {
    # warning("IN ilogdensity, omega is not poitive definite matrix!")
    return(ll=-Inf)
  }

  t=idata[,2]
  yy=scale(as.matrix(idata[,-c(1,2)]),scale=F)
  p=nrow(omega)
  n=length(t)
  if (n==1){
    a1=(1/2)*log(det(omega))
    a2=-t(yy)%*%omega%*%yy/2
    a=a1+a2
  }else{
    a=0.5*log(det(omega))-t(yy[1,])%*%omega%*%yy[1,]/2
    for (i in 2:n) {
      cc=exp(-tau*(abs(t[i]-t[i-1]))^type)
      b=0.5*log(det(1/(1-cc^2)*omega))-t((yy[i,]-cc*yy[i-1,]))%*%(1/(1-cc^2)*omega)%*%(yy[i,]-cc*yy[i-1,])/2
      a=a+b
    }
  }
  return(a)
}



#' Value of likelihood function at given parameter
#'
#' @param data Data matrix  in which the first column is subject id, the second column is
#' the time points of observation.  Columns 2 to (p+2) is the observations for p variables.
#' @param omega Precision matrix
#' @param tau  Correlation parameter
#' @param type Type of correlation function, which can take either  "abs" or "qua".
#' @author Jie Zhou
#' @export
#' @return  Value of likelihood function at given omega and tau
ll_homo=function(data,omega,tau,type){
  id=unique(data[,1])
  a=0
  for (i in 1:length(id)) {
    idata=data[which(data[,1]==id[i]),]
    ai=lli_homo(idata = idata,omega = omega,tau = tau,type = type)
    a=a+ai
  }
  return(a)
}



#' Estiamte of precision matrix and autocorrelaton parameter for homogeneous model
#' @param data Data matrix  in which the first column is subject id, the second column is
#' the time points of observation.  Columns 2 to (p+2) is the observations for p variables.
#' @param type Type of correlation function, which can take either  "abs" or "qua".
#' @param rho Tuning parameter for graphical lasso
#' @param tole  Error tolerance for determination of convergence of EM algorithm
#' @param lower Lower bound for prediction of correlation parameter tau
#' @param upper Upper bound for prediction of correlation parameter tau
#' @author Jie Zhou
#' @return A list for estimates of precision matrix and correlation parameter for given tuning parameter
homolongraph=function(data, rho,type, tole,lower,upper){
  omega0=diag(ncol(data)-2)
  subject=data[,1]
  subjectid=unique(subject)
  m=length(subjectid)
  p=ncol(data)-2
  is=vector("list",m)
  x=seq(lower,upper,length=100)
  k=1
  while(1){
    z=sapply(x, ll_homo,data=data,omega=omega0,type=type)
    tau1=x[min(which(z==max(z)))]
    #print(paste("Iteration ",k,":","tau=",tau1))
    for (i in 1:m) {
      idata=data[which(data[,1]==subjectid[i]),]
      is[[i]]=iss(idata = idata,itau = tau1,type = type)
    }
    s=Reduce("+",is)/nrow(data)
    omega1=glasso::glasso(s=s,rho = rho)$wi
    if (max(abs(omega0-omega1))<tole){
      break
    }
    if (k>=15){
      warning("the algorithm does not converge")
      break
    }else{
      omega0=omega1
      k=k+1
    }
  }

  result=list(tau=tau1,omega=omega1,ll=max(z))
  return(result)
}



#' Estiamte of precision matrix and autocorrelaton parameter for homogeneous model
#' @param data Data matrix  in which the first column is subject id, the second column is
#' the time points of observation.  Columns 2 to (p+2) is the observations for p variables.
#' @param omega The maximum likelihood estiamte of precision matrix
#' @param type Type of correlation function, which can take either  "abs" or "qua".
#' @param lower Lower bound for prediction of correlation parameter tau
#' @param upper Upper bound for prediction of correlation parameter tau
#' @author Jie Zhou
#' @return A list for estimates of precision matrix and correlation parameter for given tuning parameter
mle_tau=function(data, omega, type,lower,upper){
  x=seq(lower,upper,length=50)
  z=sapply(x, ll_homo,data=data,omega=omega,type=type)
  tau1=x[min(which(z==max(z)))]
  return(tau1)
}












#' Estimates of correlation parameters and precision matrix
#'
#' @param data Data matrix  in which the first column is subject id, the second column is
#' the time points of observation.  Columns 2 to (p+2) is the observations for p variables.
#' @param type  Type of correlation function, which can take either  "abs" or "sqr".
#' @param rho Tuning parameter used in graphical lasso
#' @param tole Error tolerance for determination of convergence of EM algorithm
#' @param lower Lower bound for prediction of correlation parameter tau
#' @param upper Upper bound for prediction of correlation parameter tau
#' @author Jie Zhou
#' @importFrom glasso glasso
#' @return S list with three components which are the final estimate of alpha, tau and precision matrix omega
heterlongraph=function(data,rho, type,tole, lower,upper){
  omega0=diag(ncol(data)-2)
  subject=data[,1]
  subjectid=unique(subject)
  m=length(subjectid)
  alpha0=rep(0.01,m)
  p=ncol(data)-2
  tau_em=matrix(1/alpha0,nrow = 1,ncol = m)
  is=vector("list",m)
  tau1=rep(1,m)
  iidata=vector("list",m)
  k=1
  while(1){
    #  print(paste("Iteration ",k,":", " Mean value of tau",":", mean(tau1),sep=""))
    ll=0
    for (i in 1:m) {
      idata=data[which(data[,1]==subjectid[i]),]
      iidata[[i]]=idata
      if (nrow(idata)==1){
        tau1[i]=NA
        is[[i]]=iss(idata = idata,itau = 1,type = type)
        z=0.5*log(det(omega0))-0.5*t(idata)%*%omega0%*%idata
      }else{
        x=seq(lower,upper,length=50)
        z=sapply(x, ilogdensity,idata=idata,omega=omega0,alpha=alpha0[i],type=type)
        tau1[i]=x[min(which(z==max(z)))]
        if (is.na(tau1[i])){
          print(z)
          stop(paste0("tau1[[",i, "]] is missing!"))
        }
        is[[i]]=iss(idata = idata,itau = tau1[i],type = type)
      }
      ll=ll+max(z)
    }
    tau_em=rbind(tau_em,tau1)
    s=Reduce("+",is)/nrow(data)

    #save(is,file="is.rd")
    #save(s,file="netwrok.rd")
    #save(rho,file="rho.rd")
    omega1=glasso(s=s,rho = rho)$wi
    #G=glasso(s=s,rho = rho)$wi
    #omega1=mle_net(data=data,priori=G)
    alpha1=1/mean(tau1,na.rm=T)
    if (max(abs(tau_em[nrow(tau_em),]-tau_em[nrow(tau_em)-1,]),na.rm = T)<tole){
      break
    }
    if (k>=20){
      warning("the algorithm does not converge")
      break
    }else{
      alpha0=rep(1/mean(tau1,na.rm=T),m)
      tau0=tau1
      omega0=omega1
      k=k+1
    }

  }

  result=list(alpha=alpha1,tau=tau1,omega=omega1,ll=ll)
  return(result)
}






#' Maximum likelihood estimate of correlation parameter for given structure of precision matrix
#' @param data Data matrix  in which the first column is subject id, the second column is
#' the time points of observation.  Columns 2 to (p+2) is the observations for p variables.
#' @param alpha0 Initial value for the parameter in exponential distribution
#' @param omega Fixed value for precision matrix
#' @param type Type of correlation function, which can take either  "abs" or "qua".
#' @param tole  Error tolerance for determination of convergence of EM algorithm
#' @param lower Lower bound for prediction of correlation parameter tau
#' @param upper Upper bound for prediction of correlation parameter tau
#' @author Jie Zhou
#' @importFrom glasso glasso

mle_alpha=function(data,alpha0,omega, type, tole, lower,upper){
  subject=data[,1]
  subjectid=unique(subject)
  m=length(subjectid)
  p=ncol(data)-2
  tau_em=matrix(1/alpha0,nrow = 1,ncol = m)
  is=vector("list",m)
  tau1=rep(0,m)

  while(1){
    for (i in 1:m) {
      idata=data[which(data[,1]==subjectid[i]),]
      x=seq(lower,upper,length=50)
      z=sapply(x, ilogdensity,idata=idata,omega=omega,alpha=alpha0,type=type)
      tau1[i]=x[min(which(z==max(z)))]
      is[[i]]=iss(idata = idata,itau = tau1[i],type = type)
    }

    tau_em=rbind(tau_em,tau1)
    s=Reduce("+",is)/nrow(data)
    alpha1=1/mean(tau1)
    if (max(abs(tau_em[nrow(tau_em),]-tau_em[nrow(tau_em)-1,]))<tole){
      break
    }else{
      alpha0=alpha1
      tau0=tau1
    }

  }

  result=list(alpha=alpha1,tau=tau1)
  return(result)
}



#' Negative likelihood function used in EM algorithm of heterogeneous marginal graphical lasso model
#'
#' @param idata Data matrix for the subject i in which the first column is id for subject, the second column is
#' the time points of observation.  Columns 2 to (p+2) is the observations for p variables.
#' @param omega Precision matrix
#' @param tau Correlation parameter
#' @param alpha Parameter in exponential distribution
#' @param type Type of correlation function, which can take either  "abs" or "qua".
#' @author Jie Zhou
#' @export
#' @return Value of complete likelihood function at given value of omega, tau and alpha
ilogdensity=function(idata,omega,tau,alpha,type){
  if (det(omega)<=0) {
    # warning("IN ilogdensity, omega is not poitive definite matrix!")
    return(ll=-Inf)
  }
  t=idata[,2]
  yy=scale(as.matrix(idata[,-c(1,2)]),scale = F)
  p=nrow(omega)
  n=length(t)
  if (n==1){
    a=0.5*log(det(omega))-0.5*t(yy)%*%omega%*%yy
  }else{
    a=0.5*log(det(omega))-t(yy[1,])%*%omega%*%yy[1,]/2
    for (i in 2:n) {
      coe=exp(-tau*(abs(t[i]-t[i-1]))^type)
      b=0.5*log(det(1/(1-coe^2)*omega))-t((yy[i,]-coe*yy[i-1,]))%*%(1/(1-coe^2)*omega)%*%(yy[i,]-coe*yy[i-1,])/2
      a=a+b
    }
    ll=a+log(alpha)-alpha*tau
  }
  return(ll)
}



#' Estimates of correlation parameters and precision matrix
#'
#' @param data Data matrix  in which the first column is subject id, the second column is
#' the time points of observation.  Columns 2 to (p+2) is the observations for p variables.
#' @param x Data frame representing the covariates affecting the auto correlations. The frist column is the ids
#' @param type  Type of correlation function, which can take either  "abs" or "sqr".
#' @param rho Tuning parameter used in graphical lasso
#' @param tole Error tolerance for determination of convergence of EM algorithm
#' @param lower Lower bound for prediction of correlation parameter tau
#' @param upper Upper bound for prediction of correlation parameter tau
#' @author Jie Zhou
#' @importFrom glasso glasso
#' @return S list with three components which are the final estimate of alpha, tau and precision matrix omega
covaheterlongraph=function(data,covariates,rho, type,tole, lower,upper){
  coe0=rep(1,ncol(covariates))
  omega0=diag(ncol(data)-2)
  subject=data[,1]
  subjectid=unique(subject)
  m=length(subjectid)
  if (nrow(covariates)!=m){
    stop("the number of subjects in data does not match that in x!")
  }
  coid=unique(covariates[,1])
  if (!setequal(subjectid,coid)){
    stop("The id in data are not same as the id in covariates!")
  }
  p=ncol(data)-2
  tau_em=matrix(1,nrow = 1,ncol = m)
  is=vector("list",m)
  tau1=rep(0,m)
  k=1
  while(1){
    ll=0
    # print(paste("Iteration ",k))
    for (i in 1:m) {
      idata=data[which(data[,1]==subjectid[i]),]
      icovariates=c(1,covariates[which(covariates[,1]==subjectid[i]),-1])
      if (nrow(idata)==1){
        tau1[i]=NA
        is[[i]]=iss(idata = idata,itau = 1,type = type)
        z=icovalogdensity(idata = idata,icovariates = icovariates,omega = omega0,tau = tau1,coe = coe0,type = type)
      }else{
        x=seq(lower,upper,length=50)
        z=sapply(x, icovalogdensity,idata=idata,icovariates=icovariates,omega=omega0,coe=coe0,type=type)
        tau1[i]=x[min(which(z==max(z)))]
        is[[i]]=iss(idata = idata,itau = tau1[i],type = type)
      }
      ll=ll+max(z)
    }

    tau_em=rbind(tau_em,tau1)
    s=Reduce("+",is)/nrow(data)
    omega1=glasso::glasso(s=s,rho = rho)$wi
    #G=glasso::glasso(s=s,rho = rho)$wi
    #omega1=mle_net(data=data,priori=G)
    alikehood=function(coe){
      index=which(!is.na(tau1))
      alpha=covariates[,-1]%*%coe[-1]+coe[1]
      a1=sum(alpha)
      a2=t(exp(alpha)[index])%*%tau1[index]
      c=-a1+a2
    }


    coe1=optim(rep(1,ncol(covariates)),alikehood)$par
    if (max(abs(tau_em[nrow(tau_em),]-tau_em[nrow(tau_em)-1,]),na.rm = T)<tole){
      break
    }
    if (k>=20){
      warning("the algorithm does not converge")
      break
    }else{
      coe0=coe1
      tau0=tau1
      omega0=omega1
      k=k+1
    }

  }

  result=list(coe=coe1,tau=tau1,omega=omega1,ll=ll)
  return(result)
}






#' Maximum likelihood estimate of correlation parameter for given structure of precision matrix
#' @param data Data matrix  in which the first column is subject id, the second column is
#' the time points of observation.  Columns 2 to (p+2) is the observations for p variables.
#' @param alpha0 Initial value for the parameter in exponential distribution
#' @param omega Fixed value for precision matrix
#' @param type Type of correlation function, which can take either  "abs" or "qua".
#' @param tole  Error tolerance for determination of convergence of EM algorithm
#' @param lower Lower bound for prediction of correlation parameter tau
#' @param upper Upper bound for prediction of correlation parameter tau
#' @author Jie Zhou
#' @importFrom glasso glasso

covamle_alpha=function(data,covariates,alpha0,omega, type, tole, lower,upper){
  coe0=rep(1,ncol(covariates))
  subject=data[,1]
  subjectid=unique(subject)
  m=length(subjectid)
  p=ncol(data)-2
  tau_em=matrix(1/alpha0,nrow = 1,ncol = m)
  is=vector("list",m)
  tau1=rep(0,m)
  coid=unique(covariates[,1])
  if (!setequal(subjectid,coid)){
    stop("The id in data are not same as the id in covariates!")
  }
  while(1){
    for (i in 1:m) {
      idata=data[which(data[,1]==subjectid[i]),]
      icovariates=covariates[which(covariates[,1]==subjectid[i]),-1]
      x=seq(lower,upper,length=50)
      z=sapply(x, covalogdensity,idata=idata,covariates=icovariates,omega=omega,coe=coe0,type=type)
      tau1[i]=x[min(which(z==min(z)))]
      is[[i]]=iss(idata = idata,itau = tau1[i],type = type)
    }

    tau_em=rbind(tau_em,tau1)

    alikehood=function(coe){
      index=which(!is.na(tau1))
      alpha=covariates[,-1]%*%coe[-1]+coe[1]
      a1=sum(alpha)
      a2=t(exp(alpha)[index])%*%tau1[index]
      c=-a1+a2
    }


    coe1=optim(rep(1,ncol(covariates)),alikehood)$par
    if (max(abs(tau_em[nrow(tau_em),]-tau_em[nrow(tau_em)-1,]))<tole){
      break
    }else{
      coe0=coe1
      tau0=tau1
    }

  }

  result=list(coe=coe1,tau=tau1)
  return(result)
}


#' Negative likelihood function used in EM algorithm of heterogeneous marginal graphical lasso model
#'
#' @param idata Data matrix for the subject i in which the first column is id for subject, the second column is
#' the time points of observation.  Columns 2 to (p+2) is the observations for p variables.
#' @param x Vector standing for the covariates affecting auto correlations
#' @param omega Precision matrix
#' @param tau Correlation parameter
#' @param alpha Parameter in exponential distribution
#' @param type Type of correlation function, which can take either  "abs" or "qua".
#' @author Jie Zhou
#' @return Value of complete likelihood function at given value of omega, tau and alpha
icovalogdensity=function(idata,icovariates,omega,tau,coe,type){
  if (det(omega)<=10^(-20)) {
    warning("In icovalogdensity, omega is not poitive definite matrix!")
    return(ll=-Inf)
  }
  if (length(icovariates)!= length(coe)){
    stop(("The length of icovariate should be equal to coeffiient coe!"))
  }
  t=idata[,2]
  yy=scale(as.matrix(idata[,-c(1,2)]),scale = F)
  p=nrow(omega)
  n=length(t)

  if (n==1){
    a=0.5*log(det(omega))-0.5*t(yy)%*%omega%*%yy
  }else{
    a=0.5*log(det(omega))-t(yy[1,])%*%omega%*%yy[1,]/2
    for (i in 2:n) {
      cc=exp(-tau*(abs(t[i]-t[i-1]))^type)
      b=0.5*log(det(1/(1-cc^2)*omega))-t((yy[i,]-cc*yy[i-1,]))%*%(1/(1-cc^2)*omega)%*%(yy[i,]-cc*yy[i-1,])/2
      a=a+b
    }
    a3=exp(t(coe[-1])%*%icovariates[-1]+coe[1])
    a=a+log(a3)-a3*tau
  }
  return(ll=a)
}








