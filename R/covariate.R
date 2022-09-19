
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
    stop("In icovalogdensity, omega is not poitive definite matrix!")
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
  return(a)
}






