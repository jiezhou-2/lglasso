
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
save(tau1,file="tau1.rd")
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
      z=sapply(x, logdensity,idata=idata,omega=omega,alpha=alpha0,type=type)
      tau1[i]=x[min(which(z==min(z)))]
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
  if (det(omega)<=10^(-20)) {
    warning("IN ilogdensity, omega is not poitive definite matrix!")
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



