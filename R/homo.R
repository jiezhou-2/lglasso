


#' full log likelihood used in EBIC computation
#' @param idata Data matrix for the subject i in which the first column is id for subject, the second column is
#' the time points of observation.  Columns 2 to (p+2) is the observations for p variables.
#' @param omega Precision matrix
#' @param tau Correlation parameter
#' @param type Type of correlation function, which can take either  "abs" or "qua".
#' @author Jie Zhou
#' @return Value of likelihood function for subject i at given omega and tau

lli_homo=function(idata,omega,tau,type){
  t=idata[,2]
  yy=as.matrix(idata[,-c(1,2)])
  y=c(t(scale(yy, scale = F)))
  p=nrow(omega)
  n=length(t)
  if (n*p!=length(y)){
    stop("the dimensions do not match")
  }

  phi=phifunction(t=t,tau = tau,type = type)
  a0=(n/2)*log(det(omega))
  a1=(-p/2)*log(det(phi))
  a2=-t(y)%*%kronecker(solve(phi),omega)%*%y/2
  a=a0+a1+a2
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
#' @return  Value of likelihood function at given omega and tau
ll_homo=function(data,omega,tau,type){
  id=unique(data[,1])
  a=0
  for (i in id) {
    idata=data[which(data[,1]==i),]
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
  tau0=1
  omega0=diag(ncol(data)-2)
  subject=data[,1]
  subjectid=unique(subject)
  m=length(subjectid)
  p=ncol(data)-2
  tau=tau0
  is=vector("list",m)
  k=1
  #tau1=rep(0,m)
  while(1){
    x=seq(lower,upper,length=50)
    z=sapply(x, ll_homo,data=data,omega=omega0,type=type)
    tau1=x[min(which(z==max(z)))]
    print(paste("Iteration ",k,":","tau=",tau1))
    for (i in 1:m) {
      idata=data[which(data[,1]==subjectid[i]),]
      is[[i]]=iss(idata = idata,itau = tau1,type = type)
    }
    s=Reduce("+",is)/nrow(data)
    omega1=glasso::glasso(s=s,rho = rho)$wi
    if (abs(tau0-tau1)<tole & max(abs(omega0-omega1))<tole){
      break
    }else{
      tau0=tau1
      omega0=omega1
      k=k+1
    }
  }

  result=list(tau=tau1,omega=omega1)
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











