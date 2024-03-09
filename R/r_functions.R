
#' Quasi covariance matrix for subject i
#' @param idata Data matrix for the subject i in which the first column is subject (cluster) id, the second column stands for
#' the time points () of observation.  Columns 2 to (p+2) is the observations for p variables respectively.
#' @param itau Correlation parameter
#' @param type  Type of correlation function, which typically take either  0, 1 or 2.
#' @author Jie Zhou
#' @export
#' @noRd
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


#' Construct the temporal component fo correlation function
#'
#' @param t Time points of observations
#' @param tau correlation parameter
#' @param type The type of correlation function, which typically take either 0,1 or 2.
#' @author Jie Zhou
#' @export
#' @noRd
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

#' Title
#'
#' @param data
#' @param rho
#' @param ty
#' @param tole
#' @param lower
#' @param upper
#' @export
#' @return
#' @noRd
heterlongraph=function(data,rho, ty,tole, lower,upper){
  alpha0=1
  omega0=diag(ncol(data)-2)
  subject=data[,1]
  subjectid=unique(subject)
  m=length(subjectid)
  p=ncol(data)-2
  tau_em=matrix(1/alpha0,nrow = 1,ncol = m)
  is=vector("list",m)
  tau1=rep(0,m)
  ll=rep(0,m)
  k=1
  while(1){
    print(paste("Iteration ",k,":", " Mean value of tau's","=", mean(tau1),sep=""))
    for (i in 1:m) {
      idata=data[which(data[,1]==subjectid[i]),]
      if (nrow(idata)==1){
        tau1[i]=NA
        is[[i]]=iss(idata = idata,itau = 1,ty = ty)
      }else{
        x=seq(lower,upper,length=50)
        z=sapply(x, logdensity,idata=idata,omega=omega0,alpha=alpha0,ty=ty)
        tau1[i]=x[min(which(z==min(z)))]
        ll[i]=min(z)
        is[[i]]=iss(idata = idata,itau = tau1[i],ty = ty)
      }
    }
    ll=sum(ll)
    tau_em=rbind(tau_em,tau1)
    s=Reduce("+",is)/nrow(data)
    omega1=glasso::glasso(s=s,rho = rho)$wi
    alpha1=1/mean(tau1,na.rm=T)
    if (max(abs(tau_em[nrow(tau_em),]-tau_em[nrow(tau_em)-1,]),na.rm = T)<tole){
      break
    }else{
      tau0=tau1
      omega0=omega1
      k=k+1
    }
  }
  result=list(alpha=alpha1,tau=tau1,omega=omega1,loglikehood=ll)
  return(result)
}







#' Title
#'
#' @param data
#' @param alpha0
#' @param omega
#' @param ty
#' @param tole
#' @param lower
#' @param upper
#'
#' @return
#' @export
#' @noRd
mle_alpha=function(data,alpha0,omega, ty, tole, lower,upper){
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
      z=sapply(x, logdensity,idata=idata,omega=omega,alpha=alpha0,ty=ty)
      tau1[i]=x[min(which(z==max(z)))]
      is[[i]]=iss(idata = idata,itau = tau1[i],ty = ty)
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




#' Title
#'
#' @param idata
#' @param omega
#' @param tau
#' @param alpha
#' @param ty
#' @export
#' @return
#' @noRd
logdensity=function(idata,omega,tau,alpha,ty){
  t=idata[,2]
  yy=as.matrix(idata[,-c(1,2)])
  y=c(t(scale(yy, scale = F)))
  p=nrow(omega)
  n=length(t)
  if (n*p!=length(y)){
    stop("the dimensions do not match")
  }

  phi=phifunction(t=t,tau = tau,ty = ty)
  a0=(n/2)*log(det(omega))
  a1=(-p/2)*log(det(phi))
  a2=-t(y)%*%kronecker(solve(phi),omega)%*%y/2
  a3=-alpha*tau
  a=a1+a2+a3+a0
  return(a)
}





#' Title
#'
#' @param idata
#' @param omega
#' @param tau
#' @param ty
#' @return
#' @export
#' @noRd
lli_homo=function(idata,omega,tau,ty){
  t=idata[,2]
  yy=as.matrix(idata[,-c(1,2)])
  y=c(t(scale(yy, scale = F)))
  p=nrow(omega)
  n=length(t)
  if (n*p!=length(y)){
    stop("the dimensions do not match")
  }

  phi=phifunction(t=t,tau = tau,ty = ty)
  a0=(n/2)*log(det(omega))
  a1=(-p/2)*log(det(phi))
  a2=-t(y)%*%kronecker(solve(phi),omega)%*%y/2
  a=a0+a1+a2
  return(a)
}


#' Title
#'
#' @param data
#' @param omega
#' @param tau
#' @param ty
#'
#' @return
#' @export
#' @noRd
ll_homo=function(data,omega,tau,ty){
  id=unique(data[,1])
  a=0
  for (i in id) {
    idata=data[which(data[,1]==i),]
    ai=lli_homo(idata = idata,omega = omega,tau = tau,ty = ty)
    a=a+ai
  }
  return(a)
}




#' Title
#'
#' @param data
#' @param rho
#' @param ty
#' @param tole
#' @param lower
#' @param upper
#'
#' @return
#' @export
#' @noRd
homolongraph=function(data, rho,ty, tole,lower,upper){
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
    z=sapply(x, ll_homo,data=data,omega=omega0,ty=ty)
    tau1=x[min(which(z==max(z)))]
    print(paste("Iteration ",k,":","tau=",tau1))
    for (i in 1:m) {
      idata=data[which(data[,1]==subjectid[i]),]
      is[[i]]=iss(idata = idata,itau = tau1,ty = ty)
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
  ll=ll_homo(data=data,omega = omega1,tau=tau1,ty=ty)
  result=list(tau=tau1,omega=omega1,loglikelihood=ll)
  return(result)
}


#' Title
#'
#' @param data
#' @param omega
#' @param ty
#' @param lower
#' @param upper
#'
#' @return
#' @export
#' @noRd
mle_tau=function(data, omega, ty,lower,upper){
  x=seq(lower,upper,length=50)
  z=sapply(x, ll_homo,data=data,omega=omega,ty=ty)
  tau1=x[min(which(z==max(z)))]
  return(tau1)
}




















#' Title
#'
#' @param data
#' @param priori
#'
#' @return
#' @export
#' @noRd
mle_net=function(data,priori){
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





#' Title
#'
#' @param data
#' @param rho
#' @param heter
#' @param ty
#' @param tole
#' @param lower
#' @param upper
#'
#' @return
#' @export
#'

lglasso=function(data, rho,heter=TRUE,ty=1, tole=0.01,lower=0.01,upper=10){
  if (heter==TRUE){
    aa=heterlongraph(data=data,rho=rho,ty=ty,tole=tole,lower=lower,upper=upper)
  }else{
    if (heter==FALSE){
      aa=homolongraph(data=data,rho=rho,ty=ty,tole=tole,lower=lower,upper=upper)
    }else{
      stop("Parameter heter only accept TRUE or FALSE!")
    }
  }
  return(aa)
}





#' Likelihood function
#'
#' @param data
#' @param network
#' @param heter
#' @param ty
#' @param tole
#' @param lower
#' @param upper
#'
#' @return
#' @export
mle=function(data,network,heter=TRUE,ty=1,tole=0.01,lower=0.01,upper=10){
  mlenetwork=mle_net(data = data[,-c(1,2)],priori = network)
  if (heter==TRUE){
    mle=mle_alpha(data=data,alpha0=1,omega=mlenetwork, ty=ty, tole=tole, lower=lower,upper=upper)
  }else{
    mle=mle_tau(data=data, omega=mlenetwork, ty=ty,lower=lower,upper=upper)
  }
  if (heter==TRUE){
    result=list(network=mlenetwork,alpha=mle$alpha,tau=mle$tau)
  }
  if (heter==FALSE){
    result=list(network=mlenetwork,tau=mle)
  }
  return(result)
}


