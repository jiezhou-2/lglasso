## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL





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
      is[[i]]=iss(idata = idata,itau = tau1[i],ty = ty)
        }
    }
    tau_em=rbind(tau_em,tau1)
    s=Reduce("+",is)/nrow(data)
    omega1=glasso::glasso(s=s,rho = rho)$wi
    alpha1=1/mean(tau1,na.rm=T)
    if (max(abs(tau_em[nrow(tau_em),]-tau_em[nrow(tau_em)-1,]),na.rm = T)<tole){
      break
    }else{
      alpha0=alpha1
      tau0=tau1
      omega0=omega1
      k=k+1
    }

  }

  result=list(alpha=alpha1,tau=tau1,omega=omega1)
  return(result)
}







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
      tau1[i]=x[min(which(z==min(z)))]
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
  a1=(-p/2)*log(det(phi))
  a2=-t(y)%*%kronecker(solve(phi),omega)%*%y/2
  a3=-alpha*tau
  a=-(a1+a2+a3)
  return(a)
}



