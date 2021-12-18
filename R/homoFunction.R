## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL



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

  result=list(tau=tau1,omega=omega1)
  return(result)
}


mle_tau=function(data, omega, ty,lower,upper){
    x=seq(lower,upper,length=50)
    z=sapply(x, ll_homo,data=data,omega=omega,ty=ty)
    tau1=x[min(which(z==max(z)))]
  return(tau1)
}











