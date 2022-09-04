bicfunction=function(data,mle,long=F,heter=F,x=NULL,type=1){
  n=nrow(data)
  p=ncol(data)-2
  id=unique(data[,1])
  M=mle$network
  k=length(which(M!=0))
##independent data
  if (long==F){
    yy=scale(as.matrix(data[,-c(1,2)]))
    f=function(x,M){t(x)%*%M%*%x}
    ll=-0.5*sum(apply(yy,1,function(x)  t(x) %*% M %*% x ))-0.5*n*log(det(2*pi*M))
    networkbic=-2*ll+k*log(n)
  }else{
  if (heter==F){
##longitudinal data: homo model
    ll=0
for (i in 1:length(id)) {
index=which(data[,1]==id[i])
idata=data[index,]
t=idata[,2]
yy=as.matrix(idata[,-c(1,2)])
y=c(t(scale(yy, scale = F)))
ni=length(t)
if (ni*p!=length(y)){
  stop("the dimensions do not match")
}
tau=mle$tau
phi=phifunction(t=t,tau = tau,type = type)
a1=(-p/2)*log(det(phi))+(n/2)*log(det(omega))
a2=-t(y)%*%kronecker(solve(phi),omega)%*%y/2
ll=ll+a1+a2
}
networkbic=-2*ll+k*log(n)
  }else{
  if (is.na(x)){

  }else{

  }
}
}
  return(bic=networkbic)
}
