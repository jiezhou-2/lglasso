## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL








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

