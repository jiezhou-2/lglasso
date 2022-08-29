
#' Generate the heterogeneous data
#'
#' @param p An integer indicating dimension of the network
#' @param prob The density of the edges in the network
#' @param alpha A positive number representing the parameter for the exponential distribution
#' @param age A list with length equal to the number of subjects. Each component is a vector indicating the
#' time points at which the observations was made.
#' @return A list for the simulated data and its parameters
#' @export
sim_heter=function(p,prob,alpha,age){
  m=length(age)
  data=vector("list",m)
  tau=rexp(m,rate = alpha)
  K=BDgraph::bdgraph.sim(p=p,n=1,type="Gaussian",prob=prob)$K
  Sigma=solve(K)
  sqK=chol(Sigma)
  for (i in 1:m) {
    n=length(age[[i]])
    a=matrix(ncol = p,nrow = n)
    error1=matrix(rnorm(p*n),nrow = p)
    error2=t(t(sqK)%*%error1)
    a[1,]=error2[1,]
  for (t in 2:n) {
    coe=exp(-tau[i]*abs(age[[i]][t]-age[[i]][t-1]))
    a[t,]=a[t-1,]*coe+error2[t,]*sqrt(1-coe^2)
  }
    data[[i]]=cbind(i,age[[i]],a)
  }
  result=list(data=data,tau=tau,precision=K)
  return(result)
}



#' Generate the homogeneous data
#'
#' @param p An integer indicating dimension of the network
#' @param prob The density of the edges in the network
#' @param tau The shared correlation among observations.
#' @param age  A list with length equal to the number of subjects. Each component is a vector indicating the
#' time points at which the observations was made.
#' @return A list for simulated data
#' @export

sim_homo=function(p,prob,tau,age){
  m=length(age)
  data=vector("list",m)
  K=BDgraph::bdgraph.sim(p=p,n=1,type="Gaussian",prob=prob)$K
  Sigma=solve(K)
  sqK=chol(Sigma)
  for (i in 1:m) {
    n=length(age[[i]])
    a=matrix(ncol = p,nrow = n)
    error1=matrix(rnorm(p*n),nrow = p)
    error2=t(t(sqK)%*%error1)
    a[1,]=error2[1,]
    for (t in 2:n) {
      coe=exp(-tau*abs(age[[i]][t]-age[[i]][t-1]))
      a[t,]=a[t-1,]*coe+error2[t,]*sqrt(1-coe^2)
    }
    data[[i]]=cbind(i,age[[i]],a)
  }
  result=list(data=data,tau=tau,precision=K)
  return(result)
}


#' Comparison function
#'
#' @param real The real network
#' @param estimate The estimated network
#'
#' @return TPR and FPR
#' @export
comparison=function(real, estimate){
  real=real+t(real)
  diag(real)=1
  estimate=estimate+t(estimate)
  diag(estimate)=1
  N1=ifelse(abs(real)<=10^(-1),0,1)
  N2=ifelse(abs(estimate)<=10^(-2),0,1)
  if (any(dim(N1)!=dim(N2)))
    stop("Two matrixes should have the same dimension")
  p=dim(real)[1]
  real=(sum(N1)-p)/2
  null=p*(p-1)/2-real
  select=(sum(N2)-p)/2
  real_select=(sum(N2[N1==1])-p)/2
  fause_select=sum(N2[N1==0])/2
  if (real==0){return( list(TPR=1, FPR=0))}
  if (select==0) {return(list(TPR=0, FPR=1))}
  TPR=real_select/real
  FPR=fause_select/null
  #FPR=fause_select/select
  aa=as.numeric(c(TPR, FPR))
  return(aa)
}




#' Using regression to find the network matrix
#' @param data Matrix representing the n observations for p variable
#' @param lambda Tuning parameter for lasso
#'
#' @return Matrix representing the network
#' @export
addition=function(data,lambda){
  data=data[,-c(1,2)]
  p=ncol(data)
  n=nrow(data)
  arr=matrix(0,nrow=p,ncol=p)

  for (i in 1:p) {
    x=as.matrix(data[,-i])
    y=data[,i]
    result=glmnet::glmnet(x=x,y=y, family="gaussian",lambda = lambda)
    bb=as.matrix(result$beta)
    arr[i,-i]=bb
  }
  web=arr+t(arr)
  return(web)
}




