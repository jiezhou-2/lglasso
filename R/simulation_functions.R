




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
  print("heter data are generated")
  m=length(age)
  data=vector("list",m)
  tau=rexp(length(alpha),rate = alpha)
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
    colnames(data[[i]])[c(1,2)]=c("id","age")
  }
  result=list(data=data,tau=tau,alpha=alpha,precision=K)
  return(result)
}







#' Generate the heterogeneous data
#'
#' @param p An integer indicating dimension of the network
#' @param prob The density of the edges in the network
#' @param alpha A positive number representing the parameter for the exponential distribution
#' @param age A list with length equal to the number of subjects. Each component is a vector indicating the
#' time points at which the observations was made.
#' @return A list for the simulated data and its parameters
#' @export
sim_2heter=function(p,prob,alpha1,alpha2,age){
  print("two-community heter data are generated")
  p1=floor(p/2)
  result1=sim_heter(p=p1,prob = prob,alpha = alpha1,age=age)
  result2=sim_heter(p=p-p1,prob = prob,alpha = alpha2,age=age)
  data1=result1[[1]]
  K1=result1[[4]]
  data2=result2[[1]]
  K2=result2[[4]]
  K=bdiag(K1,K2)
  data=vector("list",length(data1))
  for (i in 1:length(data1)) {
  data[[i]]=merge(data1[[i]],data2[[i]],by=c("id","age"))
  }
  result=list(data=data,precision=K)
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
  print("homo data are generated")
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
    colnames(data[[i]])[1:2]=c("id","age")
  }
  result=list(data=data,tau=tau,precision=K)
  return(result)
}

#' Generate the two-community homogeneous data
#'
#' @param p An integer indicating dimension of the network
#' @param prob The density of the edges in the network
#' @param tau1 The shared correlation among observations in first community.
#' @param tau2 The correlation in second community
#' @param age  A list with length equal to the number of subjects. Each component is a vector indicating the
#' time points at which the observations was made.
#' @return A list for simulated data
#' @export

sim_2homo=function(p,prob,tau1,tau2,age){
  print("two-community homo data are generated")
  p1=floor(p/2)
  m=length(age)
  data1=vector("list",m)
  data2=vector("list",m)
  data=vector("list",m)
  sim1=sim_homo(p=p1, prob = prob,tau = tau1,age=age)
  sim2=sim_homo(p=p-p1,prob = prob,tau =tau2,age=age)
  data1=sim1[[1]]
  data2=sim2[[1]]
  K=bdiag(sim1[[3]],sim2[[3]])
  for (i in 1:m) {
  data[[i]]=merge(data1[[i]],data2[[i]],by=c("id","age"))
  }
  result=list(data=data,precision=K)
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
  N1=ifelse(abs(real)<=10^(-5),0,1)
  N2=ifelse(abs(estimate)<=10^(-5),0,1)
  if (any(dim(N1)!=dim(N2)))
    stop("Two matrixes should have the same dimension")
  p=nrow(real)
  real=(sum(N1)-p)/2
  null=p*(p-1)/2-real
  select=(sum(N2)-p)/2
  real_select=(sum(N2[N1==1])-p)/2
  fause_select=sum(N2[N1==0])/2
  if (real==0){
    stop("Real network has no edge. Please regenerate the real network")
    }
  if (select==0) {
    stop("Estimated network has no edge. Please decrease the tuning parameter.")
    }
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




