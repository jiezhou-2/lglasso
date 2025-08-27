
#' Construct the temporal component fo correlation function
#'
#' @param t Time points of observations
#' @param tau length 2 vector representing the correlation parameter between time points
#' @param type The type of correlation function, which typically take either 0,1 or 2.
#' @author Jie Zhou
#' @export
#' @noRd
#' @return A square matrix with dimension equal to the length of vector t
phifunction=function(t,tau){
  n=length(t)
  #print(tau)
  if (any(tau<=0)){stop("tau should be positive!")}
  if (n==1){
    M=as.matrix(1)
  }else{
    M=matrix(nrow = n, ncol = n)
    for (i in 1:n) {
      for (j in i:n){
        M[i,j]=exp(-tau[1]*(abs(t[i]-t[j])))^tau[2]
        M[j,i]=M[i,j]
      }
    }
    diag(M)=1
  }
  return(M)
}




AA=function(B,data,type=c("general","expFixed","twoPara"),expFix,maxit=100,
            tol=10^(-4),lower=c(0.01,0.1),upper=c(10,5),...){
  ### clustered data
  type=match.arg(type)
  if (type=="general"){
    m3=length(unique(data[,2]))
    p=ncol(data)-2
    nn=length(unique(data[,1]))
    subjects=unique(data[,1])
    if (ncol(B)!=p | nrow(B)!=p){stop("Inputs do not match with each other!")}
    A=Variable(m3,m3,PSD=T) # tissue wise inverse correlation matrix
    obj1=(p*nn)/2*log_det(A)



    data_sub=split(data[,-c(1,2)],data[,1])
    amatrix=Reduce("+",lapply(data_sub, function(B,xx){as.matrix(xx)%*%B%*%t(as.matrix(xx))}, B=B))
    obj2=-0.5*matrix_trace(A%*%amatrix)
    obj=-(obj1+obj2)
    constr=list(CVXR::diag(A)==1)



    prob=Problem(Minimize(obj),constr)
    results=CVXR::solve(prob)
    corMatrix=solve(results$getValue(A))
    return(list(corMatrix=corMatrix))
  }


  ### longitudinal data

  if (type == "expFixed"){
    p=ncol(data)-2
    nn=length(unique(data[,1]))
    subjects=unique(data[,1])
    if (ncol(B)!=p | nrow(B)!=p){stop("Inputs do not match with each other!")}
    #YY=Variable(p,p,PSD=TRUE)
    time_list=split(data[,2],data[,1])
    subdata_list=split(data[,-c(1,2)],data[,1])

    likefun=function(tau){
      obj1=0
      obj2=0
      for (i in 1:nn) {

        datai=as.matrix(subdata_list[[i]])
        t=time_list[[i]]
        Aa=phifunction(t=t,tau = c(tau,expFix))
        amatrix=datai%*%B%*%t(datai)
        obji1=-0.5*p*log(det(Aa))
        obji2= -0.5*sum(diag(solve(Aa)%*%amatrix))
        obj1=obj1+obji1
        obj2=obj2+obji2
      }

      obj=-(obj1+obj2)
    }

    tau=optim(c(control$tau0[1]),likefun,method = "L-BFGS-B",lower = lower,upper = upper)$par
    A=lapply(time_list,phifunction,tau=c(tau,expFix))
    return(list(corMatrixList=A,tau=tau))
  }


  ### longitudinal data
  if (type == "twoPara"){
    p=ncol(data)-2
    nn=length(unique(data[,1]))
    subjects=unique(data[,1])
    if (ncol(B)!=p | nrow(B)!=p){stop("Inputs do not match with each other!")}
    #YY=Variable(p,p,PSD=TRUE)
    time_list=split(data[,2],data[,1])
    subdata_list=split(data[,-c(1,2)],data[,1])

    likefun=function(tau){
      obj1=0
      obj2=0
      for (i in 1:nn) {
        datai=as.matrix(subdata_list[[i]])
        t=time_list[[i]]
        Aa=phifunction(t=t,tau = c(tau,1))
        amatrix=datai%*%B%*%t(datai)
        obji1=-0.5*p*log(det(Aa))
        obji2= -0.5*sum(diag(solve(Aa)%*%amatrix))
        obj1=obj1+obji1
        obj2=obj2+obji2
      }

      obj=-(obj1+obj2)
    }

    tau=optim(c(control$tau0),likefun,method = "L-BFGS-B",lower = lower,upper = upper)$par
    A=lapply(time_list,phifunction,tau=tau)
    return(list(corMatrixList=A,tau=tau))
  }

}









BB=function(A,data,lambda,type=c("general","expFixed","twoPara"),diagonal=TRUE,maxit=100,
            tol=10^(-4),lower=c(0.01,0.1),upper=c(10,5), start=c("warm","cold"),w.init=NULL,wi.init=NULL,...){
  type=match.arg(type)
  start=match.arg(start)
  if (type=="general"){
    m3=length(unique(data[,2]))
    p=ncol(data)-2
    nn=length(unique(data[,1]))
    subjects=unique(data[,1])
    if (ncol(A)!=m3 | nrow(A)!=m3){
      stop("Inputs do not match with each other!")
    }

    data_sub=split(data[,-c(1,2)],data[,1])
    amatrix=Reduce("+",lapply(data_sub, function(A,xx){t(as.matrix(xx))%*%solve(A)%*%as.matrix(xx)}, A=A))

    if (start=="warm"){
      #w.init=diag(diag(cov(data[,-c(1,2)])))
      #wi.init=diag(1/diag(w.init))
    bb=glasso::glasso(s=amatrix/(m3*nn),rho=lambda[1], penalize.diagonal = diagonal, start = start,w.init = w.init,wi.init = wi.init)
    }else{
      bb=glasso::glasso(s=amatrix/(m3*nn),rho=lambda[1], penalize.diagonal = diagonal, start = "cold")
    }
    return(list(preMatrix=bb$wi, corMatrix=bb$w))
  }
  if (type == "expFixed"){
    p=ncol(data)-2
    nn=length(unique(data[,1]))
    subjects=unique(data[,1])

    data_sub=split(data[,-c(1,2)],data[,1])
    if (length(A)!=length(data_sub)){stop("Data do not match!")}
    bmate=0
    for (i in 1:length(A)) {
      xx=as.matrix(data_sub[[i]])
      yy=solve(as.matrix(A[[i]]))
      bmate=t(xx)%*%yy%*%xx+bmate
    }

    if (start=="warm"){
      #w.init=cov(data[,-c(1,2)])
      #wi.init=diag(ncol(w.init))
      bb=glasso::glasso(s=bmat/nrow(data),rho=lambda[1], penalize.diagonal = diagonal, start = start,w.init = w.init,wi.init = wi.init)
    }else{
      bb=glasso::glasso(s=bmate/nrow(data),rho=lambda[1], penalize.diagonal = diagonal, start = "cold")
    }

    return(list(preMatrix=bb$wi,corMatrix=bb$w))

  }

  if (type == "twoPara"){
    p=ncol(data)-2
    nn=length(unique(data[,1]))
    subjects=unique(data[,1])

    data_sub=split(data[,-c(1,2)],data[,1])
    if (length(A)!=length(data_sub)){stop("Data do not match!")}
    bmate=0
    for (i in 1:length(A)) {
      xx=as.matrix(data_sub[[i]])
      yy=solve(as.matrix(A[[i]]))
      bmate=t(xx)%*%yy%*%xx+bmate
    }

    if (start=="warm"){
     # w.init=cov(data[,-c(1,2)])
    #  wi.init=diag(ncol(w.init))
      bb=glasso::glasso(s=bmate/nrow(data),rho=lambda[1], penalize.diagonal = diagonal, start = start,w.init = w.init,wi.init = wi.init)
    }else{
      bb=glasso::glasso(s=bmate/nrow(data),rho=lambda[1], penalize.diagonal = diagonal, start = "cold")
    }


    return(list(preMatrix=bb$wi,corMatrix=bb$w))

  }

}


#' @title Longitudinal graphical lasso
#' @description
#'  This is the main function of the package, which identifies the underlying  network model from clustered data. Here
#'  clustered data include longitudinal data, or spatially correlated data,
#'  e.g, metabolites in different tissues of a same subject.  .
#'  or more broadly, clustered data for given tuning parameters.
#'
#'
#' @param data a \code{n} by \code{(p+2)} data frame in which the first column is subject ID, the second column is
#' the time point for longitudinal data or tissue ID.
#' @param lambda  a vector of length 1 or 2,  which
#' is the tuning parameter for the identification of the networks. For details, see the explanations in the below.
#' @param type a string specifying which model need to be fitted. There are three models available,
#' which are refered as \code{general, expFixed} and \code{twoPara} respectively.
#'  Please see the details in the below for the meaning of each.
#' @param expFix a nemeric number used in the model specification
#' @param group a vector  of length \code{n} if supplied which specify which data
#'  points need to be grouped together to infer a network for each of them.
#' @param maxit the maximum iterations for the estimation.
#' @param tol the minimum value for  convergence criterion
#' @param lower a vector of length 1 or 2 which specifies the lower bounds for the optimization algorithm
#' @param upper a vector of length 1 or 2 which specifies the upper bounds for the optimization algorithm
#' @param ... other inputs
#' @import glasso glasso
#' @import CVXR
#' @export
#'
lglasso=function(data,lambda, type=c("general","expFixed","twoPara"),expFix=1,group=NULL,diagonal=TRUE,maxit=100,
                 tol=10^(-4),lower=c(0.01,0.1),upper=c(10,5), start=c("cold","warm"), w.init=NULL, wi.init=NULL,trace=FALSE,
                 ...)

  {


if (is.null(group) & length(lambda)!=1 | !is.null(group) & length(lambda) !=2) {
  stop("Arguments (group, lambda) do not match!")
}

  if (!all(lambda>0)){
    stop("lambda must be positive!")
  }




  type=match.arg(type)
  start=match.arg(start)


  X_bar = apply(data[,-c(1,2)], 2, mean)
  data[,-c(1,2)] = scale(data[,-c(1,2)], center = X_bar, scale = FALSE)



  if (type=="general"){
    m3=length(unique(data[,2]))
    p=ncol(data)-2
    A=diag(m3)
    B=diag(p)
    k=0

    while(1){
      k=k+1

      A1=AA(data = data,B = B, type=type,...)$corMatrix
#browser()
      B1=BB(data=data,A=A,lambda = lambda, type=type, start = start, w.init=w.init,wi.init=wi.init,...)

      d1=round(max(abs(B-B1$preMatrix)),3)
      d2=round(max(abs(A-A1)),3)
      if (trace){
      print(paste0("iteration ",k, " precision difference: ",d1 , " /correlation difference: ",d2))
      }
      if (d1<=tol & d2<= tol ){
        if (is.null(group)){
          return(list(w=B1$corMatrix,wi=B1$preMatrix, v=A1))
          }else{
          individial=heternetwork(data = data,lambda = lambda[2],homo = B1$preMatrix, group=group)
          return(list(wi=B1$preMatrix,  wiList=individial,v=A1))
        }

      }else{
        A=A1
        B=B1$preMatrix
      }
      if (k>=maxit){
        warning("Algorithm did not converge!")

        if (is.null(group)){
          return(list(w=B1$corMatrix,wi=B1$preMatrix, v=A1))
        }else{
          individial=heternetwork(data = data,lambda = lambda[2],homo = B1$preMatrix, group=group)
          return(list(wi=B1$preMatrix,wiList=individial, v=A1))
        }

      }
    }
  }




  if (type == "expFixed"){
    if (is.null(expFix) | length(expFix)!=1 | !is.numeric(expFix)){stop("Argument expFix is not correctly specified!")}
    p=ncol(data)-2
    fre=table(data[,1])
    A=vector("list",length(fre))
    for (i in 1:length(A)) {
      A[[i]]=diag(fre[i])
    }
    B=diag(p)
    k=0
tau0=control$tau0
    while(1){
      k=k+1

      #A1=AA(data = data,B = B, type=type,fix=expFix, init=tau0,lower = lower,upper = upper)
      A1=AA(data = data,B = B, type=type,expFix,...)
      B1=BB(data=data,A=A,lambda = lambda[1], type=type,start=start, w.init=w.init,wi.init=wi.init, ...)

      d1=round(max(abs(B-B1$preMatrix)),4)
      d2=round(max(abs(tau0-A1$tau)),4)
      if (trace){
      print(paste0("iteration ",k, " precision difference: ",d1 , " /correlation tau difference: ",d2))
      }
      if (d1<=tol & d2<= tol ){

        if (is.null(group)){
          return(list(wi=B1$preMatrix, vList=A1$corMatrixList, tauhat=tau0))
        }else{

          individial=heternetwork(data = data,lambda = lambda[2],homo = B1$preMatrix, group=group)

          return(list(wi=B1$preMatrix, wiList=individial, vList=A1$corMatrixList, tauhat=tau0))
        }



      }else{
        A=A1$corMatrixList
        B=B1$preMatrix
        tau0=A1$tau
      }
      if (k>=maxit){
        warning("Algorithm did not converge!")
        if (is.null(group)){
          return(list(wi=B1$preMatrix, vList=A1$corMatrixList, tauhat=tau0))
        }else{
          #browser()
          individial=heternetwork(data = data,lambda = lambda[2],homo = B1$preMatrix, group=group)
          return(list(wi=B1$preMatrix, wiList=individial, vList=A1$corMatrixList, tauhat=tau0))
        }
      }
    }
  }






  if (type == "twoPara"){
    p=ncol(data)-2
    fre=table(data[,1])
    A=vector("list",length(fre))
    for (i in 1:length(A)) {
      A[[i]]=diag(fre[i])
    }
    B=diag(p)
    k=0
tau0=control$tau0
    while(1){
      k=k+1
      A1=AA(data = data,B = B, type=type, fix=fix)
      B1=BB(data=data,A=A,lambda = lambda[1], type=type,start=start, w.init=w.init,wi.init=wi.init,...)

      d1=round(max(abs(B-B1$preMatrix)),4)
      d2=round(max(abs(tau0-A1$tau)),4)
      if (trace){
      print(paste0("iteration ",k, " precision difference: ",d1 , " /correlation tau difference: ",d2))
      }
      if (d1<= tol & d2<= tol ){

        if (is.null(group)){
          return(list(wi=B1$preMatrix, vList=A1$corMatrixList, tauhat=tau0))
        }else{

          individial=heternetwork(data = data,lambda = lambda[2],homo = B1$preMatrix, group=group)

          return(list(wi=B1$preMatrix, wiList=individial, vList=A1$corMatrixList, tauhat=tau0))
        }

      }else{
        A=A1$corMatrixList
        B=B1$preMatrix
        tau0=A1$tau
      }
      if (k>= maxit){
        warning("Algorithm did not converge!")

        if (is.null(group)){
          return(list(wi=B1, vList=A1$corMatrixList, tauhat=tau0))
        }else{
          #browser()
          individial=heternetwork(data = data,lambda = lambda[2],homo = B1$preMatrix, group=group)
          return(list(wi=B1$preMatrix,wiList=individial, vList=A1$corMatrixList, tauhat=tau0))
      }
    }
    }

  }
}





heternetwork=function(data,lambda,homo, group){
  if (length(group)!=nrow(data)){
    stop("Argument group doens not match data!")
    }
  m3=length(unique(group))
  data_list=split(data,group)
  p=ncol(data)-2
 S=vector("list",m3)
 Q=vector("list",m3)
  obj=0
  aa=0
for (i in 1:m3) {
  S[[i]]= Variable(p,p,PSD=TRUE)
  Q[[i]]=cov(data_list[[i]][,-c(1,2)])
  obj=obj+log_det(S[[i]])-matrix_trace(S[[i]]%*%Q[[i]])
  aa=aa+ sum(abs(S[[i]]-homo))
}
constr=list(aa<=1/lambda)

prob=Problem(Maximize(obj),constr)
result=CVXR::solve(prob,solver="SCS")
S_est= lapply(S, function(x) result$getValue(x))
return(individual=S_est)
}


llike=function(data,A, preMatrix, individual){

}




cvErrori=function(data,bi, K){
  n=nrow(data)
  p=ncol(data)-2
  ind = sample(n)


  if (K<=1 | K%%1 !=0){
    stop("K should be an integer greater than 1!")
  }

  #if (any(lambda)<=0) {stop("tuning parameter lambda should be positive!")}

  if (any(! bi %in% c(0,1,2))) {stop("entries of vector bi should be 0 or 1!")}
  i=which(bi==2)
  X=data[,-c(1,2)]
  cv_error=rep(0,length=K)


  for (k in 1:K) {

    leave.out = ind[(1 + floor((k - 1) * n/K)):floor(k *
                                                       n/K)]
    X.train = X[-leave.out, , drop = FALSE]
    X_bar = apply(X.train, 2, mean)
    X.train = scale(X.train, center = X_bar, scale = FALSE)
    X.valid = X[leave.out, , drop = FALSE]
    X.valid = scale(X.valid, center = X_bar, scale = FALSE)


    index=which(bi==1)
    if (length(index)==0){
      cv_error[k]=var(X.valid[,i])

    }else{
      y=X.train[,i]
      x=X.train[,index]
      coef.train=lm(y~x)$coef
      yy=X.valid[,i]
      if(length(leave.out)>1){
        xx=cbind(1,X.valid[,index])
      }else{
        xx=c(1,X.valid[index])}
      #browser()
      cv_error[k]=mean(yy-xx%*%coef.train)^2
      if (is.na(cv_error[k])) {browser()}
    }
  }
  return(mean(cv_error))
}





