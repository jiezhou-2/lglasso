
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




AA=function(B,data,type=c("general","expFixed","twoPara"),expFix,maxit=30,
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
    #browser()
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
#' @param lambda   vector of length 1 or 2,  which
#' is the tuning parameter for the identification of the networks. For details, see the explanations in the below.
#' @param type a string specifying which model need to be fitted. There are three models available,
#' which are refered as \code{general, expFixed} and \code{twoPara} respectively.
#'  Please see the details in the below for the meaning of each.
#' @param expFix  numeric number used in the model specification
#' @param group  vector  of length \code{n} if supplied which specify which data
#'  points need to be grouped together to infer the heterogeneous networks for, e.g, pre/post vaccination.
#' @param maxit the maximum iterations for the estimation.
#' @param tol the minimum value for  convergence criterion
#' @param lower  vector of length 1 or 2 which specifies the lower bounds for (\alpha_1, \alpha_2) in the optimization algorithm
#' @param upper  vector of length 1 or 2 which specifies the upper bounds for (\alpha_1, \alpha_2) in the optimization algorithm
#' @param ... other inputs
#' @import glasso
#' @import CVXR
#' @export
#' @return  \code{w} the general covariance matrix estimate
#' @return  \code{wList} list representing the individual covariance matrix estimate
#' @return \code{wi} the general precision matrix estimate
#' @return  \code{wiList} list representing the individual precision matrix estimate
#' @return \code{v} the correlation matrix between specified classes
#' @return \code{vList} list representing the individual correlation matrix
#' @return \code{tauhat} the correlation parameters for longitudinal data
#' @details
#' This function implements three statistical models for  network inference, according to how the
#'  correlations is specified between time points (or tissues or contents in some clinical studies).
#'  These three models are referred as
#' \code{general, expFixed} and \code{twoPara}. Let's say we have two time points,\code{t_i,t_j},
#' then in model \code{general}, the correlation is \tau_{ij}, while in model $expFixed$, we have
#' \tau=exp(-\alpha_1*|t_1-t_2|^{-\alpha_2}) with \alpha_2 is pre-specified (default is \alpha_2=1).
#' In model $twoPara$, both \alpha_1 and \alpha_2 is unknown and need to be inferred from the data.
#' For longitudinal data, model $expFixed$ is recommended while for omic data from different tissues or contents,
#' model \code{general} should be adopted.
lglasso=function(data,lambda, type=c("general","expFixed","twoPara"),expFix=1,group=NULL,diagonal=TRUE,maxit=30,
                 tol=10^(-3),lower=c(0.01,0.1),upper=c(10,5), start=c("cold","warm"), w.init=NULL, wi.init=NULL,trace=FALSE,
                 ...)

  {


if (is.null(group))  {
 if (length(lambda)!=1){
  stop("Arguments (group, lambda) do not match!")
 }
}

  if (!is.null(group))  {
    if (any(!is.matrix(lambda) | !ncol(lambda)!=2)){
      stop("Arguments (group, lambda) do not match!")
    }
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
      B1=BB(data=data,A=A,lambda = lambda, type=type, start = start, w.init=w.init,wi.init=wi.init,...)

      d1=round(max(abs(B-B1$preMatrix)),3)
      d2=round(max(abs(A-A1)),3)


      if (trace){
      print(paste0("iteration ",k, " precision difference: ",d1 , " /correlation difference: ",d2))
      }

      #browser()
      if (d1<=tol && d2<= tol ){
        if (is.null(group)){
          output=structure(list(w=B1$corMatrix,wi=B1$preMatrix, v=A1), class="lglasso")
          }else{
          individial=heternetwork(data = data,lambda = lambda[2],homo = B1$preMatrix, group=group)
          output=structure(list(wi=B1$preMatrix,  wiList=individial,v=A1), class="lglasso")
        }
          break
      }else{
        A=A1
        B=B1$preMatrix
      }


      if (k>=maxit){
        warning("Algorithm did not converge!")

        if (is.null(group)){
          output=structure(list(w=B1$corMatrix,wi=B1$preMatrix, v=A1), class="lglasso")
        }else{
          individial=heternetwork(data = data,lambda = lambda[2],homo = B1$preMatrix, group=group)
          output=structure(list(wi=B1$preMatrix,wiList=individial, v=A1), class="lglasso")
        }
        break
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
          output=structure(list(wi=B1$preMatrix, wiList=NULL, vList=A1$corMatrixList, tauhat=tau0), class="lglasso")
        }else{

          individial=heternetwork(data = data,lambda = lambda[2],homo = B1$preMatrix, group=group)
          output=structure(list(wi=B1$preMatrix, wiList=individial, vList=A1$corMatrixList, tauhat=tau0), class="lglasso")
        }

        break

      }else{
        A=A1$corMatrixList
        B=B1$preMatrix
        tau0=A1$tau
      }
      if (k>=maxit){
        warning("Algorithm did not converge!")
        if (is.null(group)){

          output=structure(list(wi=B1$preMatrix, wiList=NULL, vList=A1$corMatrixList, tauhat=tau0), class="lglasso")

        }else{
          individial=heternetwork(data = data,lambda = lambda[2],homo = B1$preMatrix, group=group)
          output=structure(list(wi=B1$preMatrix, wiList=individial, vList=A1$corMatrixList, tauhat=tau0), class="lglasso")
        }
        break
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
          output=structure(list(wi=B1$preMatrix, wList=NULL,vList=A1$corMatrixList, tauhat=tau0), class="lglasso")
        }else{

          individial=heternetwork(data = data,lambda = lambda[2],homo = B1$preMatrix, group=group)
          output=structure(list(wi=B1$preMatrix, wiList=individial, vList=A1$corMatrixList, tauhat=tau0), class="lglasso")
        }
        break
      }else{
        A=A1$corMatrixList
        B=B1$preMatrix
        tau0=A1$tau
      }
      if (k>= maxit){
        warning("Algorithm did not converge!")

        if (is.null(group)){
          output=structure(list(wi=B1$preMatrix, wList=NULL,vList=A1$corMatrixList, tauhat=tau0), class="lglasso")

        }else{
          #browser()
          individial=heternetwork(data = data,lambda = lambda[2],homo = B1$preMatrix, group=group)
          output=structure(list(wi=B1$preMatrix, wiList=individial, vList=A1$corMatrixList, tauhat=tau0), class="lglasso")
      }
        break
          }

    }

  }
  return(output)
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
names(S_est)=names(data_list)
return(individual=S_est)
}









cvErrorji=function(data.train,data.valid,bi){

  if (any(! bi %in% c(0,1,2))) {stop("entries of vector bi should be 0,  1 or 2!")}
  i=which(bi==2)
  cv_error=c()


    index=which(bi==1)
    if (length(index)==0){
      cv_error=var(data.valid[,i+2])

    }else{
      if (length(index)>= nrow(data.train)){
        print(paste("number of variable ", length(index)))
        stop("network is too dense for model training!")
      }

      y=data.train[,i+2]
      x=as.matrix(data.train[,index+2])
      #browser()
      coef.train=lm(y~x)$coef
      yy=data.valid[,i+2, drop=FALSE]
      # if(nrow(data.valid)>1 && ncol(data.valid>1)){
      #   xx=cbind(1,X.valid[,index])
      # }else{
      #   xx=c(1,c(X.valid[index]))}

      xx=as.matrix(cbind(1,data.valid[,index+2,drop=FALSE]))
      err=(yy-xx%*%coef.train)^2
      cv_error=c(cv_error,mean(err[,,drop=T]))
      if (any(is.na(cv_error))) {
        print("cv error is missing!")
        browser()}
    }

  return(mean(cv_error))
}



cvErrorj=function(data.train,data.valid,B){
    a= mean(apply(B, 2, function(bi) cvErrorji(data.train=data.train,data.valid=data.valid,bi=bi)))
}




#' Title
#'
#' @param data.train trainng data
#' @param data.valid testing data
#' @param BB given network
#' @param group.train group in training data
#' @param group.valid group in testing data
#'
#' @returns a matrix
#'
cvError=function(data.train,data.valid,B,group.train=NULL,group.valid=NULL){
#browser()

  if (is.matrix(B)){
    a=  cvErrorj(data.train=data.train,data.valid=data.valid,B=B)

  }

  if (is.list(B)){

  if (any(nrow(data.train)!=length(group.train) | nrow(data.valid)!=length(group.valid) )){
    stop("group does not match dat sets!")
  }

  data.train.sub=split(data.train,group.train)
  data.valid.sub=split(data.valid,group.valid)


if (any(names(data.train.sub)!=names(B)) | any(names(data.valid.sub)!= names(B))) {
  stop("the names of data sets do not match!")
}
a=c()
  for (i in 1:length(B)) {
    dd1=data.train.sub[[i]]
    dd2=data.valid.sub[[i]]
    Bi=B[[i]]
    a=c(a,cvErrorj(data.train=dd1,data.valid=dd2,B=Bi))

}
}
  return(a=mean(a))
}


#' Title
#'
#' @param data data frame
#' @param B given network
#' @param K number of cross validation
#'
#' @returns list



cvlglasso=function(type=c("general","expFixed","twoPara"), data,group=NULL,
                   lambda=NULL,nlam=10,lam.min.ratio=0.01, K, expFix=1,trace=FALSE){

  type=match.arg(type)


  if (is.null(group) && any( !is.vector(lambda) |  !all(is.numeric(lambda)) | !all(lambda>0)))
  {stop("group and lambda does not match!")}

  if (ncol(lambda)!=2 && !is.null(group)){stop("lambda should be a n by 2 matrix when group is specified!")}

  if (any(lambda<=0)) {stop("tuning parameter lambda should be positive!")}

  if (any(K<=1 | K%%1 !=0)){
    stop("K should be an integer greater than 1!")
  }


  n=length(unique(data[,1]))
  subjects=unique(data[,1])
  p=ncol(data)-2
  ind = sample(n)

  X=data[,-c(1,2)]




  S = (nrow(X) - 1)/nrow(X) * cov(X)
  # crit.cv = match.arg(crit.cv)
  # start = match.arg(start)

  Sminus = S
  diag(Sminus) = 0
  if (is.null(lambda)) {
    if (!((lam.min.ratio <= 1) && (lam.min.ratio > 0))) {
      cat("\nlam.min.ratio must be in (0, 1]... setting to 1e-2!")
      lam.min.ratio = 0.01
    }
    if (!((nlam > 0) && (nlam%%1 == 0))) {
      cat("\nnlam must be a positive integer... setting to 10!")
      nlam = 10
    }
    lam.max = max(abs(Sminus))
    lam.min = lam.min.ratio * lam.max
    lambda = 10^seq(log10(lam.min), log10(lam.max), length = nlam)
    if (!is.null(group)){
      lambda=cbind(lambda,lambda)
    }
  }
  else {
    if (is.null(group)){
      lambda = sort(lambda)
      }
  }


nnlambda=ifelse(is.null(group),length(lambda),nrow(lambda))
cv_error=matrix(0,nrow=nnlambda,ncol=K)
  if (trace) {
    progress = txtProgressBar(max = K, style = 3)
  }

  for (k in 1:K) {
      leave.out =subjects[ind[(1 + floor((k - 1) * n/K)):floor(k *
                                                         n/K)]]
      indexValid=which(data[,1] %in% leave.out)
      data.train = data[-indexValid, , drop = FALSE]
      data_bar = apply(data.train[,-c(1,2)], 2, mean)
      data.train[,-c(1,2)] = scale(data.train[,-c(1,2)], center = data_bar, scale = FALSE)
      data.valid = data[indexValid,, drop = FALSE]
      data.valid[,-c(1,2)] = scale(data.valid[,-c(1,2)], center = data_bar, scale = FALSE)
      group.train=group[-indexValid]
      group.valid=group[indexValid]
      #S.train = crossprod(data.train[,-c(1,2)])/(dim(data.train)[1])
      #S.valid = crossprod(data.valid[,-c(1,2)])/(dim(data.valid)[1])




    if (type=="general" && is.null(group)){
      aa= sapply(lambda, function(x) lglasso(data=data.train,lambda=x,type=type)$wi, simplify = FALSE)
      cc=lapply(aa, function(B){
        M=ifelse(abs(B)<=10^(-2), 0,1)
        diag(M)=2
        M
      })
      bb=lapply(cc, function(B) cvError(data.train=data.train,data.valid=data.valid,B=B))
      bb=do.call(c,bb)
    }



  if (type=="general" && !is.null(group)){

    aa=apply(lambda, 1,function(x){lglasso(data=data.train,lambda=x,type=type, group=group.train)$wiList} )


    cc=lapply(aa, function(B){
      lapply(B, function(Z){
        M=ifelse(abs(Z)<=10^(-1), 0,1)
             diag(M)=2
             M
             }
             )
    }
    )

    bb=lapply(cc, function(BB){
      cvError(data.train=data.train,data.valid=data.valid,B=BB, group.valid=group.valid, group.train = group.train)
    }
    )
    bb=do.call(c,bb)
}


  if (type=="expFixed" && is.null(group)){

   aa= sapply(lambda, function(x) lglasso(data=data.train,lambda=x,type=type, expFix = expFix)$wi, simplify = FALSE)
   cc=lapply(aa, function(B){
     M=ifelse(abs(B)<=10^(-2), 0,1)
     diag(M)=2
     M
   })
   bb=lapply(cc, function(B) cvError(data.train=data.train,data.valid=data.valid,B=B))
   bb=do.call(c,bb)
  }

  if (type=="expFixed" && !is.null(group)){
   aa= apply(lambda,1, function(x) {lglasso(data=data.train,lambda=x,type=type, expFix = expFix, group=group.train)$wiList})
   cc=lapply(aa, function(B){
     lapply(B, function(Z){
       M=ifelse(abs(Z)<=10^(-1), 0,1)
       diag(M)=2
       M
     }
     )
   }
   )

   bb=lapply(cc, function(BB) cvError(data.train=data.train,data.valid=data.valid,B=BB,
                                      group.valid =   group.valid,group.train = group.train))
   bb=do.call(c,bb)
  }

  if (type=="twoPara" && is.null(group)){
    aa= sapply(lambda, function(x) lglasso(data=data.train,lambda=x,type=type)$wi, simplify = FALSE)
    cc=lapply(aa, function(B){
      M=ifelse(abs(B)<=10^(-2), 0,1)
      diag(M)=2
      M
    })
    bb=lapply(cc, function(B) cvError(data.train=data.train,data.valid=data.valid,B=B))
    bb=do.call(c,bb)
  }

  if (type=="twoPara" && !is.null(group)){
   # browser()
    aa= apply(lambda,1, function(x) {lglasso(data=data.train,lambda=x,type=type, group=group.train)$wiList})
    cc=lapply(aa, function(B){
      lapply(B, function(Z){
        M=ifelse(abs(Z)<=10^(-1), 0,1)
        diag(M)=2
        M
      }
      )
    }
    )

    bb=lapply(cc, function(BB) cvError(data.train=data.train,data.valid=data.valid,B=BB,
                                       group.valid =   group.valid,group.train = group.train))
    bb=do.call(c,bb)
  }
cv_error[,k]=bb
if (trace) {
  setTxtProgressBar(progress,  k)
}
  }


output=structure(list(cv_error=cv_error, lambda=lambda),class="cvlglasso")
return(output)
}




#' Title
#'
#' @param x CVlglasso object
#' @param xvar character which specify the x axis of the plot
#' @param ...
#'
#' @returns If \code{group} is NULL in \code{CVlglasso}, then a line plot will produced; otherwise, a heatmap will be produced.
#' @export
#'
plot.cvlglasso=function(x, xvar=c("lambda","step"),...){
  xvar=match.arg(xvar)
  if (!inherits(x, "cvlglasso")) {
    stop("x must be an object of class 'cvlglasso'")
  }

  if (xvar == "lambda") {
    xlab_label <- "Lambda"
    x_data <- x$lambda # Assuming your cvlglasso object has a lambda component
  } else if (xvar == "step") {
    xlab_label <- "Steps"
    x_data <- seq_along(x$lambda) # Assuming lambda can represent steps
  }



  lambda=x$lambda
  err=apply(x$cv_error, 1, mean)
  if (is.vector(lambda)){
    graphics::plot(x=x_data,y=err,xlab=xlab_label,ylab="CV Error", main = "cvlglasso Fit",
                   type = "b", ...)
  }else{
lambda=as.matrix(lambda)
a1=unique(lambda[,1])
a2=unique(lambda[,2])
err_matrix=matrix(0,nrow=length(a1),ncol=length(a2),dimnames = list(round(a1,3), round(a2,3)))
for (i in 1:length(a1)) {
  for (j in 1:length(a2)) {
    index=which(apply(lambda,1,function(x) all(x==c(a1[i],a2[j]))))
    err_matrix[i,j]=err[index]
  }
}



  # heat_plot <- pheatmap::pheatmap(err_matrix,
  #                       col = brewer.pal(8, 'OrRd'), # choose a colour scale for your data
  #                       cluster_rows = F, cluster_cols = F, # set to FALSE if you want to remove the dendograms
  #                       clustering_distance_cols = 'euclidean',
  #                       clustering_distance_rows = 'euclidean',
  #                       clustering_method = 'ward.D',
  #                       #annotation_row = gene_functions_df, # row (gene) annotations
  #                       #annotation_col = ann_df, # column (sample) annotations
  #                       #annotation_colors = ann_colors, # colours for your annotations
  #                       #annotation_names_row = F,
  #                       #annotation_names_col = F,
  #                       fontsize_row = 10,          # row label font size
  #                       fontsize_col = 7,          # column label font size
  #                       angle_col = 45, # sample names at an angle
  #                       legend_breaks = c(-2, 0, 2), # legend customisation
  #                       legend_labels = c("Low", "Medium", "High"), # legend customisation
  #                       show_colnames = T, show_rownames = F, # displaying column and row names
  #                       main = "CV error") # a title for our heatmap

heat_plot <- pheatmap::pheatmap(
  err_matrix,
  col = RColorBrewer::brewer.pal(8, 'OrRd'),
  cluster_rows = FALSE, cluster_cols = FALSE,
  main = "CV error",
  # `angle_col` can be used to rotate column labels for readability
  angle_col = 45,
  # You can also customize label font sizes if needed
  fontsize_row = 8,
  fontsize_col = 8,
  annotation_names_row = F,
  annotation_names_col = F,
  ...
)

return(invisible(heat_plot))

}

}



simulate_general=function(n,p,m1,m2=0,m3,cc){
  ## true structure
  if (ncol(cc)!=m3){stop("Unmatched inputs!")}
  real_stru=matrix(0, nrow = p, ncol = p)
  real_stru[lower.tri(real_stru,diag = T)]=1
  index=which(real_stru==0,arr.ind = T)
  a=sample(1:nrow(index),m1, replace = F)
  real_stru[index[a,]]=1
  real_stru[lower.tri(real_stru,diag = T)]=0
  real_stru=real_stru+t(real_stru)+diag(p)
  mu=rep(0,p*m3)
  Precision=list()
  Structure=list()
  Sigma=list()
  Data=list()
  ## tissue specific structure

  if (m2==0){
          ## generate the precision matrices
          theta = matrix(rnorm(p^2,mean = 0,sd=2), ncol = p,nrow = p)
          theta = theta + t(theta)
          diag(theta)=1
          ### apply the required sparsity
          theta = theta * real_stru

          # force it to be positive definite
          theta=MakePositiveDefinite(theta,pd_strategy = "diagonally_dominant",scale = T)$omega
          Precision=theta
          Sigma=kronecker(cc,solve(theta))
    fullCovariance=Sigma
  }

  if (m2!=0){
    for (i in 1:m3) {
      Sigma[[i]]=list()
      Precision[[i]]=list()
      Structure[[i]]=list()
      for (j in 1:m3) {
        if (j<i){
          Precision[[i]][[j]]=Precision[[j]][[i]]
          Sigma[[i]][[j]]=Sigma[[j]][[i]]
          next()
        }else{
          distrubance=matrix(0, nrow = p, ncol = p)
          distrubance[lower.tri(distrubance,diag = T)]=1
          index=which(distrubance==0,arr.ind = T)
          a=sample(1:nrow(index),m2, replace = F)
          distrubance[index[a,]]=1
          distrubance[lower.tri(distrubance,diag = T)]=0
          distrubance=distrubance+t(distrubance)+diag(p)
          prior_stru=(real_stru+distrubance)%%2

          diag(prior_stru)=1
          #colnames(prior_stru)=paste0("metabolite",1:p)
          #rownames(prior_stru)=paste0("metabolite",1:p)
          Structure[[i]][[j]]=prior_stru
          ## generate the precision matrices
          theta = matrix(rnorm(p^2,mean = 0,sd=2), ncol = p,nrow = p)
          #theta[lower.tri(theta, diag = TRUE)] = 0
          theta = theta + t(theta)
          diag(theta)=1
          ### apply the required sparsity
          theta = theta * prior_stru

          # force it to be positive definite
          theta=MakePositiveDefinite(theta,pd_strategy = "diagonally_dominant",scale = T)$omega
          colnames(theta)=paste0("metabolite",1:p)
          rownames(theta)=paste0("metabolite",1:p)
          Precision[[i]][[j]]=theta
          Sigma[[i]][[j]]=solve(theta)*cc[i,j]
        }
      }
    }
    mmlist=list()
    for (k in 1:m3) {
      mmlist[[k]]=do.call(cbind,Sigma[[k]])
    }
    fullCovariance=do.call(rbind,mmlist)
    fullCovariance=MakePositiveDefinite(fullCovariance,pd_strategy = "diagonally_dominant",scale = T)$omega

  }


  data=mvrnorm(n,mu=mu,Sigma = fullCovariance)
  fulldata=c()

  for (k in 1:nrow(data)) {
    a=c()
    for (i in 1:m3) {
      a=as.data.frame(rbind(a,data[k,c(((i-1)*p+1) :(i*p))]))
    }
    subject=rep(paste0("subject",k),m3)
    tissue=paste0("tissue",1:m3)
    a=cbind(subject,tissue,a)
    fulldata=rbind(fulldata,a)
  }
  return(list(data=fulldata, fullCovariance=fullCovariance,corePrecision=theta))
}







simulate_long=function(n,p,m1,tau){
  ## true structure
  timepoint=vector("list",n)

  for (i in 1:n) {
    m3=sample(x=1:15,1,prob = rep(1,1,length=15))
    t1=rexp(m3)
    timepoint[[i]]=cumsum(t1)
  }
  cc=lapply(timepoint,phifunction, tau=tau)
  real_stru=matrix(0, nrow = p, ncol = p)
  real_stru[lower.tri(real_stru,diag = T)]=1
  index=which(real_stru==0,arr.ind = T)
  a=sample(1:nrow(index),m1, replace = F)
  real_stru[index[a,]]=1
  real_stru[lower.tri(real_stru,diag = T)]=0
  real_stru=real_stru+t(real_stru)+diag(p)
  Precision=list()
  Sigma=list()



  mmlist=list()
  theta = matrix(rnorm(p^2,mean = 0,sd=2), ncol = p,nrow = p)
  theta[lower.tri(theta, diag = TRUE)] = 0
  theta = theta + t(theta) + diag(p)
  theta = theta * real_stru
  theta=MakePositiveDefinite(theta,pd_strategy = "diagonally_dominant",scale = T)$omega
  #colnames(theta)=paste0("metabolite",1:p)
  #rownames(theta)=paste0("metabolite",1:p)
  sigma=solve(theta)
  fullCovariance= lapply(cc,function(x,cc) {kronecker(cc,x)},x=sigma)
  fulldata=c()
  a=c()
  for (i in 1:n) {
    ai=c()
    m3=length(timepoint[[i]])
    mu=rep(0,m3*p)
    data=MASS::mvrnorm(1,mu=mu,Sigma = fullCovariance[[i]])
    for (j in 1:m3) {
      ai=as.data.frame(rbind(ai,data[c(((j-1)*p+1) :(j*p))]))
    }
    subject=rep(paste0("subject",i),m3)
    ai=cbind(subject,timepoint[[i]],ai)
    a=rbind(a,ai)
  }

  return(list(data=a,B=theta,timepoint=timepoint))
}





simulate_heter=function(n,p,m1,alpha){
  ## true structure
  timepoint=vector("list",n)
  tau=c()
  cc=vector("list",)
  for (i in 1:n) {
    m3=sample(x=1:15,1,prob = rep(1,1,length=15))
    t1=rexp(m3)
    timepoint[[i]]=cumsum(t1)
    tau=c(tau,rexp(1,rate=alpha))
  }
  for (j in 1:length(tau)) {
    cc[[j]]=phifunction(t=timepoint[[j]],tau=c(tau[j],1))
  }

  real_stru=matrix(0, nrow = p, ncol = p)
  real_stru[lower.tri(real_stru,diag = T)]=1
  index=which(real_stru==0,arr.ind = T)
  a=sample(1:nrow(index),m1, replace = F)
  real_stru[index[a,]]=1
  real_stru[lower.tri(real_stru,diag = T)]=0
  real_stru=real_stru+t(real_stru)+diag(p)
  Precision=list()
  Sigma=list()



  mmlist=list()
  theta = matrix(rnorm(p^2,mean = 0,sd=2), ncol = p,nrow = p)
  theta = theta + t(theta)
  diag(theta)=1
  theta = theta * real_stru
  theta=MakePositiveDefinite(theta,pd_strategy = "diagonally_dominant",scale = T)$omega
  #colnames(theta)=paste0("metabolite",1:p)
  #rownames(theta)=paste0("metabolite",1:p)
  sigma=solve(theta)
  fullCovariance= lapply(cc,function(x,cc) {kronecker(cc,x)},x=sigma)
  fulldata=c()
  a=c()
  for (i in 1:n) {
    ai=c()
    m3=length(timepoint[[i]])
    mu=rep(0,m3*p)
    data=MASS::mvrnorm(1,mu=mu,Sigma = fullCovariance[[i]])
    for (j in 1:m3) {
      ai=as.data.frame(rbind(ai,data[c(((j-1)*p+1) :(j*p))]))
    }
    subject=rep(paste0("subject",i),m3)
    ai=cbind(subject,timepoint[[i]],ai)
    a=rbind(a,ai)
  }

  return(list(data=a,B=theta,timepoint=timepoint))
}



#' @title Simulate data from a given network model
#' @description
#' This function generates three types of data based on the given network type. These data can then be used to test the
#' efficiency of algorithms.
#' @param type what type of dat should be generated. There are three types of data that can be simulated through the function, which are
#' \code{general}, \code{longihomo} and \code{longiheter} respecitvely.
#' @param n number of subjects
#' @param p number of nodes
#' @param m1 number of edges
#' @param m2 number of edges between general and individual network
#' @param m3 number of correlated data within one subject
#' @param cc correlation matrix between tissues
#' @param tau parameters in correlation matrix
#' @param alpha parameter in correlation matrix
#' @importFrom fake MakePositiveDefinite
#' @importFrom MASS mvrnorm
#' @returns data frame
#' @export
Simulate=function(type=c("general","longihomo","longiheter"),n=20,p=20,m1=20,m2=10,m3=3,cc=diag(m3),tau=c(2,1),alpha=2){
  type=match.arg(type)

  if (type=="general"){
    data=simulate_general(n=n,p=p,m1=m1,m2=m2,m3=m3,cc)
  }

  if (type=="longihomo"){
    data=simulate_long(n=n,p=p,m1=m1,tau=tau)
  }

  if (type=="longiheter"){

    data=simulate_heter(n=n, p=p,m1=m1, alpha=alpha)
  }

  return(data)
}



#' Title
#'
#' @param type model type
#' @param data raw data
#' @param group group variable
#' @param lambda tuning parameter
#' @param nlam number of tuning parameter
#' @param lam.min.ratio ratio of largest lambda vs smallest lambda
#' @param K cv folds
#' @param expFix given parameter
#' @param trace whether show the process
#' @param cores parallel computing
#' @returns list
#' @import parallel foreach doParallel

cvplglasso=function(type=c("general","expFixed","twoPara"), data,group=NULL,
                     lambda=NULL,nlam=10,lam.min.ratio=0.01, K, expFix=1,trace=FALSE, cores=1){

  type=match.arg(type)


  if (is.null(group) && any( !is.vector(lambda) |  !all(is.numeric(lambda)) | !all(lambda>0)))
  {stop("group and lambda does not match!")}

  if (ncol(lambda)!=2 && !is.null(group)){stop("lambda should be a n by 2 matrix when group is specified!")}

  if (any(lambda<=0)) {stop("tuning parameter lambda should be positive!")}

  if (any(K<=1 | K%%1 !=0)){
    stop("K should be an integer greater than 1!")
  }

  num_cores=detectCores()
  if (cores > num_cores) {
    cat("\nOnly detected", paste(num_cores, "cores...", sep = " "))
  }
  if (cores > K) {
    cat("\nNumber of cores exceeds K... setting cores = K")
    cores = K
  }
  cluster = makeCluster(cores)
  registerDoParallel(cluster)
  n=length(unique(data[,1]))
  subjects=unique(data[,1])
  p=ncol(data)-2
  ind = sample(n)

  X=data[,-c(1,2)]


  S = (nrow(X) - 1)/nrow(X) * cov(X)
  # crit.cv = match.arg(crit.cv)
  # start = match.arg(start)

  Sminus = S
  diag(Sminus) = 0
  if (is.null(lambda)) {
    if (!((lam.min.ratio <= 1) && (lam.min.ratio > 0))) {
      cat("\nlam.min.ratio must be in (0, 1]... setting to 1e-2!")
      lam.min.ratio = 0.01
    }
    if (!((nlam > 0) && (nlam%%1 == 0))) {
      cat("\nnlam must be a positive integer... setting to 10!")
      nlam = 10
    }
    lam.max = max(abs(Sminus))
    lam.min = lam.min.ratio * lam.max
    lambda = 10^seq(log10(lam.min), log10(lam.max), length = nlam)
    if (!is.null(group)){
      lambda=expand.grid(lambda,lambda)
    }
  }
  else {
    if (is.null(group)){
      lambda = sort(lambda)
    }
  }


  nnlambda=ifelse(is.null(group),length(lambda),nrow(lambda))
  cv_error=matrix(0,nrow=nnlambda,ncol=K)

  k=NULL
  CV = foreach(k = 1:K, .packages = "CVglasso", .combine = "cbind",
               .inorder = FALSE) %dopar% {

                 if (trace) {
                   progress = txtProgressBar(max = K, style = 3)
                 }

                 leave.out =subjects[ind[(1 + floor((k - 1) * n/K)):floor(k *
                                                                            n/K)]]
                 indexValid=which(data[,1] %in% leave.out)
                 data.train = data[-indexValid, , drop = FALSE]
                 data_bar = apply(data.train[,-c(1,2)], 2, mean)
                 data.train[,-c(1,2)] = scale(data.train[,-c(1,2)], center = data_bar, scale = FALSE)
                 data.valid = data[indexValid,, drop = FALSE]
                 data.valid[,-c(1,2)] = scale(data.valid[,-c(1,2)], center = data_bar, scale = FALSE)
                 group.train=group[-indexValid]
                 group.valid=group[indexValid]
                 #S.train = crossprod(data.train[,-c(1,2)])/(dim(data.train)[1])
                 #S.valid = crossprod(data.valid[,-c(1,2)])/(dim(data.valid)[1])




                 if (type=="general" && is.null(group)){
                   aa= sapply(lambda, function(x) lglasso(data=data.train,lambda=x,type=type)$wi, simplify = FALSE)
                   cc=lapply(aa, function(B){
                     M=ifelse(abs(B)<=10^(-2), 0,1)
                     diag(M)=2
                     M
                   })
                   bb=lapply(cc, function(B) cvError(data.train=data.train,data.valid=data.valid,B=B))
                   bb=do.call(c,bb)
                 }



                 if (type=="general" && !is.null(group)){

                   aa=apply(lambda, 1,function(x){lglasso(data=data.train,lambda=x,type=type, group=group.train)$wiList} )


                   cc=lapply(aa, function(B){
                     lapply(B, function(Z){
                       M=ifelse(abs(Z)<=10^(-1), 0,1)
                       diag(M)=2
                       M
                     }
                     )
                   }
                   )

                   bb=lapply(cc, function(BB){
                     cvError(data.train=data.train,data.valid=data.valid,B=BB, group.valid=group.valid, group.train = group.train)
                   }
                   )
                   bb=do.call(c,bb)
                 }


                 if (type=="expFixed" && is.null(group)){

                   aa= sapply(lambda, function(x) lglasso(data=data.train,lambda=x,type=type, expFix = expFix)$wi, simplify = FALSE)
                   cc=lapply(aa, function(B){
                     M=ifelse(abs(B)<=10^(-2), 0,1)
                     diag(M)=2
                     M
                   })
                   bb=lapply(cc, function(B) cvError(data.train=data.train,data.valid=data.valid,B=B))
                   bb=do.call(c,bb)
                 }

                 if (type=="expFixed" && !is.null(group)){
                   aa= apply(lambda,1, function(x) {lglasso(data=data.train,lambda=x,type=type, expFix = expFix, group=group.train)$wiList})
                   cc=lapply(aa, function(B){
                     lapply(B, function(Z){
                       M=ifelse(abs(Z)<=10^(-1), 0,1)
                       diag(M)=2
                       M
                     }
                     )
                   }
                   )

                   bb=lapply(cc, function(BB) cvError(data.train=data.train,data.valid=data.valid,B=BB,
                                                      group.valid=group.valid, group.train = group.train))
                   bb=do.call(c,bb)
                 }

                 if (type=="twoPara" && is.null(group)){
                   aa= sapply(lambda, function(x) lglasso(data=data.train,lambda=x,type=type)$wi, simplify = FALSE)
                   cc=lapply(aa, function(B){
                     M=ifelse(abs(B)<=10^(-2), 0,1)
                     diag(M)=2
                     M
                   })
                   bb=lapply(cc, function(B) cvError(data.train=data.train,data.valid=data.valid,B=B))
                   bb=do.call(c,bb)
                 }

                 if (type=="twoPara" && !is.null(group)){
                   # browser()
                   aa= apply(lambda,1, function(x) {lglasso(data=data.train,lambda=x,type=type, group=group.train)$wiList})
                   cc=lapply(aa, function(B){
                     lapply(B, function(Z){
                       M=ifelse(abs(Z)<=10^(-1), 0,1)
                       diag(M)=2
                       M
                     }
                     )
                   }
                   )

                   bb=lapply(cc, function(BB) cvError(data.train=data.train,data.valid=data.valid,B=BB,
                                                      group.valid=group.valid, group.train = group.train))
                   bb=do.call(c,bb)
                 }
                 cv_error=bb

                 if (trace) {
                   setTxtProgressBar(progress,  k)
                 }

                 return(cv_error=cv_error)
               }

  stopCluster(cluster)
  output=structure(list(cv_error=CV, lambda=lambda),class="cvlglasso")
  return(output)
}


#' @title Cross validation for \code{lglasso}
#' @description
#' The function computes the cross validation errors for one of the three network models in \code{lglasso} command.
#' @param type underlying model type, either \code{general}, \code{longihomo} or \code{longiheter}.
#' @param data raw data
#' @param group group variable
#' @param lambda tuning parameter
#' @param nlam number of tuning parameter
#' @param lam.min.ratio ratio of largest lambda vs smallest lambda
#' @param K cv folds
#' @param expFix given parameter
#' @param trace whether show the process
#' @param cores parallel computing
#'
#' @returns list of which the first component is the cross validation errors and the second component is the corresponding
#' tuning parameters
#' @export
#' @import parallel foreach doParallel
#'
CVlglasso=function(type=c("general","expFixed","twoPara"), data,group=NULL,
                    lambda=NULL,nlam=10,lam.min.ratio=0.01, K, expFix=1,trace=FALSE, cores=NULL){
if (is.null(cores)){

  results=cvlglasso(type=type,data=data,group=group,lambda = lambda,nlam=nlam,
                    lam.min.ratio=lam.min.ratio, K=K, expFix=expFix,trace=trace)

}else{

  results=cvplglasso(type=type,data=data,group=group,lambda = lambda,nlam=nlam,
                    lam.min.ratio=lam.min.ratio, K=K, expFix=expFix,trace=trace,cores=cores)
}
  return(results)
}


