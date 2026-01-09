
#' Title
#'
#' @param t a vector specify the time points corresponding to the data
#' @param tau the damping rate parameter with length 1 or 2
#'
#' @returns a square matrix used to constrct the likelihood
#'

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
        M[i,j]= exp(-tau[1]*(abs(t[i]-t[j])))
        M[j,i]=M[i,j]
      }
    }
    diag(M)=1
  }
  return(M)
}




#' Title
#'
#' @param B a p by p given precision matrix
#' @param data a (p+2)-by-n data frame
#' @param type specify how to model the covariance matrix
#' @param expFix the parameter in variance function when the data are longitudinal
#' @param maxit the maximum iteration number
#' @param tol the minimum difference for the algorithm to be called converged
#' @param lower vector of length 1 or 2 which specifies the lower bounds for alpha_1 (and alpha_2) in the correlation matrix
#' @param upper vector of length 1 or 2 which specifies the upper bounds for alpha_1 (and alpha_2) in the correlation matrix
#' @param ... other unspecifed parameterss
#' @returns a list of matrices

AA=function(B,data,type=c("general","expFixed"),expFix,maxit=30,
            tol=10^(-4),lower=c(0.01,0.1),upper=c(10,5),...){
  ### clustered data
if (!is.list(B)){
  B=list(B)
}

if (!is.list(data)){
  data=as.list(data)
}

  if (length(B)!= length(data)){
    stop(" B should have same length as  data!")
  }

  type=match.arg(type)
  if (type=="general"){
    corMatrix=vector("list",length(B))
    for (i in 1:length(B)) {
      dd=data[[i]]
      bb=B[[i]]
    m3=length(unique(dd[,2]))
    p=ncol(dd)-2
    nn=length(unique(dd[,1]))
    #subjects=unique(dd[,1])
    #browser()
    if (ncol(bb)!=p | nrow(bb)!=p){
      stop("Inputs do not match with each other!")
      }
    A=Variable(m3,m3,PSD=TRUE) # tissue wise inverse correlation matrix
    obj1=(p*nn)/2*log_det(A)
    data_sub=split(dd[,-c(1,2)],f=factor(dd[,1],levels=unique(dd[,1])))
    amatrix=Reduce("+",lapply(data_sub, function(B,xx){as.matrix(xx)%*%B%*%t(as.matrix(xx))}, B=bb))
    obj2=-0.5*matrix_trace(A%*%amatrix)
    obj=-(obj1+obj2)
    constr=list(CVXR::diag(A)==1)
    prob=Problem(Minimize(obj),constr)
    results=CVXR::solve(prob)
    corMatrix[[i]]=solve(results$getValue(A))
  }
    return(list(corMatrix=corMatrix))
}
  ### longitudinal data

  if (type == "expFixed"){
    corMatrix=vector("list",length(B))
    Tau=c()
    for (i in 1:length(B)) {
      dd=data[[i]]
      bb=B[[i]]
      p=ncol(dd)-2
      nn=length(unique(dd[,1]))
      #subjects=unique(dd[,1])
      #browser()
      if (ncol(bb)!=p | nrow(bb)!=p){
        stop("Inputs do not match with each other!")
      }

    #YY=Variable(p,p,PSD=TRUE)
    time_list=split(dd[,2],f=factor(dd[,1],levels=unique(dd[,1])))
    subdata_list=split(dd[,-c(1,2)],factor(dd[,1],levels=unique(dd[,1])))

    likefun=function(tau){
      obj1=0
      obj2=0
      for (i in 1:nn) {
        datai=as.matrix(subdata_list[[i]])
        t=time_list[[i]]
        Aa=phifunction(t=t,tau = tau)
        amatrix=datai%*%bb%*%t(datai)
        obji1=-0.5*p*log(det(Aa))
        obji2= -0.5*sum(diag(solve(Aa)%*%amatrix))
        obj1=obj1+obji1
        obj2=obj2+obji2
      }
      obj=-(obj1+obj2)
    }

    tau=stats::optim(c(tau0[1]),likefun,method = "L-BFGS-B",lower = lower,upper = upper)$par
    Tau=c(Tau,tau)
    A=lapply(time_list,phifunction,tau=tau)
    corMatrix[[i]]=A
  }
    return(list(corMatrix=corMatrix,tau=Tau))
}
  ### longitudinal data
  # if (type == "twoPara"){
  #   p=ncol(data)-2
  #   nn=length(unique(data[,1]))
  #   subjects=unique(data[,1])
  #   if (ncol(B)!=p | nrow(B)!=p){stop("Inputs do not match with each other!")}
  #   #YY=Variable(p,p,PSD=TRUE)
  #   time_list=split(data[,2],data[,1])
  #   subdata_list=split(data[,-c(1,2)],data[,1])
  #
  #   likefun=function(tau){
  #     obj1=0
  #     obj2=0
  #     for (i in 1:nn) {
  #       datai=as.matrix(subdata_list[[i]])
  #       t=time_list[[i]]
  #       Aa=phifunction(t=t,tau = c(tau,1))
  #       amatrix=datai%*%B%*%t(datai)
  #       obji1=-0.5*p*log(det(Aa))
  #       obji2= -0.5*sum(diag(solve(Aa)%*%amatrix))
  #       obj1=obj1+obji1
  #       obj2=obj2+obji2
  #     }
  #
  #     obj=-(obj1+obj2)
  #   }
  #
  #   tau=stats::optim(c(tau0),likefun,method = "L-BFGS-B",lower = lower,upper = upper)$par
  #   A=lapply(time_list,phifunction,tau=tau)
  #   return(list(corMatrixList=A,tau=tau))
  # }

}









#' Title
#'
#' @param A
#' @param data
#' @param lambda
#' @param type
#' @param diagonal
#' @param maxit
#' @param tol
#' @param lower
#' @param upper
#' @param start
#' @param w.init
#' @param wi.init
#' @param ...
#'
#' @returns

BB=function(A,data,lambda,type=c("general","expFixed"),diagonal=TRUE,maxit=100,
            tol=10^(-4),lower=c(0.01,0.1),upper=c(10,5), start=c("warm","cold"),w.init=NULL,wi.init=NULL,...){

  type=match.arg(type)
  start=match.arg(start)



  if (!is.list(data) | !is.list(A)){
    stop("A and data must be a list!")
  }

  if (length(A)!= length(data)){
    stop(" A should have same length as  data!")
  }






  if (type=="general"){
    for (i in 1:length(A)) {
      Ai=A[[i]]
      datai=data[[i]]
      if (nrow(Ai)!=length(unique(datai[,2]))){
        stop("The format of A does not match the format of data!")
      }
    }


    m= length(A)
    B=vector("list",m)
    obj=0
    aa=0
    bb=0
     for (i in 1:m){
      Ai=A[[i]]
      dd=data[[i]]
      p=ncol(dd)-2
      data_sub=split(dd[,-c(1,2)],f=factor(dd[,1],levels=unique(dd[,1])))
      nn=length(unique(dd[,1]))
      amatrix=Reduce("+",lapply(data_sub,
      function(A,xx){t(as.matrix(xx))%*%solve(A)%*%as.matrix(xx)}, A=Ai))
          B[[i]]=Variable(p,p,PSD=TRUE) # tissue wise inverse correlation matrix
          obj=(nn*nrow(Ai))/2*log_det(B[[i]])-0.5*matrix_trace(B[[i]]%*%amatrix)+obj
          # Create a mask matrix
          mask <- matrix(1, p, p)
          diag(mask) <- 0
          aa=aa+ sum(abs(B[[i]]*mask))

          if (m>1){
          if (i>=2){
            for (j in 1:i) {
              bb=bb+sum(abs(B[[i]]-B[[j]])*mask)
            }
          }
          }
     }
#     if (m==1){
#       contr=list(aa<=1/lambda[1])
#     }else{
#       contr=list(aa<=1/lambda[1],bb<=1/lambda[2])
# }

    if (m==1){
      S_est=list(glasso::glasso(s=amatrix/(nn*nrow(Ai)), rho = lambda[1])$wi)
    }else{
obj=-obj+lambda[1]*aa+lambda[2]*bb
#obj=-obj+aa+bb
prob=Problem(Minimize(obj))
result=CVXR::solve(prob)
S_est= lapply(B, function(x) result$getValue(x))
}
return(wiList=S_est)
  }



  if (type == "expFixed"){
    m= length(A)
    for (i in 1:length(A)) {
      Ai=A[[i]]
      datai=data[[i]]
      if (length(Ai)!=length(unique(datai[,1]))){
        stop("The format of A does not match the format of data!")
      }
      subjects=unique(datai[,1])
      for (j in 1:length(Ai)) {
        Aij=Ai[[j]]
        index=which(datai[,1]==subjects[j])
        dataij=datai[index,]
        if (nrow(Aij)!= nrow(dataij))
        {stop("The format of A does not match the format of data!")}
      }

    }



    B=vector("list",m)
    obj=0
    aa=0
    bb=0
    for (i in 1:m){
      Ai=A[[i]]
      dd=data[[i]]
      nn=length(unique(dd[,1]))
      p=ncol(dd)-2
      data_sub=split(dd[,-c(1,2)],factor(dd[,1],unique(dd[,1])))
      B[[i]]=Variable(p,p,PSD=TRUE) # tissue wise inverse correlation matrix
      if (length(Ai)!=length(data_sub)){stop("Data do not match!")}
      amatrix=0
      # Create a mask matrix
      mask1 <- matrix(lambda[1], p, p)
      diag(mask1) <- 0
      mask2 <- matrix(lambda[2], p, p)
      diag(mask2) <- 0
      for (j in 1:nn) {
        xx=as.matrix(data_sub[[j]])
        yy=solve(as.matrix(Ai[[j]]))
        amatrix=t(xx)%*%yy%*%xx+amatrix
        #amatrix=t(xx)%*%xx+amatrix
      }
      obj=log_det(B[[i]])-matrix_trace(B[[i]]%*%amatrix)/nrow(dd)+obj
      aa=aa+ sum(abs(B[[i]])*mask1)

      if (m>1){
        if (i>=2){
          for (j in 1:(i-1)) {
            bb=bb+sum(abs((B[[i]]-B[[j]]))*mask2)
          }
        }
      }
    }

    obj=-obj+aa+bb
if (m==1){
  S_est=list(glasso::glasso(s=amatrix/(nrow(dd)),rho=lambda[1])$wi)
}else{
    prob=Problem(Minimize(obj))
    result=CVXR::solve(prob)
    S_est= lapply(B, function(x) result$getValue(x))
    return(wiList=S_est)
  }
  }
}
  # if (type == "twoPara"){
  #   p=ncol(data)-2
  #   nn=length(unique(data[,1]))
  #   subjects=unique(data[,1])
  #
  #   data_sub=split(data[,-c(1,2)],data[,1])
  #   if (length(A)!=length(data_sub)){stop("Data do not match!")}
  #   bmate=0
  #   for (i in 1:length(A)) {
  #     xx=as.matrix(data_sub[[i]])
  #     yy=solve(as.matrix(A[[i]]))
  #     bmate=t(xx)%*%yy%*%xx+bmate
  #   }
  #
  #   if (start=="warm"){
  #    # w.init=cov(data[,-c(1,2)])
  #   #  wi.init=diag(ncol(w.init))
  #     bb=glasso::glasso(s=bmate/nrow(data),rho=lambda[1], penalize.diagonal = diagonal, start = start,w.init = w.init,wi.init = wi.init)
  #   }else{
  #     bb=glasso::glasso(s=bmate/nrow(data),rho=lambda[1], penalize.diagonal = diagonal, start = "cold")
  #   }
  #
  #
  #   return(list(preMatrix=bb$wi,corMatrix=bb$w))
  #
  # }




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
#' which are referred as \code{general}, \code{expFixed} respectively.
#'  Please see the details in the below for the meaning of each.
#' @param expFix  numeric number used in the model specification
#' @param group  vector  of length \code{n} if supplied which specify which data
#'  points need to be grouped together to infer the heterogeneous networks for, e.g, pre/post vaccination.
#' @param maxit the maximum iterations for the estimation.
#' @param tol the minimum value for  convergence criterion
#' @param lower  vector of length 1 or 2 which specifies the lower bounds for alpha_1 (and alpha_2) in the correlation matrix
#' @param upper  vector of length 1 or 2 which specifies the upper bounds for alpha_1 (and alpha_2) in the correlation matrix
#' @param start how to start the initial values for lglasso algorithm
#' @param w.init initial value for covariance matrix
#' @param wi.init inital value for precision matrix
#' @param trace whether or not show the progress of the computation
#' @param ... other inputs
#' @import glasso CVXR
#' @importFrom "utils" "fix"
#' @export
#' @return list which include following components:
#'
#' \code{w} the general covariance matrix estimate;
#'
#' \code{wList} list representing the individual covariance matrix estimate;
#'
#' \code{wi} the general precision matrix estimate;
#'
#' \code{wiList} list representing the individual precision matrix estimate;
#'
#' \code{v} the correlation matrix between specified classes;
#'
#' \code{vList} list representing the individual correlation matrix;
#'
#' \code{tauhat} the correlation parameters for longitudinal data
#'
#' @details This function implements three statistical models for  network inference,
#' according to how the correlations is specified between time points (or tissues or
#' contents in some clinical studies). These three models are referred as
#'  \code{general}, \code{expFixed}.Let's say we have
#'  two time points,t_i,t_j, then in model \code{general}, the correlation is
#'   tau_ij, while in model *expFixed*, we have  tau=exp(-alpha_1|t_1-t_2|^(-alpha_2))
#'   with alpha_2 need to be pre-specified (default is alpha_2=1). In model \code{twoPara},
#'   both alpha_1 and alpha_2 is unknown and need to be inferred from the data.
#'    For longitudinal data, model \code{expFixed} is recommended while for omics data
#'    from different tissues or contents, model \code{general} should be adopted.
#'
#'
lglasso=function(data,lambda,group=NULL,random=FALSE,expFix=1,N=100,maxit=30,
                 tol=10^(-3),lower=c(0.01,0.1),upper=c(10,5), start=c("cold","warm"),
                 w.init=NULL, wi.init=NULL,trace=FALSE, type=c("expFixed"),
                 ...)

  {
if (type!="expFixed"){
  stop("type can only be expFixed currently!")
}

  p=ncol(data)-2
  X_bar = apply(data[,-c(1,2)], 2, mean)
  data[,-c(1,2)] = scale(data[,-c(1,2)], center = X_bar, scale = FALSE)


  if (random==FALSE){

    if (!is.null(group))  {

      if (length(group)!=nrow(data)){
        stop("group should be the same length of the columns of data!")
      }

      data=split(data,f=factor(group,levels = unique(group)))

      if (!length(lambda)==2){
        stop("Arguments (group, lambda) do not match!")
      }

    }

if (is.null(group))  {
  group=rep(1,nrow(data))
  data=list(data)
 if (length(lambda)!=1){
  stop("Arguments (group, lambda) do not match!")
 }
}





  if (!all(lambda>0)){
    stop("lambda must be positive!")
  }

  type=match.arg(type)
  start=match.arg(start)
  # Create a mask matrix
  mask <- matrix(1, p, p)
  diag(mask) <- 0




  if (type == "expFixed"){
    if (is.null(expFix) | length(expFix)!=1 | !is.numeric(expFix)){
      stop("Argument expFix is not correctly specified!")
    }
    A=vector("list",length(data))
    B=vector("list",length((data)))

    for (i in 1:length(A)) {
      dd=data[[i]]
      subjects=unique(dd[,1])
      A[[i]]=vector("list",length(subjects))
      for (j in 1:length(A[[i]])) {
        index=which(dd[,1]==subjects[j])
        A[[i]][[j]]=diag(length(index))
      }
      B[[i]]=diag(p)
    }
k=0
while(1){
  k=k+1
A1=AA(data = data,B = B, type=type,expFix,...)
B1=BB(data=data,A=A,lambda = lambda, type=type,start=start, w.init=w.init,wi.init=wi.init, ...)
d1=c()
d2=c()
for (i in 1:length(B1)) {
  d1=c(d1,round(max(abs(B[[i]]-B1[[i]])*mask),3))
  d2=c(d2,round(max(abs(tau0-A1$tau)),3))
}
if (trace){
  print(paste0("iteration ",k, " precision difference: ",max(d1) , " /correlation tau difference: ",max(d2)))
}

if (max(d1)<=tol && max(d2)<= tol ){
    output=structure(list(wi=B1, v=A1$corMatrix, tau=A1$tau), class="lglasso")
  break
}else{
  A=A1$corMatrix
  B=B1
  tau0=A1$tau
}


if (k>=maxit){
  warning("Algorithm did not converge!")
  output=structure(list(wi=B1, v=A1$corMatrix, tau=tau0), class="lglasso")
  break
}

}

    }
return(output)
  }


  if (random==TRUE){

    if (is.null(group))  {
      if (length(lambda)!=1){
        stop("Arguments (group, lambda) do not match!")
      }
    }

    if (!is.null(group))  {

      if (length(group)!=nrow(data)){
        stop("group should be the same length of the columns of data!")
      }


      if (!length(lambda)==2){
        stop("Arguments (group, lambda) do not match!")
      }


    }



    if (!all(lambda>0)){
      stop("lambda must be positive!")
    }

output=lglassoHeter(data=data,lambda=lambda,expFix=expFix,N=N,group=group,maxit=maxit,
                   tol=tol,trace=trace)
  }
  return(output)
}

#' Title
#'
#' @param tau
#' @param datai
#' @param wi
#' @param alpha
#' @param groupi
#'
#' @returns
#' @export
#'
#' @examples
conDensityTau=function(tau,expFixed, datai,wi,alpha,groupi){



    if (length(groupi)!=nrow(datai)){
      stop("group should be the same length of the columns of data!")
    }

    if (length(unique(groupi))!= length(wi)){
      stop("the format of groupi does not match that of sigmaM!")
    }

  data=split(datai,f=factor(groupi,levels = unique(groupi)))
  nn=length(unique(groupi))
  likelihood=1
  p=nrow(wi[[1]])
  for (i in 1:nn) {
    timepoints=data[[i]][,2]
    phiMi=phifunction(t=timepoints,tau=tau)
    s=t(as.matrix(data[[i]][,-c(1,2)]))%*%solve(phiMi)%*%as.matrix(data[[i]][,-c(1,2)])%*%wi[[i]]
    a=det(phiMi)^(-p/2)*exp(-0.5*sum(diag(s)))
    likelihood=a*likelihood
  }
  likelihoodi=likelihood*exp(-alpha*tau)
  return(likelihoodi=likelihoodi)
}

#' Title
#'
#' @param n the number of random samples
#' @param datai the data for subject i
#' @param wi the given precision matrix
#' @param alpha the exponential distribution with rate alpha
#' @param groupi specify how datai is grouped
#'
#' @returns a data frame for samples and their weights

importanceSample=function(n,datai,wi,alpha,groupi){
  dd1=matrix(rexp(n=n,rate=alpha),ncol=1)
  likelihood1=apply(dd1, 1, conDensityTau,datai=datai,wi=wi,alpha=alpha,groupi=groupi)
  likelihood2=apply(dd1, 1, dexp,rate=alpha)
  index1=which(!is.nan(likelihood1))
  index2=which(!is.infinite(likelihood1) )
index=intersect(index1,index2)
if (length(index)==0){stop("No valid samples are generated!")}
  weights=likelihood1[index]/likelihood2[index]
  norWeights=weights/sum(weights)
  aa=data.frame(sample=dd1[index],weight=norWeights)
  return(aa)
}

#' Title
#'
#' @param importancesSample given random samples
#' @param datai the data for subject i
#' @param groupi specify how datai is grouped
#' @returns a list of estimated

importanceEstimates=function(importancesSample,datai,groupi){

  data=split(datai,f=factor(groupi,levels = unique(groupi)))
  sample0=importancesSample[,1]
  ff=seq(1,length(sample0),1)
  sampleTau=split(sample0,f=ff)
  weightTau=importancesSample[,2]
  estimateTau=sum(sample0*weightTau)
  estimatePhi=vector("list",length(data))
  for (i in 1:length(data)) {
    timepoints=data[[i]][,2]
    phiMi=lapply(sampleTau,phifunction,t=timepoints)
    #s=lapply(phiMi,function(phiMi) {t(data[[i]][,-c(1,2)])%*%solve(phiMi)%*%data[[i]][,-c(1,2)]})
    aa=0
    for (j in 1:length(phiMi)) {
      aa=aa+phiMi[[j]]*weightTau[j]
    }
    estimatePhi[[i]]=aa
  }
  estimates=list(estimateTau=estimateTau,estimatePhi=estimatePhi)
  return(estimates)
}


#' Title
#'
#' @param data longitudinal data set
#' @param wi given precision matrix
#' @param alpha expoential distribution with rate alpha
#' @param group specify how data is grouped
#' @param l number of random samples in importance sampling
#' @returns a list for estimates of tau and AA
AAheter=function(data,wi,alpha,group,l){
  subjects=unique(data[,1])
  nn=length(unique(group))
  A=vector("list",length(subjects))
  Tau=matrix(nrow=length(subjects),ncol=1)
  simTau=matrix(rexp(n=l,rate=alpha),ncol=1)
  #browser()
  dataList=split(data,f=factor(data[,1],levels=subjects))
  groupList=split(group,f=factor(data[,1],levels=subjects))
  for (i in 1:length(subjects)) {
    datai=dataList[[i]]
    groupi=groupList[[i]]
      imSample=importanceSample(n=l,datai=datai,wi=wi,alpha=alpha,groupi =groupi )
      index=which(!is.nan(imSample[,2]))
      imSample=imSample[index,]
      imporResults=importanceEstimates(importancesSample=imSample,datai=datai,groupi = groupi)
      Tau[i,1]=imporResults$estimateTau
      A[[i]]=imporResults$estimatePhi
  }
  AA=vector("list",nn)

  for (i in 1:nn) {
    AA[[i]]=vector("list",length(subjects))
    for (j in 1:length(subjects)) {
      AA[[i]][[j]]=A[[j]][[i]]
    }
  }

  return(list(Tau=Tau,AA=AA))
}



#' Title
#'
#' @param data
#' @param lambda
#' @param expFix
#' @param group
#' @param maxit
#' @param tol
#' @param trace
#'
#' @returns
#' @export
#'
#' @examples
lglassoHeter=function(data,lambda,expFix,group,maxit,
                      tol=10^(-3),trace=FALSE,start=c("warm","cold"), w.init=NULL, wi.init=NULL, N)

{
  p=ncol(data)-2
  m=length(unique(data[,1]))
  if (is.null(group))  {
    group=rep(1,nrow(data))
    if (length(lambda)!=1){
      stop("Arguments (group, lambda) do not match!")
    }
  }

  if (!is.null(group))  {

    if (length(group)!=nrow(data)){
      stop("group should be the same length of the columns of data!")
    }


    dataList=split(data,f=factor(group,levels = unique(group)))

    if (length(lambda)!= length(unique(group))){
      stop("Arguments (group, lambda) do not match!")
    }

  }

nn=length(unique(group))

  if (!all(lambda>0)){
    stop("lambda must be positive!")
  }

  # Create a mask matrix
  mask <- matrix(1, p, p)
  diag(mask) <- 0


  if (is.null(expFix) | length(expFix)!=1 | !is.numeric(expFix)){
    stop("Argument expFix is not correctly specified!")
  }

  B=vector("list",nn)

  for (i in 1:length(B)) {
    B[[i]]=diag(p)
  }
  k=0
  tau0=rep(1,m)
  alpha0=1
  while(1){
    k=k+1
    print(k)
    print(paste0("alpha estimate: ", alpha0))
    #if (k==11){browser()}
    A1=AAheter(data=data,wi=B,alpha=alpha0,group=group,l=N)
    B1=BB(data=dataList,A=A1$AA,lambda = lambda, type="expFixed",start=start, w.init=w.init,wi.init=wi.init)
    d1=c()
    d2=c()
    for (i in 1:length(B1)) {
      d1=c(d1,round(max(abs(B[[i]]-B1[[i]])*mask),3))
      d2=c(d2,round(max(abs(alpha0-1/mean(A1$Tau))),3))
    }
    if (trace){
      print(paste0("iteration ",k, " precision difference: ",max(d1) , " /correlation alpha difference: ",max(d2)))
    }

    if (max(d1)<=tol && max(d2)<= tol ){
      output=structure(list(wi=B1, tau=A1$Tau,alpha=1/mean(A1$Tau)), class="lglasso")
      break
    }else{
      A=A1$AA
      B=B1
      tau0=A1$Tau
      alpha0=1/mean(tau0)
    }

    if (k>=maxit){
      warning("Algorithm did not converge!")
      output=structure(list(wi=B1, tau=tau0,alpha=alpha0), class="lglasso")
      break
    }

  }


  return(output)
}



#   if (type == "twoPara"){
#     p=ncol(data)-2
#     fre=table(data[,1])
#     A=vector("list",length(fre))
#     for (i in 1:length(A)) {
#       A[[i]]=diag(fre[i])
#     }
#     B=diag(p)
#     k=0
# tau0=tau0
#     while(1){
#       k=k+1
#       A1=AA(data = data,B = B, type=type, fix=fix)
#       B1=BB(data=data,A=A,lambda = lambda[1], type=type,start=start, w.init=w.init,wi.init=wi.init,...)
#
#       d1=round(max(abs(B-B1$preMatrix)),4)
#       d2=round(max(abs(tau0-A1$tau)),4)
#       if (trace){
#       print(paste0("iteration ",k, " precision difference: ",d1 , " /correlation tau difference: ",d2))
#       }
#       if (d1<= tol & d2<= tol ){
#
#         if (is.null(group)){
#           output=structure(list(wi=B1$preMatrix, wList=NULL,vList=A1$corMatrixList, tauhat=tau0), class="lglasso")
#         }else{
#
#           individial=heternetwork(data = data,lambda = lambda[2],homo = B1$preMatrix, group=group)
#           output=structure(list(wi=B1$preMatrix, wiList=individial, vList=A1$corMatrixList, tauhat=tau0), class="lglasso")
#         }
#         break
#       }else{
#         A=A1$corMatrixList
#         B=B1$preMatrix
#         tau0=A1$tau
#       }
#       if (k>= maxit){
#         warning("Algorithm did not converge!")
#
#         if (is.null(group)){
#           output=structure(list(wi=B1$preMatrix, wList=NULL,vList=A1$corMatrixList, tauhat=tau0), class="lglasso")
#
#         }else{
#           #browser()
#           individial=heternetwork(data = data,lambda = lambda[2],homo = B1$preMatrix, group=group)
#           output=structure(list(wi=B1$preMatrix, wiList=individial, vList=A1$corMatrixList, tauhat=tau0), class="lglasso")
#       }
#         break
#           }
#
#     }
#
#   }












cvErrorji=function(data.train,data.valid,bi){

  if (any(! bi %in% c(0,1,2))) {stop("entries of vector bi should be 0,  1 or 2!")}
  i=which(bi==2)
  cv_error=c()


    index=which(bi==1)
    if (length(index)==0){
      cv_error=stats::var(data.valid[,i+2])

    }else{
      if (length(index)>= nrow(data.train)){
        print(paste("number of variable ", length(index)))
        stop("network is too dense for model training!")
      }

      y=data.train[,i+2]
      x=as.matrix(data.train[,index+2])
      #browser()
      coef.train=stats::lm(y~x)$coef
      yy=data.valid[,i+2, drop=FALSE]
      # if(nrow(data.valid)>1 && ncol(data.valid>1)){
      #   xx=cbind(1,X.valid[,index])
      # }else{
      #   xx=c(1,c(X.valid[index]))}

      xx=as.matrix(cbind(1,data.valid[,index+2,drop=FALSE]))
      err=(yy-xx%*%coef.train)^2
      cv_error=c(cv_error,mean(err[,,drop=TRUE]))
      if (any(is.na(cv_error))) {
        print("cv error is missing!")
        #browser()}
    }

  return(mean(cv_error))
}

}

cvErrorj=function(data.train,data.valid,B){
    a= mean(apply(B, 2, function(bi) cvErrorji(data.train=data.train,data.valid=data.valid,bi=bi)))
}




#' Title
#'
#' @param data.train trainng data
#' @param data.valid testing data
#' @param B given network (or network list)
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

  data.train.sub=split(data.train,f=factor(group.train,levels=unique(group.train)))
  data.valid.sub=split(data.valid,f=factor(group.valid,levels=unique(group.valid)))


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



cvlglasso=function(type=c("expFixed"), data,group=NULL,
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




  S = (nrow(X) - 1)/nrow(X) * stats::cov(X)
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
    progress = utils::txtProgressBar(max = K, style = 3)
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











  if (type=="expFixed" && is.null(group)){

   aa= sapply(lambda, function(x){lglasso(data=data.train,lambda=x,type=type, expFix = expFix)$wi}, simplify = FALSE)
   cc=lapply(aa, function(B){
     M=ifelse(abs(B)<=10^(-2), 0,1)
     diag(M)=2
     M
   })
   bb=lapply(cc, function(B){cvError(data.train=data.train,data.valid=data.valid,B=B)})
   bb=do.call(c,bb)
  }

  if (type=="expFixed" && !is.null(group)){
   aa= apply(lambda,1, function(x){lglasso(data=data.train,lambda=x,type=type, expFix = expFix, group=group.train)$wiList})
   cc=lapply(aa, function(B){
     lapply(B, function(Z){
       M=ifelse(abs(Z)<=10^(-1), 0,1)
       diag(M)=2
       M
     }
     )
   }
   )

   bb=lapply(cc, function(BB){cvError(data.train=data.train,data.valid=data.valid,B=BB,
                                      group.valid =   group.valid,group.train = group.train)})
   bb=do.call(c,bb)
  }

cv_error[,k]=bb
if (trace) {
  utils::setTxtProgressBar(progress,  k)
}
  }

output=structure(list(cv_error=cv_error, lambda=lambda),class="cvlglasso")
return(output)
}




#' Title
#'
#' @param x CVlglasso object
#' @param xvar character which specify the x axis of the plot
#' @param ... other plot arguments
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
for (i in 1:length(a1)){
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

cvplglasso=function(type=c("expFixed"), data,group=NULL,
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
  subjects=unique(data[,1])
  n=length(subjects)
  p=ncol(data)-2
  ind = sample(n)

  X=data[,-c(1,2)]


  S = (nrow(X) - 1)/nrow(X) * stats::cov(X)
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
     CV = foreach(k = 1:K, .packages = "lglasso", .combine = "cbind",
                  .inorder = FALSE) %dopar% {

                 if (trace) {
                   progress = utils::txtProgressBar(max = K, style = 3)
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


                 if (is.null(group)){

                   aa= sapply(lambda, function(x) {lglasso(data=data.train,lambda=x)$wi})
                   cc=lapply(aa, function(B){
                     M=ifelse(abs(B)<=10^(-2), 0,1)
                     diag(M)=2
                     M
                   })
                   bb=lapply(cc, function(B) {cvError(data.train=data.train,data.valid=data.valid,B=B)})
                   bb=do.call(c,bb)
                 }
                 if (!is.null(group)){
                   aa= apply(lambda,1, function(x) {lglasso(data=data.train,lambda=x, expFix = expFix, group=group.train)$wiList})
                   cc=lapply(aa, function(B){
                     lapply(B, function(Z){
                       M=ifelse(abs(Z)<=10^(-1), 0,1)
                       diag(M)=2
                       M
                     }
                     )
                   }
                   )

                   bb=lapply(cc, function(BB){cvError(data.train=data.train,data.valid=data.valid,B=BB,
                                                      group.valid=group.valid, group.train = group.train)})
                   bb=do.call(c,bb)
                 }

                 cv_error=bb

                 if (trace) {
                   utils::setTxtProgressBar(progress,  k)
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
CVlglasso=function(type=c("expFixed"), data,group=NULL,
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


