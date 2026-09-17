
#' Title
#'
#' @param t a vector specify the time points corresponding to the data
#' @param tau the damping rate parameter with length 1 or 2
#' @param expFix a scalar specifying the form of correlation function
#'
#' @returns a square matrix used to construct the likelihood
#'

phifunction=function(t,tau,expFix=1){
  n=length(t)
  if (length(tau)>1){stop("Tau should be a scalar!")}
  if (tau<=0){
    stop("tau should be positive!")
    }
  if (n==1) return(matrix(1,1,1))
   d=(abs(outer(t,t,"-")))^{expFix}
   M=exp(-tau*d)
  diag(M)=1
   M
}




#' Find the tau's for homogeneous model
#'
#' @param B a list of length 1 or 2.  If B is of length 1, then its entry is a
#'   p by p given precision matrix representing the whole data points.
#'   If B is of length 2, then they represent the pre- and post-treatment network.
#' @param data a (p+2)-by-n data frame
#' @param expFix the parameter in variance function when the data are longitudinal
#' @param maxit the maximum of iteration number
#' @param tol the minimum difference of algorithm convergence
#' @param lower vector of length 1 or 2 which specifies the lower bounds for alpha_1 (and alpha_2) in the correlation matrix
#' @param upper vector of length 1 or 2 which specifies the upper bounds for alpha_1 (and alpha_2) in the correlation matrix
#' @returns a list of matrices

AA=function(B,data,expFix=1,maxit=30,
            tol=10^(-4),lower=c(0.01,0.1),upper=c(10,5)){
  ### clustered data
if (!is.list(B)){
  B=list(B)
}

if (is.data.frame(data)){
  data=list(data)
}

  if (length(B)!= length(data)){
    stop(" B should have same length as  data!")
  }

  ### longitudinal data
    corMatrix=vector("list",length(B))
    Tau=c()
    for (i in 1:length(B)){
      dd=data[[i]]
      bb=B[[i]]
      p=ncol(dd)-2
      nn=length(unique(dd[,1]))
      if (ncol(bb)!=p | nrow(bb)!=p){
        stop("Precision matrix has to be a square one!")
      }
    time_list=split(dd[,2],f=factor(dd[,1],levels=unique(dd[,1])))
    subdata_list=split(dd[,-c(1,2)],factor(dd[,1],levels=unique(dd[,1])))

    likefun=function(tau){
      obj1=0
      obj2=0
      for (i in 1:nn) {
        datai=as.matrix(subdata_list[[i]])
        t=time_list[[i]]
        if (nrow(datai)==0 | length(t)==0){next}
        Aa=phifunction(t=t,tau = tau)
        amatrix=datai%*%bb%*%t(datai)
        obji1=-0.5*p*log(det(Aa))
        obji2= -0.5*sum(diag(solve(Aa)%*%amatrix))
        obj1=obj1+obji1
        obj2=obj2+obji2
      }
      obj=-(obj1+obj2)
    }
    tau=stats::optim(c(1),likefun,method = "L-BFGS-B",lower = lower,upper = upper)$par
    Tau=c(Tau,tau)
    A=lapply(time_list,phifunction,tau=tau)
    corMatrix[[i]]=A
  }
    return(list(corMatrix=corMatrix,tau=Tau))

}

#' Title
#'
#' @param A a list of length 1 or 2 corresponding to the number of stages. The each entry of A is a
#'  list representing all phi matrices before or after the treatment.
#' @param data a list of (p+2)-by-ni data frame
#' @param lambda given tuning parameter(s)
#' @param random logical variable indicating the type of the model
#' @param tau scalar if *random* is FALSE and a vector if *random* is TRUE
#' @returns a list with the same length as A

BB=function(A,data,lambda,random=FALSE,tau){
  if (!is.list(data) | !is.list(A)){
    stop("A and data must be lists!")
  }
  featureNames=colnames(data[[1]])[-c(1,2)]


  if (length(A)!= length(data)){
    stop(" List A should have same length as list  data!")
  }

  # if (length(A)!=length(tau)){
  #   stop(" List A should have same length as tau!")
  # }
  glev=names(A)

    m= length(A)
    for (i in 1:m){
      Ai=A[[i]]
      datai=data[[i]]
      # if (length(Ai)!=length(unique(datai[,1]))){
      #   stop("The format of A does not match the format of data!")
      # }
      subjects=unique(datai[,1])
      for (j in 1:length(subjects)) {
        Aij=Ai[[subjects[j]]]
        if (is.null(Aij)){next}
        index=which(datai[,1]==subjects[j])
        dataij=datai[index,]
        #browser()
        if (nrow(Aij)!= nrow(dataij))
        {stop("The format of A does not match the format of data!!")}
      }

    }



    B=vector("list",m)
    names(B)=glev
    likeli=0
    aa=0
    bb=0
    amatrix=list()
    sz=c()
    for (i in 1:m){
      Ai=A[[i]]
      dd=data[[i]]
      sz=c(sz,nrow(dd))
      nn=length(unique(dd[,1]))
      p=ncol(dd)-2
      data_sub=split(dd[,-c(1,2)],factor(dd[,1],unique(dd[,1])))
      B[[i]]=Variable(p,p,PSD=TRUE)
      if (length(Ai)!=length(data_sub)){stop("Data do not match!")}
      amatrix[[i]]=0
      # Create a mask matrix
      mask1 <- matrix(lambda[1], p, p)
      diag(mask1) <- 0
      mask2 <- matrix(lambda[2], p, p)
      diag(mask2) <- 0
      extra=0
      for (j in 1:nn) {
        xx=as.matrix(data_sub[[j]])
        yy=solve(as.matrix(Ai[[j]]))
        if (is.null(yy)){next}
        amatrix[[i]]=t(xx)%*%yy%*%xx+amatrix[[i]]
        extra=extra+p*log(det(as.matrix(Ai[[j]])))
      }
      if (!random){
      likeli=-extra/nrow(dd)+log_det(B[[i]])-matrix_trace(B[[i]]%*%amatrix[[i]])/nrow(dd)+likeli
    }else{
      likeli=-extra/nrow(dd)+log_det(B[[i]])-matrix_trace(B[[i]]%*%amatrix[[i]])/nrow(dd)
      -2*nn*log(mean(tau))/nrow(dd)+likeli
    }
      aa=aa+ sum(abs(B[[i]])*mask1)

      if (m>1){
        if (i>=2){
          for (j in 1:(i-1)) {
            bb=bb+sum(abs((B[[i]]-B[[j]]))*mask2)
          }
        }
      }
    }

    obj=-likeli+aa+bb
if (m==1){
  results=glasso::glasso(s=amatrix[[1]]/(nrow(dd)),rho=lambda[1])
  S_est=list(results$wi)
  colnames(S_est[[1]])=featureNames
  rownames(S_est[[1]])=featureNames
  likelihood=ifelse(!random,
                    -extra/nrow(dd)+log(det(results$wi))-sum(diag(results$wi%*%amatrix[[1]]))/(nrow(dd)),
    -extra/nrow(dd)+log(det(results$wi))-sum(diag(results$wi%*%amatrix[[1]]))/nrow(dd)-2*nn*log(mean(tau))/nrow(dd)
    )
  #likelihood1=1
}else{
    prob=Problem(Minimize(obj))
    result=CVXR::psolve(prob)
    S_est= lapply(B, function(x) result$getValue(x))
    likelihood=-(result$value-sum(abs(mask1*S_est[[1]])+abs(mask1*S_est[[2]]))-sum(mask2*abs(S_est[[1]]-S_est[[2]])))
    # likelihood1=ifelse(!random,
    #                   -extra/nrow(data[[1]])+log(det(S_est[[1]]))-sum(diag(S_est[[1]]%*%amatrix[[1]]))/nrow(data[[1]])
    #                   -extra/nrow(data[[2]])+log(det(S_est[[2]]))-sum(diag(S_est[[2]]%*%amatrix[[2]]))/nrow(data[[2]])
    #                   ,
    #                   -extra/nrow(data[[1]])+log(det(S_est[[1]]))-sum(diag(S_est[[1]]%*%amatrix[[1]]))/nrow(data[[1]])
    #                   -2*nn*log(mean(tau))/nrow(data[[1]])
    #                   -extra/nrow(data[[2]])+log(det(S_est[[2]]))-sum(diag(S_est[[2]]%*%amatrix[[2]]))/nrow(data[[2]])
    #                   -2*nn*log(mean(tau))/nrow(data[[2]])
    #                   )
    for (i in 1:length(S_est)) {
      colnames(S_est[[i]])=featureNames
      rownames(S_est[[i]])=featureNames
    }
}
    names(S_est)=glev
    return(list(wiList=S_est,ll=likelihood/2))

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
#' @param N a integer specifying the number of sampling for heterogeneous model
#' @param random a logical variable specifying the type of the model
#' @param ... other inputs
#' @import glasso glasso
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
                 tol=10^(-2),lower=c(0.01,0.1),upper=c(10,5), start=c("cold","warm"),
                 w.init=NULL, wi.init=NULL,trace=FALSE,...)

  {
  p=ncol(data)-2
  X_bar = apply(data[,-c(1,2)], 2, mean)
  data[,-c(1,2)] = scale(data[,-c(1,2)], center = X_bar, scale = FALSE)
data[,1]=as.character(data[,1])
  if (random==FALSE){

    if (!is.null(group))  {
group=as.character(group)
      if (length(group)!=nrow(data)){
        stop("group should be equal to number of the rows of data!")
      }

      data=split(data,f=factor(group,levels = unique(group)))

      if (!length(lambda)==2){
        stop("Arguments (group, lambda) do not match!")
      }

    }

if (is.null(group))  {
  group=as.character(rep(1,nrow(data)))
   data=list(data)
 if (length(lambda)!=1){
  stop("Arguments (group, lambda) do not match!!")
 }
}
    glev=unique(group)

  if (!all(lambda>0)){
    stop("Tuning parameter lambda must be positive!")
  }

  start=match.arg(start)
  # Create a mask matrix
  mask <- matrix(1, p, p)
  diag(mask) <- 0

    if (is.null(expFix)  | !is.numeric(expFix)){
      stop("Argument expFix is not correctly specified!")
    }
    A=vector("list",length(data))
    names(A)=glev
    B=vector("list",length((data)))
names(B)=glev
    for (i in 1:length(A)) {
      dd=data[[i]]
      subjects=unique(dd[,1])
      A[[i]]=vector("list",length(subjects))
      names(A[[i]])=subjects
      for (j in 1:length(A[[i]])) {
        index=which(dd[,1]==subjects[j])
        # if (length(index)==0){
        #   A[[i]][[j]]=NULL
        # }else{
        A[[i]][[j]]=diag(length(index))
        #}
      }
      B[[i]]=diag(p)
    }
k=0
tau0=0.5
while(1){
  k=k+1
A1=AA(data = data,B = B,expFix=expFix)
B1=BB(data=data,A=A,lambda = lambda,tau=A1$tau)

d1=c()
d2=c()
for (i in 1:length(B1[[1]])) {
    d1=c(d1,round(max(mask*abs(B[[i]]-B1$wiList[[i]])),3))
    d2=c(d2,round(abs(tau0-A1$tau),3))
  }
if (trace){
  print(paste0("iteration ",k, " precision difference: ",max(d1) , " /correlation tau difference: ",max(d2)))
}

if (max(d1)<=tol && max(d2)<= tol ){
    output=structure(list(wi=B1$wiList, v=A1$corMatrix, tau=A1$tau,ll=B1$ll), class="lglasso")
  break
}else{
  A=A1$corMatrix
  B=B1$wiList
  tau0=A1$tau
}

if (k>=maxit){
  message("Algorithm reached the maximum iteration!")
  output=structure(list(wi=B1$wiList, v=A1$corMatrix, tau=tau0,ll=B1$ll), class="lglasso")
  break
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

#' density function in EM algorithm
#'
#' @param tau the dampening rate
#' @param datai the data set for subject i
#' @param wi the given precision matrices
#' @param alpha the  rate in exponential distribution
#' @param groupi the data point indices
#' @param expFix a scalar specifying the form of the correlation function
#' @returns a numeric standing for the likelihood for a given subject

conDensityTau=function(tau,expFix=1, datai,wi,alpha,groupi){
    if (length(groupi)!=nrow(datai)){
      stop("group should be the same length of the columns of data!")
    }
  # if (!all(is.numeric(groupi))){
  #   stop("groupi need to be numeric!" )
  # }
  #browser()
  groupi=as.character(groupi)
  glev=unique(groupi)
  p=nrow(wi[[1]])
  #browser()
  if (length(glev)==1){
    timepoints=datai[,2]
    phiMi=phifunction(t=timepoints,tau=tau,expFix=expFix)
    s=t(as.matrix(datai[,-c(1,2)]))%*%solve(phiMi)%*%as.matrix(datai[,-c(1,2)])%*%wi[[glev]]
    likelihood=log(det(phiMi))*(-p/2)-0.5*sum(diag(s))
  }else{
    likelihood=0
  data=split(datai,f=factor(groupi,levels = glev))
  for (i in glev) {
    timepoints=data[[i]][,2]
    phiMi=phifunction(t=timepoints,tau=tau,expFix=expFix)
    s=t(as.matrix(data[[i]][,-c(1,2)]))%*%solve(phiMi)%*%as.matrix(data[[i]][,-c(1,2)])%*%wi[[i]]
    a=log(det(phiMi))*(-p/2)-0.5*sum(diag(s))
    likelihood=a+likelihood
  }
  }
  likelihoodi=log(alpha)+likelihood-alpha*tau
  return(likelihoodi=likelihoodi)
}

#' Function for generating the samples from posterior distribution in EM algorithm
#'
#' @param n the number of random samples
#' @param datai the data for subject i
#' @param wi the given precision matrix
#' @param alpha the exponential distribution with rate alpha
#' @param groupi specify how datai is grouped
#' @param expFix a scalar specifying the form of the correlation function
#' @returns a data frame for samples and their weights

importanceSample=function(n,datai,wi,alpha,groupi,expFix=1){
  dd1=matrix(rexp(n=n,rate=alpha),ncol=1)
  likelihood1=apply(dd1, 1, conDensityTau,datai=datai,wi=wi,alpha=alpha,groupi=groupi)
  likelihood2=log(apply(dd1, 1, dexp,rate=alpha))
  index1=which(!is.nan(likelihood1))
  index2=which(!is.infinite(likelihood1))
index=intersect(index1,index2)
if (length(index)==0){stop("No valid samples are generated!")}
  weights=likelihood1[index]-likelihood2[index]
  weightNew=weights-(log(sum(exp(weights-max(weights))))+max(weights))
  norWeights=exp(weightNew)
  aa=data.frame(sample=dd1[index],weight=norWeights)
  return(aa)
}

#' Function for computing estimates in EM algorithm
#'
#' @param importancesSample the random samples
#' @param datai the data for subject i
#' @param groupi specify how datai is grouped
#' @param expFix a scalar specifying the form of the correlation function
#' @returns a list of estimated

importanceEstimates=function(importancesSample,datai,groupi,expFix=1){
  lev=as.character(unique(groupi))
  # if (length(lev)==1){
  #   estimatePhi=vector("list",1)
  #   tt=datai[,2]
  #   sample0=importancesSample[,1]
  #   weightTau=importancesSample[,2]
  #   estimateTau=sum(sample0*weightTau)
  #   estimatePhi[[1]]=phifunction(t=tt,tau=estimateTau, expFix=expFix)
  # }else{
  tt=split(datai[,2],f=factor(groupi,levels = lev))
  sample0=importancesSample[,1]
  weightTau=importancesSample[,2]
  estimateTau=sum(sample0*weightTau)
  estimatePhi=lapply(tt, phifunction,tau=estimateTau, expFix=expFix)
  names(estimatePhi)=names(tt)
  estimates=list(estimateTau=estimateTau,estimatePhi=estimatePhi)
  return(estimates)
}


#' Estimate the phimatrix in heterogeneous model
#' @param data longitudinal data set
#' @param wi given precision matrix
#' @param alpha exponential distribution with rate alpha
#' @param group specify how data is grouped
#' @param l number of random samples in importance sampling
#' @param expFix a scalar specifying the the form of the correlation function.
#' @param ... other arguments used in the downstream analysis
#' @returns a list for estimates of tau and AA
AAheter=function(data,wi,alpha,group,l=5000,expFix=1,...){
  data[,1]=as.character(data[,1])
  subjects=unique(data[,1])
  group=as.character(group)
  glev=unique(group)
  nn=length(glev)
  A=vector("list",length(subjects))
  names(A)=subjects
  Tau=matrix(nrow=length(subjects),ncol=1)
  rownames(Tau)=subjects
  simTau=matrix(rexp(n=l,rate=alpha),ncol=1)
  dataList=split(data,f=factor(data[,1],levels=subjects))
  groupList=split(group,f=factor(data[,1],levels=subjects))
  for (i in 1:length(subjects)) {
    datai=dataList[[i]]
    groupi=groupList[[i]]
      imSample=importanceSample(n=l,datai=datai,wi=wi,alpha=alpha,groupi =groupi,expFix=expFix )
      index=which(!is.nan(imSample[,2]))
      imSample=imSample[index,]
      imporResults=importanceEstimates(importancesSample=imSample,datai=datai,
                                       groupi = groupi,expFix=expFix,...)
      Tau[i,1]=imporResults$estimateTau
      A[[i]]=imporResults$estimatePhi
  }

  AA=vector("list",nn)
  names(AA)=glev
#if (nn==1){AA[[1]]=A}
  for(i in 1:nn) {
    index=which(group==glev[i])
    subjects=unique(data[index,1])
    AA[[glev[i]]]=vector("list",length(subjects))
    names(AA[[glev[i]]])=subjects
    for (j in 1:length(subjects)) {
      if (is.null(A[[subjects[j]]][[glev[i]]])){next}
        else{
      AA[[glev[i]]][[subjects[j]]]=A[[subjects[j]]][[glev[i]]]
      }
    }
  # for (i in 1:nn) {
  #   for (j in 1:length(subjects)) {
  #     AA[[i]][[j]]=A[[j]][[i]]
  #   }
  # }
}
  return(list(Tau=Tau,AA=AA))
}



#' Title
#'
#' @param data a n by (p+2) data frame representing the longitudinal data
#' @param lambda tuning parameters
#' @param group vector indicating the membership of each data point.
#' @param maxit the maximum number of the iterations
#' @param tol the lower bound which determine when the algorithm is thought to reach  convergence.
#' @param trace  a logical variable specifying how the output is displayed on the screen
#' @param start a binary variable specifying how the initial value of the algorithm is chosen.
#' @param w.init the initial value for the covariance matrix
#' @param wi.init the initial value for the precision matrix
#' @param N the number of sampling for heterogeneous model
#' @param expFix a scalar specifying the form of the correlation function.
#' @param ... other arguments
#' @returns a list of length 4 representing the final outcome

lglassoHeter=function(data,lambda,group,maxit,
                      tol=10^(-3),trace=FALSE,start=c("warm","cold"),
                      w.init=NULL, wi.init=NULL, N,expFix=1,...)

{
  p=ncol(data)-2
  m=length(unique(data[,1]))
  if (!is.null(group)){
    group=as.character(group)
  }else{
         group=as.character(rep(1,nrow(data)))
         if (length(lambda)!=1){
             stop("Arguments (group, lambda) do not match!")
           }
       }
data[,1]=as.character(data[,1])
         dataList=split(data,f=factor(group,levels = unique(group)))

         if (length(lambda)!= length(unique(group))){
             stop("Arguments (group, lambda) do not match!")
           }

glev=unique(group)

  if (!all(lambda>0)){
    stop("lambda must be positive!")
  }

  # Create a mask matrix
  mask <- matrix(1, p, p)
  diag(mask) <- 0


  if (is.null(expFix) | !is.numeric(expFix)){
    stop("Argument expFix is not correctly specified!")
  }

  B=vector("list",length(glev))
names(B)=glev
  for (i in 1:length(glev)) {
    B[[i]]=diag(p)
  }
  k=0
  tau0=rep(1,m)
  alpha0=1
  while(1){
    k=k+1
    A1=AAheter(data=data,wi=B,alpha=alpha0,group=group,expFix=expFix,l=N,...)
    B1=BB(data=dataList,A=A1$AA,lambda = lambda,random=T,tau=A1$Tau,...)
    d1=c()
    d2=round(abs(alpha0-1/mean(A1$Tau)),3)
    for (i in 1:length(B1[[1]])) {
      d1=c(d1,round(max(mask*abs(B[[i]]-B1$wiList[[i]])),3))
    }
    if (trace){
      print(paste0("alpha estimate: ", alpha0))
      print(paste0("iteration ",k, " precision difference: ",max(d1) , " /correlation alpha difference: ",max(d2)))
    }

    if (max(d1)<=tol && d2<= tol ){
      output=structure(list(wi=B1$wiList, tau=A1$Tau,alpha=1/mean(A1$Tau),ll=B1$ll), class="lglasso")
      break
    }else{
      A=A1$AA
      B=B1$wiList
      tau0=A1$Tau
      alpha0=1/mean(tau0)
    }

    if (k>=maxit){
      message("Algorithm reached the maximum iteration!")
      output=structure(list(wi=B1$wiList, tau=tau0,alpha=alpha0,ll=B1$ll), class="lglasso")
      break
    }
  }
  return(output)
}



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
      coef.train=stats::lm(y~x)$coef
      yy=data.valid[,i+2, drop=FALSE]
      xx=as.matrix(cbind(1,data.valid[,index+2,drop=FALSE]))
      err=(yy-xx%*%coef.train)^2
      cv_error=c(cv_error,mean(err[,,drop=TRUE]))
      if (any(is.na(cv_error))) {
        print("cv error is missing!")
    }

  return(mean(cv_error))
}

}

cvErrorj=function(data.train,data.valid,B){
    a= mean(apply(B, 2, function(bi) cvErrorji(data.train=data.train,data.valid=data.valid,bi=bi)))
}




#' Compute the cross validation error
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









#' Plot function for CVlglasso
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
  #err=apply(x$cv_error, 1, mean)
  err=x$cv_error
  if (is.vector(lambda)){
    graphics::plot(x=x_data,y=err,xlab=xlab_label,ylab="CV Error", main = "cvlglasso Fit",
                   type = "b", ...)
  }else{
lambda=as.matrix(lambda)
a1=sort(unique(lambda[,1]),decreasing = F)
a2=sort(unique(lambda[,2]),decreasing = F)
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







#' @title Cross validation for \code{lglasso}
#' @description
#' The function computes the cross validation errors for one of the three network models in \code{lglasso} command.
#' @param data raw data
#' @param group group variable
#' @param lambda tuning parameter
#' @param nlam number of tuning parameter
#' @param lam.min.ratio ratio of largest lambda vs smallest lambda
#' @param K cv folds
#' @param expFix given parameter
#' @param trace whether show the process
#' @param NN the number of sampling
#' @param random a logical variable indicating the type of the model
#' @returns list of which the first component is the cross validation errors and the second component is the corresponding
#' tuning parameters
#' @export
#'
CVlglasso=function(data,group=NULL,random=FALSE,
                    lambda=NULL,nlam=10,lam.min.ratio=0.01, K, expFix=1,trace=FALSE,NN=500){

  results=cvlglassofull(data=data,group=group,lambda = lambda,nlam=nlam,random = random,
                        lam.min.ratio=lam.min.ratio, K=K, expFix=expFix,trace=trace, NN=NN)

  return(results)
}



#' Cross validation for lglasso
#'
#' @param data raw data
#' @param group group variable
#' @param lambda tuning parameter
#' @param nlam number of tuning parameter
#' @param lam.min.ratio ratio of largest lambda vs smallest lambda
#' @param K cv folds
#' @param expFix given parameter
#' @param trace whether show the process
#' @param random a logical variable specifying the type of the model
#' @param NN a integer specifying the number of the sampling
#' @returns list
#' @import parallel foreach doParallel

cvlglassofull=function(data,group=NULL,
                    lambda=NULL,random=FALSE,nlam=10,lam.min.ratio=0.01,
                    K, expFix=1,trace=FALSE,NN){

if (!is.null(lambda)){
  if (is.null(group) && !is.vector(lambda))
  {stop("group and lambda does not match!")}

  if (!is.null(group) && ncol(lambda)!=2)
  {stop("group and lambda does not match!")}
  if (any(lambda<=0)) {stop("tuning parameter lambda should be positive!")}

}

  if (any(K<=1 | K%%1 !=0)){
    stop("K should be an integer greater than 1!")
  }

  cores=detectCores()


  if (is.null(lambda)){
    if (is.null(group)){
      N=K*nlam
    }else{
      N=K*nlam^2
    }
  }else{
N=ifelse(is.vector(lambda),K*length(lambda),K*nrow(lambda))
  }



  if (cores > N) {cores = N}
  cat("\nNumber of cores used =", cores, "\n")

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
  }else {
    if (is.null(group)){
      lambda = sort(lambda,decreasing = FALSE)
    }
  }


  nnlambda=ifelse(is.null(group),length(lambda),nrow(lambda))
  cv_error=matrix(0,nrow=nnlambda,ncol=K)
  crossData=vector("list",K)
  names(crossData)=paste0("CVdata",1:K)
  for (k in seq_len(K)) {
    leave.out =subjects[ind[(1 + floor((k - 1) * n/K)):floor(k *
                                                               n/K)]]
    indexValid=which(data[,1] %in% leave.out)
    data.train = data[-indexValid, , drop = FALSE]
    data_bar = apply(data.train[,-c(1,2)], 2, mean)
    data.train[,-c(1,2)] = scale(data.train[,-c(1,2)], center = data_bar, scale = FALSE)
    data.valid = data[indexValid,, drop = FALSE]
    data.valid[,-c(1,2)] = scale(data.valid[,-c(1,2)], center = data_bar, scale = FALSE)


    if(is.null(group)){
      crossData[[k]]=list(train=data.train,
                        validation=data.valid)
    }else{
   crossData[[k]]=list(train=data.train,
                          validation=data.valid,
                          trainGroup=group[-indexValid],
                          validationGroup=group[indexValid])
      }
  }



crossDataLambda=vector("list",N)
  for (j1 in 1:nnlambda) {
    for (j2 in 1:K){
      i=(j1-1)*K+j2
      if (is.null(group)){
        LL=lambda[j1]
      }else{
        LL=unlist(lambda[j1,])
      }
      crossDataLambda[[i]]=list(lambda=LL,crossData=crossData[[j2]])
    }
  }


  k=NULL
        CV = foreach(k = 1:length(crossDataLambda), .packages = "lglasso", .combine = "cbind",
                     .inorder = TRUE) %dopar% {
                 if (trace) {
                   progress = utils::txtProgressBar(max = N, style = 3)
                 }
                 if (is.null(group)){
                   if (random==FALSE){
                   aa= lglasso(data=crossDataLambda[[k]]$crossData$train,lambda=crossDataLambda[[k]]$lambda,N=NN)$wi[[1]]
                   }else{
                    aa= lglasso(data=crossDataLambda[[k]]$crossData$train,lambda=crossDataLambda[[k]]$lambda,random = TRUE,N=NN)$wi[[1]]
                   }
                   cc=ifelse(abs(aa)<=10^(-2), 0,1)
                   diag(cc)=2
                   bb= cvError(data.train=crossDataLambda[[k]]$crossData$train,
                               data.valid=crossDataLambda[[k]]$crossData$validation,
                               B=cc)
                 }
                 if (!is.null(group)){
                   if (random==FALSE){
                   aa=lglasso(data=crossDataLambda[[k]]$crossData$train,
                              lambda=crossDataLambda[[k]]$lambda,
                              expFix = expFix,
                              group=crossDataLambda[[k]]$crossData$trainGroup,N=NN)$wi
                   }else{
                     aa=lglasso(data=crossDataLambda[[k]]$crossData$train,
                                lambda=crossDataLambda[[k]]$lambda,
                                expFix = expFix,
                                group=crossDataLambda[[k]]$crossData$trainGroup,
                                random = TRUE,N=NN)$wi
                   }
                   aa[[1]]=ifelse(abs(aa[[1]])<=10^(-1), 0,1)
                   diag(aa[[1]])=2
                   aa[[2]]=ifelse(abs(aa[[2]])<=10^(-1), 0,1)
                   diag(aa[[2]])=2
                   cc=aa
                   bb=cvError(data.train=crossDataLambda[[k]]$crossData$train,
                              data.valid=crossDataLambda[[k]]$crossData$validation,
                              B=cc,
                              group.valid=crossDataLambda[[k]]$crossData$validationGroup,
                              group.train = crossDataLambda[[k]]$crossData$trainGroup)
                 }

                 cv_error=bb

                 if (trace) {
                   utils::setTxtProgressBar(progress,  k)
                 }

                 return(cv_error=cv_error)
                 }
    aa=c()
    for (i in 1:nnlambda) {
     aa=c(aa, mean(CV[((i-1)*K+1):(i*K)]))
    }
    names(aa)=paste0("lambda",1:nnlambda)

  stopCluster(cluster)
  output=structure(list(cv_error=aa, lambda=lambda),class="cvlglasso")
  return(output)
}

