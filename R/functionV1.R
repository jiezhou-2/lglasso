


#' @title Longitudinal graphical lasso
#' @description
#'  This function estimates  precision matrices(networks) and random effects from longitudinal
#'  high-dimensional data under normality assumption.
#' @param data \code{n} by \code{(p+2)} data frame in which the first column is subject IDs, the second column is
#' the time points of longitudinal data.
#' @param lambda   numerical vector of tuning parameters. For one-stage model, \code{lambda} is a scalar controlling the
#' sparsity of the network. For two-stage model,
#' \code{lambda} is a vector of length 2, in which the first entry controls the sparsity of both
#' pre-treatment and post-treatment networks,
#' while the second entry controls the overlap of the two networks.
#' @param group  factor  of length \code{n} if supplied. If \code{group} is a  one-level factor,
#' then a one-stage model is fitted. If \code{group} is a two-level factor, then a two-stage model is fitted.
#'  Data points that are before (after) the treatment (exposure) share the same
#'   level in \code{group}. If NULL, then a one-stage model is fitted.
#' @param random a logical variable. If TRUE, then a heterogeneous model is fitted.
#' Otherwise, a homogeneous model is fitted.
#' @param expFix  numerical number specifying the exponent in the covariance function of
#' the longitudinal data. Default is 1.
#' @param maxit integer specifying  the maximum iterations for the algorithms.
#' @param tol   a small number determining whether the algorithms converged and should stop.
#' @param lower  vector of length 1 or 2 which specifies the lower bounds for temporal correlation \code{tau}.
#' It is of length 1 for one-stage model and 2 for two-stage model
#' @param upper  vector of length 1 or 2 which specifies the upper bounds for temporal correlation \code{tau}
#' It is of length 1 for one-stage model and 2 for two-stage model
#' @param w.init initial value for covariance matrix. Default is identity matrix.
#' @param wi.init initial value for precision matrix. Default is identity matrix.
#' @param trace whether or not show the progress of the computation
#' @param N a integer specifying the number of sampling for heterogeneous model
#' @param ... other inputs
#' @export
#' @example inst/examples.R
#' @return A list which includes:
#'
#' \code{w:} the list of the estimates for covariance matrices
#'
#' \code{wi:} the list of the estimates for precision matrices
#'
#' \code{tau:} the estimate of dampening rate \code{tau}. For heterogeneous models, the
#' output is a vector called random effects. For homogeneous model, the output is a scalar.
#'
#' \code{alpha:} a scalar representing the estimate of the parameter in exponential distribution of \code{tau} for heterogeneous models.
#' The output is NULL for homogeneous model
#'
#'\code{ll:} the likelihood.  \code{ll} is used to compute extended QIC for tuning parameter selection
#'
#' @details \code{lglasso} is developed to estimate precision matrices, or networks, and random effects \eqn{\tau_i}
#' from high-dimensional longitudinal data.
#'  The one-stage model in \code{lglasso} is proposed in Zhou *et al* (2024).
#'  Currently, it contains two network identification models,
#'   *i.e.,* one-stage  and two-stage model.
#'  One-stage model assume a common network underlying
#'   the longitudinal data for all the subjects. Consequently,
#'   the function outputs a single network as the estimate in this case.
#'   In two-stage models, a treatment is applied at time \eqn{t_i} during the time interval
#'     for subject \eqn{i}.
#'  Therefore, there are two networks, *i.e*., pre- and post-treatment networks,  that need to be estimated simultaneously.
#'  For details, please check the reference paper
#'    and the online resources.
lglasso=function(data,lambda,group=NULL,random=FALSE,expFix=1,N=100,maxit=50,
                 tol=10^(-2),lower=c(0.01,0.01),upper=c(10,10),
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
        A[[i]][[j]]=diag(length(index))
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
    d1=c(d1,round(max(mask*abs(B[[i]]-B1$wi[[i]])),3))
    d2=c(d2,round(abs(tau0-A1$tau),3))
  }
if (trace){
  print(paste0("iteration ",k, " precision difference: ",max(d1) , " /correlation tau difference: ",max(d2)))
}

if (max(d1)<=tol && max(d2)<= tol ){
    output=structure(list(wi=B1$wi,w=B1$w, tau=A1$tau,alpha=NA,ll=B1$ll), class="lglasso")
  break
}else{
  A=A1$corMatrix
  B=B1$wi
  tau0=A1$tau
}

if (k>=maxit){
  message("Algorithm reached the maximum iteration!")
  output=structure(list(wi=B1$wi, w=B1$w, tau=tau0,alpha=NA,ll=B1$ll), class="lglasso")
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

output=lglassoHeter(data=data,lambda=lambda,expFix=expFix,
                    N=N,group=group,maxit=maxit,trace=trace)
  }
  return(output)
}

#' Title
#'
#' @param data a n by (p+2) data frame representing the longitudinal data
#' @param lambda tuning parameters
#' @param group vector indicating the membership of each data point.
#' @param maxit the maximum number of the iterations
#' @param tol the lower bound which determine when the algorithm is thought to reach  convergence.
#' @param trace  a logical variable specifying how the output is displayed on the screen
#' @param w.init the initial value for the covariance matrix
#' @param wi.init the initial value for the precision matrix
#' @param N the number of sampling for heterogeneous model
#' @param expFix a scalar specifying the form of the correlation function.
#' @param ... other arguments
#' @noRd
#' @returns a list of length 4 representing the final outcome

lglassoHeter=function(data,lambda,group,maxit=50,
                      tol=10^(-1),trace=FALSE,
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
      d1=c(d1,round(max(mask*abs(B[[i]]-B1$wi[[i]])),3))
    }
    if (trace){
      print(paste0("alpha estimate: ", alpha0))
      print(paste0("iteration ",k, " precision difference: ",max(d1) , " /correlation alpha difference: ",max(d2)))
    }

    if (max(d1)<=tol && d2<= tol ){
      output=structure(list(wi=B1$wi,w=B1$w, tau=A1$Tau,alpha=1/mean(A1$Tau),ll=B1$ll), class="lglasso")
      break
    }else{
      A=A1$AA
      B=B1$wi
      tau0=A1$Tau
      alpha0=1/mean(tau0)
    }

    if (k>=maxit){
      message("Algorithm reached the maximum iteration!")
      output=structure(list(wi=B1$wi,w=B1$w, tau=tau0,alpha=alpha0,ll=B1$ll), class="lglasso")
      break
    }
  }
  return(output)
}


#' @noRd
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
      }else{
      y=data.train[,i+2]
      x=as.matrix(data.train[,index+2])
      coef.train=stats::lm(y~x)$coef
      yy=data.valid[,i+2, drop=FALSE]
      xx=as.matrix(cbind(1,data.valid[,index+2,drop=FALSE]))
      err=(yy-xx%*%coef.train)^2
      cv_error=mean(err[,,drop=TRUE])
      }
      if (any(is.na(cv_error))) {
        print("cv error is missing!")
    }

  return(cv_error)
}

}

#' @noRd
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
#' @noRd
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
#' @param ... other plot arguments
#' @noRd
#' @returns If \code{group} is NULL in \code{CVlglasso}, then a line plot will produced; otherwise, a heatmap will be produced.
#' @export
plot.cvlglasso=function(x,...){

  if (!inherits(x, "cvlglasso")) {
    stop("x must be an object of class 'cvlglasso'")
  }

  lambda=x$lambda
  err=x$cv_error
  if (is.vector(lambda)){
    graphics::plot(x=lambda,y=err,xlab="tuning parameter",ylab="CV Error", main = "cvlglasso Fit",
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
#' The function computes the cross validation errors  in \code{lglasso}.
#' @param data same as in \code{lglasso}
#' @param K fold number for the cross validation
#' @param group same as in \code{lglasso}
#' @param lambda tuning parameter. For one-stage model, lambda is a vector. For two-stage model,
#'  lambda is a 2-column matrix. The first column  is tuning parameter controlling the sparsity,
#'  the second column is tuning parameter controlling the
#'  similarity of two networks.
#' @param nlam If \code{lambda} is NULL, then \code{nlam} set the number of tuning parameter
#' @param lam.min.ratio ratio of largest lambda vs smallest lambda
#' @param expFix same as in \code{lglasso}
#' @param trace logical variable. Whether show the computation process
#' @param random same as in \code{lglasso}
#' @export
#' @returns list of length 2.  The first component is the cross validation errors, the second component is the corresponding
#' tuning parameters


CVlglasso=function(data,K,group=NULL,random=FALSE,
                    lambda=NULL,nlam=10,lam.min.ratio=0.01, expFix=1,trace=FALSE){

  results=cvlglassofull(data=data,group=group,lambda = lambda,nlam=nlam,random = random,
                        lam.min.ratio=lam.min.ratio, K=K, expFix=expFix,trace=trace)

  return(results)
}



#' Cross validation for lglasso
#'
#' @param data data used in lglasso
#' @param group group variable used in lglasso
#' @param lambda tuning parameter. For one-stage model, lambda is a vector. For two-stage model,
#'  lambda is a matrix. The first entry is to control the sparsity, the second entry is to control the
#'  similarity of two networks.
#' @param nlam If lambda is NULL, then nlam set the number of tuning parameter
#' @param lam.min.ratio ratio of largest lambda vs smallest lambda
#' @param K cv folds
#' @param expFix given parameter
#' @param trace whether show the process
#' @param random a logical variable specifying the type of the model
#' @returns list
#' @noRd
#' @import parallel foreach doParallel

cvlglassofull=function(data,group=NULL,
                    lambda=NULL,random=FALSE,nlam=10,lam.min.ratio=0.01,
                    K, expFix=1,trace=FALSE){

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
                   aa= lglasso(data=crossDataLambda[[k]]$crossData$train,lambda=crossDataLambda[[k]]$lambda)$wi[[1]]
                   }else{
                    aa= lglasso(data=crossDataLambda[[k]]$crossData$train,lambda=crossDataLambda[[k]]$lambda,random = TRUE)$wi[[1]]
                   }
                   cc=ifelse(abs(aa)<=10^(-5), 0,1)
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
                              group=crossDataLambda[[k]]$crossData$trainGroup)$wi
                   }else{
                     aa=lglasso(data=crossDataLambda[[k]]$crossData$train,
                                lambda=crossDataLambda[[k]]$lambda,
                                expFix = expFix,
                                group=crossDataLambda[[k]]$crossData$trainGroup,
                                random = TRUE)$wi
                   }
                   aa[[1]]=ifelse(abs(aa[[1]])<=10^(-5), 0,1)
                   diag(aa[[1]])=2
                   aa[[2]]=ifelse(abs(aa[[2]])<=10^(-5), 0,1)
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



#' Simulate longitudinal data from one-stage/two-stage model
#' @description This function randomly generates precision matrices (networks)
#'  and then  simulates normal data that follows these network structures.
#' @param type which type of data you are generating. There are two options. One is \code{homo} which generates subjects with identical
#' temporal correlation parameter. The other is \code{heter} which generates subjects with different temporal correlation parameter.
#' @param n the number of subjects in the data set
#' @param p the dimension of the normal distribution
#' @param m1 the number of edges in true networks
#' @param m2 the edge difference between two networks
#' @param tt the average time points for each subject
#' @param tau the true dampening rate in homogeneous models
#' @param alpha the true parameter in exponential distribution of tau when \code{type}
#' is \code{heter}
#' @param group a scalar of 1 or 2,  indicating one-stage or two-stage model.
#' @export
#' @returns a data list. It includes the data generated, true networks underlying the data,true tau.
#' If \code{type} is \code{heter}, then true parameter \code{alpha} is included as well.

Simulate=function(type=c("homo","heter"),n,p,m1,
                  m2,tt,tau,alpha=1,group=1){

  type=match.arg(type)
  if (type=="homo"){
    data=simulate_long(n=n,p=p,m1=m1,m2=m2,tau=tau,tt=tt)
  }

  if (type=="heter"){

    data=simulate_randomTau(n=n, p=p,m1=m1,m2=m2,tt=tt,
                            alpha=alpha,group = group)
  }

  return(data)
}






