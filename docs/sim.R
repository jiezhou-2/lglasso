# library(fake)
# library(CVXR)
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
        M[i,j]=exp(-tau[1]*(abs(t[i]-t[j])))
        M[j,i]=M[i,j]
      }
    }
    diag(M)=1
  }
  return(M)
}

sim_stru=function(p,m1,m2){
  real_stru=matrix(0, nrow = p, ncol = p)
  real_stru[lower.tri(real_stru,diag = TRUE)]=1
  index=which(real_stru==0,arr.ind = TRUE)
  a=sample(1:nrow(index),m1, replace = F)
  real_stru[index[a,]]=1
  real_stru[lower.tri(real_stru,diag = TRUE)]=0
  real_stru1=real_stru+t(real_stru)+diag(p)


  distrubance=matrix(0, nrow = p, ncol = p)
  distrubance[lower.tri(distrubance,diag = TRUE)]=1
  index=which(distrubance==0,arr.ind = TRUE)
  a=sample(1:nrow(index),m2, replace = F)
  distrubance[index[a,]]=1
  distrubance[lower.tri(distrubance,diag = TRUE)]=0
  distrubance=distrubance+t(distrubance)+diag(p)
  real_stru2=(real_stru1+distrubance)%%2

  theta = matrix(stats::rnorm(p^2,mean = 0,sd=2), ncol = p,nrow = p)
  theta[lower.tri(theta, diag = TRUE)] = 0
  theta = theta + t(theta) + diag(p)
  theta1 = theta * real_stru1
  theta2 = theta * real_stru2
  theta1=MakePositiveDefinite(theta1,pd_strategy = "diagonally_dominant",scale = TRUE)$omega
  theta2=MakePositiveDefinite(theta2,pd_strategy = "diagonally_dominant",scale = TRUE)$omega
  sigma1=solve(theta1)
  sigma2=solve(theta2)
  return(list(precision1=theta1,precision2=theta2,covmat1=sigma1,covmat2=sigma2))
}

sim_timepoints=function(n,tt){
    ## true structure
  timepoint1=vector("list",n)
  timepoint2=vector("list",n)
  for (i in 1:n) {
    m3=sample(x=1:tt[1],1,prob = rep(1,tt[1]))
    m4=sample(x=1:tt[2],1,prob = rep(1,tt[2]))
    t1=stats::rexp(m3,rate=10)
    t2=stats::rexp(m4,rate=10)
    timepoint1[[i]]=cumsum(t1)[1:m3]
    timepoint2[[i]]=cumsum(t2)[1:m4]
  }
  return(list(timepoint1=timepoint1,timepoint2=timepoint2))
}



sim_phi=function(timepoints,tau){
  phi1=lapply(timepoints[[1]],phifunction, tau=tau[1])
  phi2=lapply(timepoints[[2]],phifunction, tau=tau[2])
  return(list(phi1=phi1,phi2=phi2))
}


sim_data_old=function(covmat,timepoints,tau){
  print("homo data are generated")
  if (length(timepoints[[1]])!=length(timepoints[[2]])){
    stop("pretreatment and posttreatment should have same length time points!")
  }
  p=ncol(covmat[[1]])
  m=length(timepoints[[1]])
  data1=vector("list",m)
  data2=vector("list",m)
  sqK1=chol(covmat[[1]])
  sqK2=chol(covmat[[2]])
  age1=timepoints[[1]]
  age2=timepoints[[2]]
  for (i in 1:m) {
    n=length(age1[[i]])
    a=matrix(ncol = p,nrow = n)
    error1=matrix(rnorm(p*n),nrow = p)
    error2=t(t(sqK1)%*%error1)
    a[1,]=error2[1,]
    if (n>1){
      for (t in 2:n) {
        coe=exp(-tau[1]*abs(age1[[i]][t]-age1[[i]][t-1]))
        #coe=0
        a[t,]=a[t-1,]*coe+error2[t,]*sqrt(1-coe^2)
      }
    }
    #dd=cenfunction(a,zirate = zirate)
    data1[[i]]=cbind(i,age1[[i]],a)
    colnames(data1[[i]])[1:2]=c("subject","time")
  }
  dd1=as.data.frame(do.call(rbind,data1))

  for (i in 1:m) {
    n=length(age2[[i]])
    a=matrix(ncol = p,nrow = n)
    error1=matrix(rnorm(p*n),nrow = p)
    error2=t(t(sqK2)%*%error1)
    a[1,]=error2[1,]
    if (n>1){
      for (t in 2:n) {
        coe=exp(-tau[2]*abs(age2[[i]][t]-age2[[i]][t-1]))
        a[t,]=a[t-1,]*coe+error2[t,]*sqrt(1-coe^2)
      }
    }
    #dd=cenfunction(a,zirate = zirate)
    data2[[i]]=cbind(i,age2[[i]],a)
    colnames(data2[[i]])[1:2]=c("subject","time")
  }
  dd2=as.data.frame(do.call(rbind,data2))
  return(list(data=list(pre=dd1,post=dd2)))
}





sim_data=function(phi1,phi2,covmat, timepoints){

  fullCovariance1= lapply(phi1,function(x,cc) {kronecker(cc,x)},x=covmat[[1]])
    fullCovariance2= lapply(phi2,function(x,cc) {kronecker(cc,x)},x=covmat[[2]])

  fulldata=c()
  a1=c()
  a2=c()
  if (length(phi1)!=length(phi2)){stop("the length of phi1 and phi2 should be same!")}
  n=length(phi1)
  for (i in 1:n) {
    ai=c()
    m3=nrow(fullCovariance1[[i]])
    mu=rep(0,m3)
    data1=MASS::mvrnorm(1,mu=mu,Sigma = fullCovariance1[[i]])
    for (j in 1:(m3/p)) {
      ai=as.data.frame(rbind(ai,data1[c(((j-1)*p+1) :(j*p))]))
    }
    subject=paste0("subject",i)
    ai=cbind(subject,timepoints[[1]][[i]],ai)
    a1=rbind(a1,ai)
  }
  colnames(a1)[2]="time"


    for (i in 1:n) {
      ai=c()
      m3=nrow(fullCovariance2[[i]])
      mu=rep(0,m3)
      data2=MASS::mvrnorm(1,mu=mu,Sigma = fullCovariance2[[i]])
      for (j in 1:(m3/p)) {
        ai=as.data.frame(rbind(ai,data2[c(((j-1)*p+1) :(j*p))]))
      }

      subject=paste0("subject",i)
      ai=cbind(subject,timepoints[[2]][[i]],ai)
      a2=rbind(a2,ai)
    }
    colnames(a2)[2]="time"
    return(list(data=list(pre=a1,post=a2)))
}


BB=function(A,data,lambda,type="expFixed",diagonal=TRUE,maxit=100,
            tol=10^(-4),lower=c(0.01,0.1),upper=c(10,5), start=c("warm","cold"),w.init=NULL,wi.init=NULL,...){

  type=match.arg(type)
  start=match.arg(start)


  if (!is.list(data) | !is.list(A)){
    stop("A and data must be a list!")
  }

  if (length(A)!= length(data)){
    stop(" A should have same length as  data!")
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
        #browser()
        amatrix=t(xx)%*%xx+amatrix
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
      return(wi=S_est)
    }
  }
}


BB=function(A,data,lambda,random=FALSE,tau=rep(1,length(A))){
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

  m= length(A)
  for (i in 1:length(A)){
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
      {stop("The format of A does not match the format of data!!")}
    }

  }



  B=vector("list",m)
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
      amatrix[[i]]=t(xx)%*%yy%*%xx+amatrix[[i]]
      extra=extra+p*log(det(as.matrix(Ai[[j]])))
    }
    if (!random){
      likeli=-extra/nrow(dd)+log_det(B[[i]])-matrix_trace(B[[i]]%*%amatrix[[i]])/nrow(dd)+likeli
    }else{
      likeli=-extra/nrow(dd)+log_det(B[[i]])-matrix_trace(B[[i]]%*%amatrix[[i]])/nrow(dd)
      -2*nn*log(mean(tau))+likeli
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
                      -extra/nrow(dd)+log(det(results$wi)-sum(diag(results$wi%*%(amatrix[[1]]/(nrow(dd)))))),
                      -extra/nrow(dd)+log(det(results$wi)-sum(diag(results$wi%*%(amatrix[[1]]/(nrow(dd))))))-2*nn*log(mean(tau))
    )
  }else{
    prob=Problem(Minimize(obj))
    result=CVXR::solve(prob)
    S_est= lapply(B, function(x) result$getValue(x))
    likelihood=-(result$value
                 -lambda[1]*(sum(abs(S_est[[1]])+abs(S_est[[2]])))
                 -lambda[2]*(sum(abs(S_est[[1]]-S_est[[2]]))))
    for (i in 1:length(S_est)) {
      colnames(S_est[[i]])=featureNames
      rownames(S_est[[i]])=featureNames
    }
  }
  return(list(wiList=S_est,ll=likelihood))

}




cvBBfull=function(data,group=NULL,
                       lambda=NULL,random=FALSE,nlam=10,lam.min.ratio=0.01, K,
                       expFix=1,trace=FALSE,NN){

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
  #browser()
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



# BB=function(data,lambda,diagonal=TRUE,lower=c(0.01,0.1),upper=c(10,5)){
#     obj=0
#     aa=0
#     bb=0
#     B=vector("list",2)
#     p= ncol(data[[1]])-2
#     # Create a mask matrix
#     mask1 <- matrix(lambda[1], p, p)
#     diag(mask1) <- 0
#     mask2 <- matrix(lambda[2], p, p)
#     diag(mask2) <- 0
#     B[[1]]=Variable(p,p,PSD=TRUE) # tissue wise inverse correlation matrix
#     B[[2]]=Variable(p,p,PSD=TRUE)
#     for (i in 1:2){
#       dd=data[[i]]
#       nn=length(unique(dd[,1]))
#       data_sub=split(dd[,-c(1,2)],factor(dd[,1],unique(dd[,1])))
#       amatrix=0
#
#       for (j in 1:nn) {
#         xx=as.matrix(data_sub[[j]])
#         amatrix=t(xx)%*%xx+amatrix
#       }
#       obj=log_det(B[[i]])-matrix_trace(B[[i]]%*%amatrix)/nrow(dd)+obj
#       aa=aa+ sum(abs(B[[i]])*mask1)
#     }
#
#     bb=bb+sum(abs((B[[1]]-B[[2]]))*mask2)
#     obj=-obj+aa+bb
#       prob=Problem(Minimize(obj))
#       result=CVXR::solve(prob)
#       S_est= lapply(B, function(x) result$getValue(x))
#       return(wi=S_est)
# }

tunefinder_lglasso=function(phiM,sigmaM,timepoints,rho,
                         zirate=c(0.2,0),random,group,K,tau)
{
  if (length(rho_lglasso)!=2){
    stop("rho_lglasso should be length 2!")
  }



  simData=sim_data_old(covmat = sigmaM[3:4],timepoints=timepoints, tau = tau)
  #simData=sim_data(phi1 = phiM[[1]],phi2=phiM[[2]], covmat=sigmaM[3:4],timepoints = timepoints)
  preData=simData$data$pre
  postData=simData$data$post
  fullData=rbind(preData,postData)

  groupIndex=c(rep(1,nrow(preData)),rep(2,nrow(postData)))

  X_bar = apply(fullData[,-c(1,2)], 2, mean)
  fullData[,-c(1,2)] = scale(fullData[,-c(1,2)], center = X_bar, scale = FALSE)
  fullData3=list(pre=fullData[1:nrow(preData),],post=fullData[(1+nrow(preData)):nrow(fullData),])


  edgeProb11=c()
  edgeProb21=c()
  edgeProb31=c()

  edgeProb12=c()
  edgeProb22=c()
  edgeProb32=c()

  ## lglasso

  aa1=CVlglasso(data = fullData, nlamb =  10,random=random,group=group,K=K)
  dd11=aa1$wi[[1]][upper.tri(aa1$wi[[1]],diag=F)]
  edgeProb11=ifelse(abs(dd11)<=10^(-5),0,1)
  dd12=aa1$wi[[2]][upper.tri(aa1$wi[[2]],diag=F)]
  edgeProb12=ifelse(abs(dd12)<=10^(-5),0,1)




  ## glasso

  ##estiamte the network based on glasso
  s1=cov(preData[,-c(1,2)])
  s2=cov(postData[,-c(1,2)])
  aa21=CVglasso(s=s1,nlam = nlam)$Tuning

  aa22=CVglasso:: CVglasso(s=s2,nlam= nlam)$Tuning


  ## EstimateGroupNetwork

  aa3=tuningBB(data=fullData3,A=phiM,nlam=nlam)
  dd31=aa3[[1]][upper.tri(aa3[[1]],diag=F)]
  edgeProb31=ifelse(abs(dd31)<=10^(-5),0,1)
  dd32=aa3[[2]][upper.tri(aa3[[2]],diag=F)]
  edgeProb32=ifelse(abs(dd32)<=10^(-5),0,1)

  ## true networks
  networkPre=ifelse(abs(sigmaM[[1]])<=10^(-5),0,1)
  edgeProbPre=networkPre[upper.tri(networkPre,diag = F)]
  networkPost=ifelse(abs(sigmaM[[2]])<=10^(-5),0,1)
  edgeProbPost=networkPost[upper.tri(networkPost,diag = F)]

  edgeProb=as.data.frame(rbind(edgeProbPre,edgeProbPost))
  rownames(edgeProb)=c()
  edgeProb1=as.data.frame(rbind(edgeProb11,edgeProb12))
  rownames(edgeProb1)=c()
  edgeProb2=as.data.frame(rbind(edgeProb21,edgeProb22))
  rownames(edgeProb2)=c()
  edgeProb3=as.data.frame(rbind(edgeProb31,edgeProb32))
  rownames(edgeProb3)=c()

  mprediction=cbind(method="true",type=c("pre","post"),edgeProb)
  mprediction1=cbind(method="lglasso",type=c("pre","post"),edgeProb1)
  #colnames(mprediction1)[-c(1,2)]=c()
  mprediction2=cbind(method="glasso",type=c("pre","post"),edgeProb2)
  #colnames(mprediction2)[-c(1,2)]=c()
  mprediction3=cbind(method="jgl",type=c("pre","post"),edgeProb3)
  #colnames(mprediction3)[-c(1,2)]=c()

  prediction=rbind(mprediction,mprediction1,mprediction2,mprediction3)
  return(prediction)
}

power_compare1=function(phiM,sigmaM,timepoints,rho_lglasso,
                        rho_glasso,rho_jgl,zirate=c(0.2,0),tau=c(1,1))
  {
  ## generate the network data




  if (length(rho_lglasso)!=2){
    stop("rho_lglasso should be length 2!")
  }

  if (length(rho_jgl)!=2){
    stop("rho_jgl should be length 2!")
  }

  if (length(rho_glasso)!=2){
    stop("rho_glasso should be length 2!")
  }


    simData=sim_data_old(covmat = sigmaM[3:4],timepoints=timepoints, tau = tau)
    #simData=sim_data(phi1 = phiM[[1]],phi2=phiM[[2]], covmat=sigmaM[3:4],timepoints = timepoints)
    preData=simData$data$pre
    postData=simData$data$post
    fullData=rbind(preData,postData)

    groupIndex=c(rep(1,nrow(preData)),rep(2,nrow(postData)))

    X_bar = apply(fullData[,-c(1,2)], 2, mean)
    fullData[,-c(1,2)] = scale(fullData[,-c(1,2)], center = X_bar, scale = FALSE)
    fullData3=list(pre=fullData[1:nrow(preData),],post=fullData[(1+nrow(preData)):nrow(fullData),])


    edgeProb11=c()
    edgeProb21=c()
    edgeProb31=c()

    edgeProb12=c()
    edgeProb22=c()
    edgeProb32=c()

    ## lglasso

    aa1=lglasso(data = fullData, lambda =  rho_lglasso,type="expFixed",group=groupIndex)
    dd11=aa1$wi[[1]][upper.tri(aa1$wi[[1]],diag=F)]
    edgeProb11=ifelse(abs(dd11)<=10^(-5),0,1)
    dd12=aa1$wi[[2]][upper.tri(aa1$wi[[2]],diag=F)]
    edgeProb12=ifelse(abs(dd12)<=10^(-5),0,1)

    ## glasso


      ##estiamte the network based on glasso
      s1=cov(preData[,-c(1,2)])
      s2=cov(postData[,-c(1,2)])
      aa21=glasso(s=s1,rho= rho_glasso[1])$wi
      dd21=aa21[upper.tri(aa21,diag=F)]
      edgeProb21=ifelse(abs(dd21)<=10^(-5),0,1)

      aa22=glasso(s=s2,rho= rho_glasso[2])$wi
      dd22=aa22[upper.tri(aa22,diag=F)]
      edgeProb22=ifelse(abs(dd22)<=10^(-5),0,1)



    ## jgl

      aa3=BB(data=fullData3,A=phiM,lambda=rho_jgl,type="expFixed")
      dd31=aa3[[1]][upper.tri(aa3[[1]],diag=F)]
      edgeProb31=ifelse(abs(dd31)<=10^(-5),0,1)
      dd32=aa3[[2]][upper.tri(aa3[[2]],diag=F)]
      edgeProb32=ifelse(abs(dd32)<=10^(-5),0,1)

    ## true networks
      networkPre=ifelse(abs(sigmaM[[1]])<=10^(-5),0,1)
      edgeProbPre=networkPre[upper.tri(networkPre,diag = F)]
      networkPost=ifelse(abs(sigmaM[[2]])<=10^(-5),0,1)
      edgeProbPost=networkPost[upper.tri(networkPost,diag = F)]

  edgeProb=as.data.frame(rbind(edgeProbPre,edgeProbPost))
  rownames(edgeProb)=c()
  edgeProb1=as.data.frame(rbind(edgeProb11,edgeProb12))
  rownames(edgeProb1)=c()
  edgeProb2=as.data.frame(rbind(edgeProb21,edgeProb22))
  rownames(edgeProb2)=c()
  edgeProb3=as.data.frame(rbind(edgeProb31,edgeProb32))
  rownames(edgeProb3)=c()

mprediction=cbind(method="true",type=c("pre","post"),edgeProb)
  mprediction1=cbind(method="lglasso",type=c("pre","post"),edgeProb1)
  #colnames(mprediction1)[-c(1,2)]=c()
  mprediction2=cbind(method="glasso",type=c("pre","post"),edgeProb2)
  #colnames(mprediction2)[-c(1,2)]=c()
  mprediction3=cbind(method="jgl",type=c("pre","post"),edgeProb3)
  #colnames(mprediction3)[-c(1,2)]=c()

  prediction=rbind(mprediction,mprediction1,mprediction2,mprediction3)
  return(prediction)
}


ebic_selection=function(phiM,sigmaM,timepoints,rho_lglasso,
                        rho_glasso,rho_jgl,zirate=c(0.2,0),tau=c(1,1))
{
  ## generate the network data




  if (length(rho_lglasso)!=2){
    stop("rho_lglasso should be length 2!")
  }

  if (length(rho_jgl)!=2){
    stop("rho_jgl should be length 2!")
  }

  if (length(rho_glasso)!=2){
    stop("rho_glasso should be length 2!")
  }


  simData=sim_data_old(covmat = sigmaM[3:4],timepoints=timepoints, tau = tau)
  #simData=sim_data(phi1 = phiM[[1]],phi2=phiM[[2]], covmat=sigmaM[3:4],timepoints = timepoints)
  preData=simData$data$pre
  postData=simData$data$post
  fullData=rbind(preData,postData)

  groupIndex=c(rep(1,nrow(preData)),rep(2,nrow(postData)))

  X_bar = apply(fullData[,-c(1,2)], 2, mean)
  fullData[,-c(1,2)] = scale(fullData[,-c(1,2)], center = X_bar, scale = FALSE)
  fullData3=list(pre=fullData[1:nrow(preData),],post=fullData[(1+nrow(preData)):nrow(fullData),])


  edgeProb11=c()
  edgeProb21=c()
  edgeProb31=c()

  edgeProb12=c()
  edgeProb22=c()
  edgeProb32=c()

  ## lglasso

  aa1=lglasso(data = fullData, lambda =  rho_lglasso,type="expFixed",group=groupIndex)
  dd11=aa1$wi[[1]][upper.tri(aa1$wi[[1]],diag=F)]
  edgeProb11=ifelse(abs(dd11)<=10^(-5),0,1)
  dd12=aa1$wi[[2]][upper.tri(aa1$wi[[2]],diag=F)]
  edgeProb12=ifelse(abs(dd12)<=10^(-5),0,1)
  ebic_lglasso=-aa1$ll+(log(nrow(preData))
                            +4*0.5*log(ncol(preData)-2))*sum(edgeProb11)
                        +(log(nrow(postData))
                          +4*0.5*log(ncol(postData)-2))*sum(edgeProb12)



  ## glasso


  ##estiamte the network based on glasso
  s1=cov(preData[,-c(1,2)])
  s2=cov(postData[,-c(1,2)])
  aa21=glasso(s=s1,rho= rho_glasso[1])$wi
  dd21=aa21[upper.tri(aa21,diag=F)]
  edgeProb21=ifelse(abs(dd21)<=10^(-5),0,1)
  ebic_glasso1=-(log(det(aa21))-sum(diag(s1%*%aa21)))
  +(log(nrow(preData))+4*0.5*log(ncol(preData)-2))*sum(edgeProb21)

  aa22=glasso(s=s2,rho= rho_glasso[2])$wi
  dd22=aa22[upper.tri(aa22,diag=F)]
  edgeProb22=ifelse(abs(dd22)<=10^(-5),0,1)
  ebic_glasso2=-(log(det(aa22))-sum(diag(s2%*%aa22)))
  +(log(nrow(postData))+4*0.5*log(ncol(postData)-2))*sum(edgeProb22)
  ebic_glasso=ebic_glasso1+ebic_glasso2


  ## jgl

  aa3=BB(data=fullData3,A=phiM,lambda=rho_jgl)$wiList
  dd31=aa3[[1]][upper.tri(aa3[[1]],diag=F)]
  edgeProb31=ifelse(abs(dd31)<=10^(-5),0,1)
  ebic_jgl1=-(log(det(aa3[[1]]))-sum(diag(s1%*%aa3[[1]])))
  +(log(nrow(preData))+4*0.5*log(ncol(preData)-2))*sum(edgeProb31)

  dd32=aa3[[2]][upper.tri(aa3[[2]],diag=F)]
  edgeProb32=ifelse(abs(dd32)<=10^(-5),0,1)
  ebic_jgl2=-(log(det(aa3[[2]]))-sum(diag(s2%*%aa3[[2]])))
  +(log(nrow(postData))+4*0.5*log(ncol(postData)-2))*sum(edgeProb32)
ebic_jgl=ebic_jgl1+ebic_jgl2

ebic=data.frame(lglasso=ebic_lglasso,glasso=ebic_glasso,jgl=ebic_jgl)
  return(ebic=ebic)
}





ff_glasso=function(index, preData,postData,rho_glasso){
  s1=cov(preData[[index[1]]][,-c(1,2)])
  s2=cov(postData[[index[1]]][,-c(1,2)])
  aa21=glasso(s=s1,rho= rho_glasso[index[2],1])$wi
  dd21=aa21[upper.tri(aa21,diag=F)]
  edgeProb21=ifelse(abs(dd21)<=10^(-5),0,1)

  aa22=glasso(s=s2,rho= rho_glasso[index[2],2])$wi
  dd22=aa22[upper.tri(aa22,diag=F)]
  edgeProb22=ifelse(abs(dd22)<=10^(-5),0,1)
  edgeProb2=as.data.frame(rbind(edgeProb21,edgeProb22))
  rownames(edgeProb2)=c()
  mprediction2=cbind(method="glasso",type=c("pre","post"),edgeProb2)

}


ff_jgl=function(index, preData,postData,phiM,rho_jgl){
  fullData3=list(preData[[index[1]]],postData[[index[1]]])
  aa3=BB(data=fullData3,A=phiM,lambda=rho_jgl[index[2],],type="expFixed")
  dd31=aa3[[1]][upper.tri(aa3[[1]],diag=F)]
  edgeProb31=ifelse(abs(dd31)<=10^(-5),0,1)
  dd32=aa3[[2]][upper.tri(aa3[[2]],diag=F)]
  edgeProb32=ifelse(abs(dd32)<=10^(-5),0,1)


  edgeProb3=as.data.frame(rbind(edgeProb31,edgeProb32))
  rownames(edgeProb3)=c()
  mprediction3=cbind(method="jgl",type=c("pre","post"),edgeProb3)
}


ff_lglasso=function(index, preData,postData,rho_lglasso){
  fullData=rbind(preData[[index[1]]],postData[[index[1]]])
  groupIndex=c(rep(1,nrow(preData[[index[1]]])), rep(2,nrow(postData[[index[1]]])))
  aa1=lglasso(data = fullData, lambda =rho_lglasso[index[2],],
              type="expFixed",group=groupIndex)
  dd11=aa1$wi[[1]][upper.tri(aa1$wi[[1]],diag=F)]
  edgeProb11=ifelse(abs(dd11)<=10^(-5),0,1)
  dd12=aa1$wi[[2]][upper.tri(aa1$wi[[2]],diag=F)]
  edgeProb12=ifelse(abs(dd12)<=10^(-5),0,1)

  edgeProb1=as.data.frame(rbind(edgeProb11,edgeProb12))
  rownames(edgeProb1)=c()
  mprediction1=cbind(method="lglasso",type=c("pre","post"),edgeProb1)
}

ff_true=function(sigmaM){
  ## true networks
  networkPre=ifelse(abs(sigmaM[[1]])<=10^(-5),0,1)
  edgeProbPre=networkPre[upper.tri(networkPre,diag = F)]
  networkPost=ifelse(abs(sigmaM[[2]])<=10^(-5),0,1)
  edgeProbPost=networkPost[upper.tri(networkPost,diag = F)]
  edgeProb=as.data.frame(rbind(edgeProbPre,edgeProbPost))
  rownames(edgeProb)=c()
  mprediction=cbind(method="true",type=c("pre","post"),edgeProb)
}

