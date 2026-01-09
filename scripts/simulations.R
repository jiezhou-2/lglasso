




simulate_general=function(n,p,m1,m2=0,m3,cc){
  ## true structure
  if (ncol(cc)!=m3){stop("Unmatched inputs!")}
  real_stru=matrix(0, nrow = p, ncol = p)
  real_stru[lower.tri(real_stru,diag = TRUE)]=1
  index=which(real_stru==0,arr.ind = TRUE)
  a=sample(1:nrow(index),m1, replace = F)
  real_stru[index[a,]]=1
  real_stru[lower.tri(real_stru,diag = TRUE)]=0
  real_stru=real_stru+t(real_stru)+diag(p)
  mu=rep(0,p*m3)
  Precision=list()
  Structure=list()
  Sigma=list()
  Data=list()
  ## tissue specific structure

  if (m2==0){
          ## generate the precision matrices
          theta = matrix(stats::rnorm(p^2,mean = 0,sd=2), ncol = p,nrow = p)
          theta = theta + t(theta)
          diag(theta)=1
          ### apply the required sparsity
          theta = theta * real_stru

          # force it to be positive definite
          theta=MakePositiveDefinite(theta,pd_strategy = "diagonally_dominant",scale = TRUE)$omega
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
          distrubance[lower.tri(distrubance,diag = TRUE)]=1
          index=which(distrubance==0,arr.ind = TRUE)
          a=sample(1:nrow(index),m2, replace = F)
          distrubance[index[a,]]=1
          distrubance[lower.tri(distrubance,diag = TRUE)]=0
          distrubance=distrubance+t(distrubance)+diag(p)
          prior_stru=(real_stru+distrubance)%%2

          diag(prior_stru)=1
          #colnames(prior_stru)=paste0("metabolite",1:p)
          #rownames(prior_stru)=paste0("metabolite",1:p)
          Structure[[i]][[j]]=prior_stru
          ## generate the precision matrices
          theta = matrix(stats::rnorm(p^2,mean = 0,sd=2), ncol = p,nrow = p)
          #theta[lower.tri(theta, diag = TRUE)] = 0
          theta = theta + t(theta)
          diag(theta)=1
          ### apply the required sparsity
          theta = theta * prior_stru

          # force it to be positive definite
          theta=MakePositiveDefinite(theta,pd_strategy = "diagonally_dominant",scale = TRUE)$omega
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
    fullCovariance=MakePositiveDefinite(fullCovariance,pd_strategy = "diagonally_dominant",scale = TRUE)$omega

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







simulate_long=function(n,p,m1,tt=5,m2=0,tau){
  ## true structure
  timepoint1=vector("list",n)
  timepoint2=vector("list",n)

  for (i in 1:n) {
    m3=sample(x=1:tt,1,prob = rep(1,1,tt))
    if (length(tau)==1){
    t1=stats::rexp(m3)
    timepoint1[[i]]=cumsum(t1)[1:m3]
    }else{
      t1=stats::rexp(2*m3)
      timepoint1[[i]]=cumsum(t1)[1:m3]
      timepoint2[[i]]=cumsum(t1)[(m3+1):(2*m3)]
    }
  }
  cc1=lapply(timepoint1,phifunction, tau=tau[1])
  if (length(tau)==2){
  cc2=lapply(timepoint2,phifunction, tau=tau[2])
  }
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


  Precision=list()
  Sigma=list()
  mmlist=list()
  theta = matrix(stats::rnorm(p^2,mean = 0,sd=2), ncol = p,nrow = p)
  theta[lower.tri(theta, diag = TRUE)] = 0
  theta = theta + t(theta) + diag(p)
  theta1 = theta * real_stru1
  theta2 = theta * real_stru2
  theta1=MakePositiveDefinite(theta1,pd_strategy = "diagonally_dominant",scale = TRUE)$omega
  theta2=MakePositiveDefinite(theta2,pd_strategy = "diagonally_dominant",scale = TRUE)$omega
  #colnames(theta)=paste0("metabolite",1:p)
  #rownames(theta)=paste0("metabolite",1:p)
  sigma1=solve(theta1)
  sigma2=solve(theta2)
  fullCovariance1= lapply(cc1,function(x,cc) {kronecker(cc,x)},x=sigma1)
  if (length(tau)==2){
    fullCovariance2= lapply(cc2,function(x,cc) {kronecker(cc,x)},x=sigma2)
  }
  fulldata=c()
  a1=c()
  a2=c()
  for (i in 1:n) {
    ai=c()
    m3=nrow(fullCovariance1[[i]])
    mu=rep(0,m3)
    data1=MASS::mvrnorm(1,mu=mu,Sigma = fullCovariance1[[i]])
    for (j in 1:(m3/p)) {
      ai=as.data.frame(rbind(ai,data1[c(((j-1)*p+1) :(j*p))]))
    }
    subject=paste0("subject",i)
    ai=cbind(subject,timepoint1[[i]],ai)
    a1=rbind(a1,ai)
  }
  colnames(a1)[2]="time"

  if (length(tau)==2){
    for (i in 1:n) {
    ai=c()
    m3=nrow(fullCovariance2[[i]])
    mu=rep(0,m3)
    data2=MASS::mvrnorm(1,mu=mu,Sigma = fullCovariance2[[i]])
    for (j in 1:(m3/p)) {
      ai=as.data.frame(rbind(ai,data2[c(((j-1)*p+1) :(j*p))]))
    }

    subject=paste0("subject",i)
    ai=cbind(subject,timepoint2[[i]],ai)
    a2=rbind(a2,ai)
    }
    colnames(a2)[2]="time"
    return(list(data=list(pre=a1,post=a2),network=list(pre=real_stru1,post=real_stru2)))
  }

  return(list(data=a1,network=real_stru1))
}


simulate_randomTau=function(n,p,m1,tt,m2,alpha,group){
  ## true structure
  timepoint1=vector("list",n)
  timepoint2=vector("list",n)
  cc1=vector("list",n)
  cc2=vector("list",n)
trueTau=rexp(n=n,rate=alpha)
  for (i in 1:n) {
    m3=sample(x=1:tt,1,prob = rep(1,1,tt))
    if (group==1){
      t1=stats::rexp(m3)
      timepoint1[[i]]=cumsum(t1)[1:m3]
      cc1[[i]]=phifunction(t=timepoint1[[i]],tau=trueTau[i])
    }else{
      t1=stats::rexp(2*m3)
      timepoint1[[i]]=cumsum(t1)[1:m3]
      timepoint2[[i]]=cumsum(t1)[(m3+1):(2*m3)]
      cc1[[i]]=phifunction(t=timepoint1[[i]],tau=trueTau[i])
      cc2[[i]]=phifunction(t=timepoint2[[i]],tau=trueTau[i])
    }
  }



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


  Precision=list()
  Sigma=list()



  mmlist=list()
  theta = matrix(stats::rnorm(p^2,mean = 0,sd=2), ncol = p,nrow = p)
  theta[lower.tri(theta, diag = TRUE)] = 0
  theta = theta + t(theta) + diag(p)
  theta1 = theta * real_stru1
  theta2 = theta * real_stru2
  theta1=MakePositiveDefinite(theta1,pd_strategy = "diagonally_dominant",scale = TRUE)$omega
  theta2=MakePositiveDefinite(theta2,pd_strategy = "diagonally_dominant",scale = TRUE)$omega
  #colnames(theta)=paste0("metabolite",1:p)
  #rownames(theta)=paste0("metabolite",1:p)
  sigma1=solve(theta1)
  sigma2=solve(theta2)
  fullCovariance1= lapply(cc1,function(x,cc) {kronecker(cc,x)},x=sigma1)
  if (group==2){
    fullCovariance2= lapply(cc2,function(x,cc) {kronecker(cc,x)},x=sigma2)
  }
  fulldata=c()
  a1=c()
  a2=c()
  for (i in 1:n) {
    ai=c()
    m3=nrow(fullCovariance1[[i]])
    mu=rep(0,m3)
    data1=MASS::mvrnorm(1,mu=mu,Sigma = fullCovariance1[[i]])
    for (j in 1:(m3/p)) {
      ai=as.data.frame(rbind(ai,data1[c(((j-1)*p+1) :(j*p))]))
    }
    subject=paste0("subject",i)
    ai=cbind(subject,timepoint1[[i]],ai)
    a1=rbind(a1,ai)
  }
  colnames(a1)[2]="time"

  if (group==2){
    for (i in 1:n) {
      ai=c()
      m3=nrow(fullCovariance2[[i]])
      mu=rep(0,m3)
      data2=MASS::mvrnorm(1,mu=mu,Sigma = fullCovariance2[[i]])
      for (j in 1:(m3/p)) {
        ai=as.data.frame(rbind(ai,data2[c(((j-1)*p+1) :(j*p))]))
      }

      subject=paste0("subject",i)
      ai=cbind(subject,timepoint2[[i]],ai)
      a2=rbind(a2,ai)
    }
    colnames(a2)[2]="time"
    return(list(data=list(pre=a1,post=a2),network=list(pre=real_stru1,post=real_stru2),tau=trueTau))
  }

  return(list(data=list(pre=a1),network=list(pre=real_stru1),tau=trueTau))
}



simulate_heter=function(n,p,m1,alpha){
  ## true structure
  timepoint=vector("list",n)
  tau=c()
  cc=vector("list",)
  for (i in 1:n) {
    m3=sample(x=1:15,1,prob = rep(1,1,length=15))
    t1=stats::rexp(m3)
    timepoint[[i]]=cumsum(t1)
    tau=c(tau,stats::rexp(1,rate=alpha))
  }
  for (j in 1:length(tau)) {
    cc[[j]]=phifunction(t=timepoint[[j]],tau=c(tau[j],1))
  }

  real_stru=matrix(0, nrow = p, ncol = p)
  real_stru[lower.tri(real_stru,diag = TRUE)]=1
  index=which(real_stru==0,arr.ind = TRUE)
  a=sample(1:nrow(index),m1, replace = F)
  real_stru[index[a,]]=1
  real_stru[lower.tri(real_stru,diag = TRUE)]=0
  real_stru=real_stru+t(real_stru)+diag(p)
  Precision=list()
  Sigma=list()



  mmlist=list()
  theta = matrix(stats::rnorm(p^2,mean = 0,sd=2), ncol = p,nrow = p)
  theta = theta + t(theta)
  diag(theta)=1
  theta = theta * real_stru
  theta=MakePositiveDefinite(theta,pd_strategy = "diagonally_dominant",scale = TRUE)$omega
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



Simulate=function(type=c("general","longihomo","longiheter"),n=20,p=20,m1=20,m2=1,m3=3,tt=5,cc=diag(m3),tau=c(2,1),alpha=2,group){
  type=match.arg(type)

  if (type=="general"){
    data=simulate_general(n=n,p=p,m1=m1,m2=m2,m3=m3,cc)
  }

  if (type=="longihomo"){
    data=simulate_long(n=n,p=p,m1=m1,m2=m2,tau=tau,tt=tt)
  }

  if (type=="longiheter"){

    data=simulate_randomTau(n=n, p=p,m1=m1,m2=m2,tt=tt, alpha=alpha,group = group)
  }

  return(data)
}



power_compare1=function(sigmaM,timepoints,rho_lglasso,
                        rho_glasso,rho_jgl,tau,tt)
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

  if (length(tau)!=length(timepoints[[1]])){
    stop("tau should have the same length as timepoints!")
  }

  #simData=sim_data_old(covmat = sigmaM[3:4],timepoints=timepoints, tau = tau)
  simData=Simulate(type="longihomo",n=n,p=p,m1=m1,m2=m2,tt=tt,tau=tau,group=2)
  preData=simData$data$pre
  postData=simData$data$post
  fullData=rbind(preData,postData)
  groupIndex=c(rep(1,nrow(preData)),rep(2,nrow(postData)))


  X_bar = apply(fullData[,-c(1,2)], 2, mean)
  fullData[,-c(1,2)] = scale(fullData[,-c(1,2)], center = X_bar, scale = FALSE)
  fullData3=list(pre=fullData[1:nrow(preData),],post=fullData[(1+nrow(preData)):nrow(fullData),])
  n=length(unique(fullData[,1]))

  edgeProb11=c()
  edgeProb21=c()
  edgeProb31=c()

  edgeProb12=c()
  edgeProb22=c()
  edgeProb32=c()

  ## lglasso
  #browser()
  aa1=lglasso(data = fullData, lambda =  rho_lglasso,type="expFixed",group=groupIndex,random=random)
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



  ## EstimateGroupNetwork

  aa3=BBB(data=fullData3,lambda=rho_jgl,type="expFixed")
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
  tau=data.frame(matrix(c(aa1$tau,rep(0,ncol(edgeProb)-length(aa1$tau))),nrow=1))
  mtau=cbind(method="estiamte",type="tau",tau)
  colnames(mtau)=colnames(mprediction)
  prediction=rbind(mprediction,mprediction1,mprediction2,mprediction3,mtau)
  return(prediction=prediction)
}



power_compare_heter=function(n,p,m1,m2,tt,tau,group,rho_lglasso,
                             rho_glasso,rho_jgl)
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



  simData=Simulate(type="longiheter",n=n,p=p,m1=m1,m2=m2,tt=5,tau=tau,group=2)
  preData=simData$data$pre
  postData=simData$data$post
  fullData=rbind(preData,postData)
  groupIndex=c(rep(1,nrow(preData)),rep(2,nrow(postData)))


  X_bar = apply(fullData[,-c(1,2)], 2, mean)
  fullData[,-c(1,2)] = scale(fullData[,-c(1,2)], center = X_bar, scale = FALSE)
  fullData3=list(pre=fullData[1:nrow(preData),],post=fullData[(1+nrow(preData)):nrow(fullData),])
  n=length(unique(fullData[,1]))

  edgeProb11=c()
  edgeProb21=c()
  edgeProb31=c()

  edgeProb12=c()
  edgeProb22=c()
  edgeProb32=c()

  ## lglasso
  #browser()
  aa1=lglasso(data = fullData, lambda =  rho_lglasso,type="expFixed",group=groupIndex,random=TRUE)
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



  ## EstimateGroupNetwork

  aa3=BBB(data=fullData3,lambda=rho_jgl,type="expFixed")
  dd31=aa3[[1]][upper.tri(aa3[[1]],diag=F)]
  edgeProb31=ifelse(abs(dd31)<=10^(-5),0,1)
  dd32=aa3[[2]][upper.tri(aa3[[2]],diag=F)]
  edgeProb32=ifelse(abs(dd32)<=10^(-5),0,1)

  ## true networks
  networkPre=ifelse(abs(simData$network$pre)<=10^(-5),0,1)
  edgeProbPre=networkPre[upper.tri(networkPre,diag = F)]
  networkPost=ifelse(abs(simData$network$post)<=10^(-5),0,1)
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
  tau=data.frame(matrix(c(simData$tau,rep(0,ncol(edgeProb)-length(simData$tau))),nrow=1))
  mtau=cbind(method="estiamte",type="tau",tau)
  colnames(mtau)=colnames(mprediction)
  prediction=rbind(mprediction,mprediction1,mprediction2,mprediction3,mtau)
  return(prediction=prediction)
}

