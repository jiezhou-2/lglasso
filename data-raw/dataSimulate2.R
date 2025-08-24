
simulate=function(n,p,m1,m2=0,m3,cc,homo=TRUE){
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
    if (homo==TRUE){
    mmlist=list()
    theta = matrix(rnorm(p^2,mean = 0,sd=2), ncol = p,nrow = p)
    theta[lower.tri(theta, diag = TRUE)] = 0
    theta = theta + t(theta) + diag(p)
    theta = theta * real_stru
    theta=MakePositiveDefinite(theta,pd_strategy = "diagonally_dominant",scale = T)$omega
    colnames(theta)=paste0("metabolite",1:p)
    rownames(theta)=paste0("metabolite",1:p)
    sigma=solve(theta)
    fullCovariance= kronecker(cc,sigma)
    }

    if (homo==FALSE){
      for (i in 1:m3) {
        Sigma[[i]]=list()
        Precision[[i]]=list()
        for (j in 1:m3) {
          if (j<i){
            Precision[[i]][[j]]=Precision[[j]][[i]]
            Sigma[[i]][[j]]=Sigma[[j]][[i]]
            next()
          }else{
              ## generate the precision matrices
              theta = matrix(rnorm(p^2,mean = 0,sd=2), ncol = p,nrow = p)
              theta[lower.tri(theta, diag = TRUE)] = 0
              theta = theta + t(theta) + diag(p)
              ### apply the required sparsity
              theta = theta * real_stru

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
    colnames(prior_stru)=paste0("metabolite",1:p)
    rownames(prior_stru)=paste0("metabolite",1:p)
    Structure[[i]][[j]]=prior_stru
    ## generate the precision matrices
    theta = matrix(rnorm(p^2,mean = 0,sd=2), ncol = p,nrow = p)
    theta[lower.tri(theta, diag = TRUE)] = 0
    theta = theta + t(theta) + diag(p)
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









