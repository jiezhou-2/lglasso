library(fake)
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



sim_data=function(n,tt=5,tau=c(1,1),covmat){
  ## true structure
  timepoint1=vector("list",n)
  timepoint2=vector("list",n)
  p=nrow(covmat[[1]])
  for (i in 1:n) {
    m3=sample(x=1:tt[1],1,prob = rep(1,1,tt[1]))
    m4=sample(x=1:tt[2],1,prob = rep(1,1,tt[2]))
      t1=stats::rexp(m3,rate=2)
      t2=stats::rexp(m4,rate=2)
      timepoint1[[i]]=cumsum(t1)[1:m3]
      timepoint2[[i]]=cumsum(t2)[1:m4]

  }
  cc1=lapply(timepoint1,phifunction, tau=tau[1])
  cc2=lapply(timepoint2,phifunction, tau=tau[2])

  fullCovariance1= lapply(cc1,function(x,cc) {kronecker(cc,x)},x=covmat[[1]])
    fullCovariance2= lapply(cc2,function(x,cc) {kronecker(cc,x)},x=covmat[[2]])

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
    return(list(data=list(pre=a1,post=a2)))
}



BB=function(data,lambda,diagonal=TRUE,lower=c(0.01,0.1),upper=c(10,5)){
    obj=0
    aa=0
    bb=0
    B=vector("list",2)
    p= p=ncol(data[[1]])-2
    # Create a mask matrix
    mask1 <- matrix(lambda[1], p, p)
    diag(mask1) <- 0
    mask2 <- matrix(lambda[2], p, p)
    diag(mask2) <- 0
    B[[1]]=Variable(p,p,PSD=TRUE) # tissue wise inverse correlation matrix
    B[[2]]=Variable(p,p,PSD=TRUE)
    for (i in 1:2){
      dd=data[[i]]
      nn=length(unique(dd[,1]))
      p=ncol(dd)-2
      data_sub=split(dd[,-c(1,2)],factor(dd[,1],unique(dd[,1])))
      amatrix=0

      for (j in 1:nn) {
        xx=as.matrix(data_sub[[j]])
        amatrix=t(xx)%*%xx+amatrix
      }
      obj=log_det(B[[i]])-matrix_trace(B[[i]]%*%amatrix)/nrow(dd)+obj
      aa=aa+ sum(abs(B[[i]])*mask1)
    }

    bb=bb+sum(abs((B[[1]]-B[[2]]))*mask2)
    obj=-obj+aa+bb
      prob=Problem(Minimize(obj))
      result=CVXR::solve(prob)
      S_est= lapply(B, function(x) result$getValue(x))
      return(wi=S_est)
}
