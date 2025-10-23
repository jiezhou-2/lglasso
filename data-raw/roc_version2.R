
sim_homo=function(p,K,tau,age,zirate=c(0.2,0)){
  print("homo data are generated")
  m=length(age)
  data=vector("list",m)
  Sigma=solve(K)
  sqK=chol(Sigma)
  for (i in 1:m) {
    n=length(age[[i]])
    a=matrix(ncol = p,nrow = n)
    error1=matrix(rnorm(p*n),nrow = p)
    error2=t(t(sqK)%*%error1)
    a[1,]=error2[1,]
    if (n>1){
    for (t in 2:n) {
      coe=exp(-tau*abs(age[[i]][t]-age[[i]][t-1]))
      a[t,]=a[t-1,]*coe+error2[t,]*sqrt(1-coe^2)
    }
    }
    dd=cenfunction(a,zirate = zirate)
    data[[i]]=cbind(i,age[[i]],dd)
    colnames(data[[i]])[1:2]=c("id","age")
  }
  result=list(data=data,tau=tau,precision=K)
  return(result)
}

#' Data censor
#'
#' @param data Data matrix to be censored
#' @param zirate Censor parameter, in which
#' zirate[1] is the lower percentage of data to be censored, zirate[2] is
#' the censor probability for the selected ones.
#' @return A censored matrix
#' @export


cenfunction=function(data,zirate=c(0.2,0)){
  if (prod(zirate)==0){
    return(dd=data)
  }
  n=nrow(data)
  p=ncol(data)-2
  M=data[,-c(1,2)]
  for (j in 1:p) {
    qq=quantile(M[,j],probs=zirate[1])
    index=which(M[,j]<=qq)
    cenmark=sample(x=c(1,0),size=length(index),replace = T,prob = c(zirate[2],1-zirate[2]))
    M[index,j]=ifelse(cenmark==1,qq,M[index,j])
  }
  dd=cbind(data[,c(1,2)],M)
  return(dd)
}


power_compare1=function(n,tt=5,p,m1,m2,tau,l,rho,zirate=c(0.2,0)){


mprediction1=c()
mprediction2=c()
mprediction3=c()

## generate the network data

sigmaM=sim_stru(p=p,m1=m1,m2=m2)




  for (h in 1:l){
    print(paste0("The ",h,"th"," simulation:"))

    simData=sim_data(n=n,tt=tt,tau=tau,covmat=sigmaM[3:4])
    preData=simData$data$pre
    postData=simData$data$post
    fullData=rbind(rbind(preData,postData))
    fullData3=list(pre=preData,post=postData)
    groupIndex=c(rep(1,nrow(preData)),rep(2,nrow(postData)))
    #covariate=cbind(id,x)
    #Tem=log(p)/(2*log(1/prob-1))
    #Tem=2

    aa1=vector("list",length(rho[[1]]))
    aa21=vector("list",length(rho[[2]]))
    aa22=vector("list",length(rho[[2]]))
    aa3=vector("list",length(rho[[3]]))


    edgeProb11=c()
    edgeProb21=c()
    edgeProb31=c()

    edgeProb12=c()
    edgeProb22=c()
    edgeProb32=c()

    ## lglasso

    for (j in 1:ncol(rho[[1]])) {
      #print(paste0("tuning parameter is ", rho[[1]][j]))
      #browser()
      aa1[[j]]=lglasso(data = fullData, lambda =  rho[[1]][,j],type="expFixed",group=groupIndex)
      dd11=aa1[[j]]$wi[[1]][upper.tri(aa1[[j]]$wi[[1]],diag=F)]
edgeProb11=rbind(edgeProb11,ifelse(abs(dd11)<=10^(-5),0,1))
dd12=aa1[[j]]$wi[[2]][upper.tri(aa1[[j]]$wi[[2]],diag=F)]
edgeProb12=rbind(edgeProb12,ifelse(abs(dd12)<=10^(-5),0,1))
#results[[1]][[j]]=as.numeric(comparison(graph,aa1[[j]]$omega))
    }



    ## glasso

    for (j in 1:ncol(rho[[2]])) {
      ##estiamte the network based on glasso
      s1=cov(preData[,-c(1,2)])
      s2=cov(postData[,-c(1,2)])
      aa21[[j]]=glasso(s=s1,rho= rho[[2]][1,j])$wi
      dd21=aa21[[j]][upper.tri(aa21[[j]],diag=F)]
      edgeProb21=rbind(edgeProb21,ifelse(abs(dd21)<=10^(-5),0,1))

      aa22[[j]]=glasso(s=s2,rho= rho[[2]][2,j])$wi
      dd22=aa22[[j]][upper.tri(aa22[[j]],diag=F)]
      edgeProb22=rbind(edgeProb22,ifelse(abs(dd22)<=10^(-5),0,1))
      #results[[2]][[j]]=as.numeric(comparison(graph,aa2[[j]]))
    }


    ## EstimateGroupNetwork
    for (j in 1:ncol(rho[[3]])) {
      #browser()
      aa3[[j]]=BB(data=fullData3,lambda=rho[[3]][,j])
      dd31=aa3[[j]][[1]][upper.tri(aa3[[j]][[1]],diag=F)]
      edgeProb31=rbind(edgeProb31,ifelse(abs(dd31)<=10^(-5),0,1))
      dd32=aa3[[j]][[2]][upper.tri(aa3[[j]][[2]],diag=F)]
      edgeProb32=rbind(edgeProb32,ifelse(abs(dd32)<=10^(-5),0,1))
      #results[[3]][[j]]=as.numeric(comparison(graph,aa3[[j]]))
    }

  }
edgeProb1=rbind(edgeProb11,edgeProb12)
edgeProb2=rbind(edgeProb21,edgeProb22)
edgeProb3=rbind(edgeProb31,edgeProb32)


mprediction1=colSums(edgeProb1)/l
mprediction2=colSums(edgeProb2)/l
mprediction3=colSums(edgeProb3)/l

prediction=data.frame(net_lglasso=mprediction1,net_glasso=mprediction2, net_jgl=mprediction3)
  return(prediction)
}




power_compare2=function(m,nn=15,p,coe,l,rho,K,heter,community2=F,uu=c(0,0),zirate=c(0.2,0)){
  results=vector("list",length = 5) # container for the final FPR and TPR
  results[[1]]=vector("list",length = length(rho[[1]]))
  results[[2]]=vector("list",length = length(rho[[2]]))
  results[[3]]=vector("list",length = length(rho[[3]]))

  RR=vector("list",length = 5)
  RR[[1]]=matrix(nrow=l,ncol = 2)
  RR[[2]]=matrix(nrow=l,ncol = 2)
  RR[[3]]=matrix(nrow=l,ncol = 2)


  for (i in 1:5) {
    for (j in 1:length(rho[[i]])) {
      results[[i]][[j]]= matrix(nrow = l,ncol = 2)
    }
  }
  GG=vector("list",l)
  mprediction1=c()
  mprediction2=c()
  mprediction3=c()
  for (h in 1:l){
    print(paste0("The ",h,"th"," simulation:"))
    # subject level covariates
    x1=sample(x=c(0,1),size=m,prob=c(0.5,0.5),replace = T)
    x2=runif(m,min = 0,max = 1)
    if (community2==T){
      alpha=c(exp(uu[1]),exp(uu[2]))
    }else{
      alpha=exp(coe[1]+coe[2]*x1+coe[3]*x2)
    }
    #x=cbind(x1,x2)
    ## generate the observation time for each subjects
    # browser()
    age=vector("list",m)
    for (k in 1:m) {
      #browser()
      n=rpois(1,nn)+2
      a1=rpois(n=1,lambda = 1) # the space between observations

      age[[k]][1]=max(a1,0.5)
      if (n>1) {
        for (i in 2:n) {
          ai=rpois(n=1,lambda = 1)
          age[[k]][i]=age[[k]][i-1]+max(ai,0.5)
          # if (length(age[[k]])!= length(unique(age[[k]]))){
          #   browser()
          #   }
        }
      }
    }

    ## generate the network data

    if (heter==TRUE & community2==F){
      ss=sim_heter(p = p,K=K,alpha = alpha,age = age,zirate=zirate)
    }
    if (heter==TRUE & community2==T){
      ss=sim_2heter(p=p,K1=K[[1]],K2=K[[2]],alpha1=exp(uu[1]),alpha2=exp(uu[2]),age=age,zirate = zirate)
    }
    if (heter==F & community2==F){
      ss=sim_homo(p = p,K=K,tau = 1/alpha,age = age,zirate=zirate)
    }

    if (heter==F & community2==T){
      ss=sim_2homo(p=p,K1=K[[1]],K2=K[[2]],tau1=1/alpha[1],tau2=1/alpha[2],age=age,zirate = zirate)
    }

    simdata=ss$data
    lower=0.01
    upper=20
    graph=ss$precision
    dd=do.call(rbind,simdata)
    id=unique(dd[,1])
    #covariate=cbind(id,x)
    #Tem=log(p)/(2*log(1/prob-1))
    Tem=2

    aa1=vector("list",length(rho[[1]]))
    aa2=vector("list",length(rho[[2]]))
    aa3=vector("list",length(rho[[3]]))
    aa4=vector("list",length(rho[[4]]))
    aa5=vector("list",length(rho[[5]]))
    bic1=c()
    bic2=c()
    bic3=c()
    bic4=c()
    bic5=c()
    edgeProb1=c()
    edgeProb2=c()
    edgeProb3=c()
    edgeProb4=c()
    edgeProb5=c()
    ## lglasso

    for (j in 1:length(rho[[1]])) {
      print(paste0("tuning parameter is ", rho[[1]][j]))
      aa1[[j]]=lglasso(data = dd, rho =0.5*rho[[1]][j],heter = heter)
      # if (heter==TRUE){
      #     aa[[j]]=lglasso(data = dd, rho =0.5*rho[[1]][j])$omega
      # }else{
      #   aa[[j]]=lglasso(data = dd, rho = 0.5*rho[[1]][j],heter = F)$omega
      # }
      #browser()
      dd1=aa1[[j]]$omega[upper.tri(aa1[[j]]$omega,diag=F)]
      edgeProb1=rbind(edgeProb1,ifelse(abs(dd1)<=10^(-5),0,1))
      bb1=-2*aa1[[j]]$ll +  0.5*length(which(aa1[[j]]$omega!=0))*log(nrow(dd))
      +0.5*length(which(aa1[[j]]$omega!=0))*log(p)/Tem
      bic1=c(bic1,bb1)
      results[[1]][[j]][h,]=as.numeric(comparison(graph,aa1[[j]]$omega))
    }

    mprediction1=rbind(mprediction1,colSums(edgeProb1)/nrow(edgeProb1))
    ## glasso

    for (j in 1:length(rho[[2]])) {
      ##estiamte the network based on glasso
      s=cov(dd[,-c(1,2)])
      aa2[[j]]=glasso(s=s,rho=0.5*rho[[2]][j])$wi
      dd2=aa2[[j]][upper.tri(aa2[[j]],diag=F)]
      edgeProb2=rbind(edgeProb2,ifelse(abs(dd2)<=10^(-5),0,1))
      bb2=bicfunction(data=dd,G=as.matrix(aa2[[j]]),Tem=2)
      bic2=c(bic2,bb2)
      results[[2]][[j]][h,]=as.numeric(comparison(graph,aa2[[j]]))
    }

    mprediction2=rbind(mprediction2,colSums(edgeProb2)/nrow(edgeProb2))
    ## nh

    for (j in 1:length(rho[[3]])) {
      prior=addition(data=dd,lambda=0.5*rho[[3]][j])
      aa3[[j]]=mle_net(data=dd,prior=prior)
      dd3=aa3[[j]][upper.tri(aa3[[j]],diag=F)]
      edgeProb3=rbind(edgeProb3,ifelse(abs(dd3)<=10^(-5),0,1))
      bb3=bicfunction(data=dd,G=aa3[[j]],Tem=2)
      bic3=c(bic3,bb3)
      results[[3]][[j]][h,]=as.numeric(comparison(graph,aa3[[j]]))

    }
    mprediction3=rbind(mprediction3,colSums(edgeProb3)/nrow(edgeProb3))
  }
}



power_compare3=function(n,tt=5,p,m1,m2,tau1,tau2,rho11,rho12,rho21,rho22,rho3,zirate1=0.2,zirate2=0){


  mprediction1=c()
  mprediction2=c()
  mprediction3=c()

  ## generate the network data

  sigmaM=sim_stru(p=p,m1=m1,m2=m2)




  for (h in 1:l){
    print(paste0("The ",h,"th"," simulation:"))

    simData=sim_data(n=n,tt=tt,tau=tau,covmat=sigmaM[3:4])
    preData=simData$data$pre
    postData=simData$data$post
    fullData=rbind(rbind(preData,postData))
    fullData3=list(pre=preData,post=postData)
    groupIndex=c(rep(1,nrow(preData)),rep(2,nrow(postData)))
    #covariate=cbind(id,x)
    #Tem=log(p)/(2*log(1/prob-1))
    #Tem=2

    aa1=vector("list",length(rho[[1]]))
    aa21=vector("list",length(rho[[2]]))
    aa22=vector("list",length(rho[[2]]))
    aa3=vector("list",length(rho[[3]]))


    edgeProb11=c()
    edgeProb21=c()
    edgeProb31=c()

    edgeProb12=c()
    edgeProb22=c()
    edgeProb32=c()

    ## lglasso

    for (j in 1:ncol(rho[[1]])) {
      #print(paste0("tuning parameter is ", rho[[1]][j]))
      #browser()
      aa1[[j]]=lglasso(data = fullData, lambda =  rho[[1]][,j],type="expFixed",group=groupIndex)
      dd11=aa1[[j]]$wi[[1]][upper.tri(aa1[[j]]$wi[[1]],diag=F)]
      edgeProb11=rbind(edgeProb11,ifelse(abs(dd11)<=10^(-5),0,1))
      dd12=aa1[[j]]$wi[[2]][upper.tri(aa1[[j]]$wi[[2]],diag=F)]
      edgeProb12=rbind(edgeProb12,ifelse(abs(dd12)<=10^(-5),0,1))
      #results[[1]][[j]]=as.numeric(comparison(graph,aa1[[j]]$omega))
    }



    ## glasso

    for (j in 1:ncol(rho[[2]])) {
      ##estiamte the network based on glasso
      s1=cov(preData[,-c(1,2)])
      s2=cov(postData[,-c(1,2)])
      aa21[[j]]=glasso(s=s1,rho= rho[[2]][1,j])$wi
      dd21=aa21[[j]][upper.tri(aa21[[j]],diag=F)]
      edgeProb21=rbind(edgeProb21,ifelse(abs(dd21)<=10^(-5),0,1))

      aa22[[j]]=glasso(s=s2,rho= rho[[2]][2,j])$wi
      dd22=aa22[[j]][upper.tri(aa22[[j]],diag=F)]
      edgeProb22=rbind(edgeProb22,ifelse(abs(dd22)<=10^(-5),0,1))
      #results[[2]][[j]]=as.numeric(comparison(graph,aa2[[j]]))
    }


    ## EstimateGroupNetwork
    for (j in 1:ncol(rho[[3]])) {
      #browser()
      aa3[[j]]=BB(data=fullData3,lambda=rho[[3]][,j])
      dd31=aa3[[j]][[1]][upper.tri(aa3[[j]][[1]],diag=F)]
      edgeProb31=rbind(edgeProb31,ifelse(abs(dd31)<=10^(-5),0,1))
      dd32=aa3[[j]][[2]][upper.tri(aa3[[j]][[2]],diag=F)]
      edgeProb32=rbind(edgeProb32,ifelse(abs(dd32)<=10^(-5),0,1))
      #results[[3]][[j]]=as.numeric(comparison(graph,aa3[[j]]))
    }

  }
  edgeProb1=rbind(edgeProb11,edgeProb12)
  edgeProb2=rbind(edgeProb21,edgeProb22)
  edgeProb3=rbind(edgeProb31,edgeProb32)


  mprediction1=colSums(edgeProb1)/l
  mprediction2=colSums(edgeProb2)/l
  mprediction3=colSums(edgeProb3)/l

  prediction=data.frame(net_lglasso=mprediction1,net_glasso=mprediction2, net_jgl=mprediction3)
  return(prediction)
}

