

#' Main simulation function comparing five methods for network selection.
#'
#' @param m number of subjects
#' @param n  number of observations per subjects
#' @param p the dimension of the data to be generated
#' @param coe coefficients for covariates
#' @param l the simulation scale
#' @param rho1 tuning parameter for glasso,nh and lglasso
#' @param rho2 tuning parameter for
#' @return list with length equal to 3.
#' @export

power_compare2=function(m,n,p,coe,l,rho,prob,heter,nu){
  results=vector("list",length = 5) # container for the final FPR and TPR
  results[[1]]=vector("list",length = length(rho[[1]]))
  results[[2]]=vector("list",length = length(rho[[2]]))
  results[[3]]=vector("list",length = length(rho[[3]]))
  results[[4]]=vector("list",length = length(rho[[4]]))
  results[[5]]=vector("list",length = length(rho[[5]]))
  age=vector("list",m) #generate the time points container



  bb=vector("list",length = 5)
  bb[[1]]=matrix(nrow = l,ncol=2)
  bb[[2]]=matrix(nrow = l,ncol=2)
  bb[[3]]=matrix(nrow = l,ncol=2)
  bb[[4]]=matrix(nrow = l,ncol=2)
  bb[[5]]=matrix(nrow = l,ncol=2)
  names(bb)=c("lglasso","glasso","nh","co","la")
  for (h in 1:l){
    print(paste0("The ",h,"th"," simulation:"))
    # subject level covariates
    x1=sample(x=c(0,1),size=m,prob=c(0.5,0.5),replace = T)
    x2=runif(m,min = 0,max = 1)
    alpha=exp(coe[1]+coe[2]*x1+coe[3]*x2)
    x=cbind(x1,x2)
    ## generate the observation time for each subjects
    for (k in 1:m) {
      a1=rpois(n=n,lambda = 1) # the space between observations
      age[[k]][1]=max(a1[1],0.5)
      for (i in 2:n) {
        age[[k]][i]=age[[k]][i-1]+max(a1[i-1],0.5)
      }
    }

    ## generate the network data
    if (heter==TRUE){
      ss=sim_heter(p = p,prob=prob,alpha = alpha,age = age)
    }else{
      ss=sim_homo(p = p,prob=prob,tau = 1/alpha[1],age = age)
    }
    simdata=ss$data
    tau=ss$tau
    lower=min(abs(tau))
    upper=max(abs(tau))
    graph=ss$precision
    dd=do.call(rbind,simdata)
    id=unique(dd[,1])
    covariates=cbind(id,x)

bic1=c()
bic2=c()
bic3=c()
bic4=c()
bic5=c()
##estiamte the network based on lglasso
    for (j in 1:length(rho[[1]])) {
      if (heter==TRUE){
        if (length(which(coe==0))==2){
          results[[1]][[j]]=lglasso(data = dd, rho = 0.5*rho[[1]][j],heter=T)
        }else{
          results[[1]][[j]]=lglasso(data = dd,x=covariates, rho = 0.5*rho[[1]][j],heter=T)
        }
      }else{
        results[[1]][[j]]=lglasso(data = dd, rho = 0.5*rho[[1]][j],heter=F)
      }
bb1=-2*results[[1]][[j]]$ll +  0.5*length(which(results[[1]][[j]]$omega!=0))*log(nrow(dd))*nu*1.2
bic1=c(bic1,bb1)
    }


##estimate the network based on glasso
s=cov(dd[,-c(1,2)])
for (j in 1:length(rho[[2]])) {
  results[[2]][[j]]=glasso(s=s,rho=rho[[2]][j])$wi
  bb2=bicfunction(data=dd,G=as.matrix(results[[2]][[j]]),nu=nu)
  bic2=c(bic2,bb2)

}



##estimate the network based on nh

for (j in 1:length(rho[[3]])) {
  GG=addition(data=dd,lambda=rho[[3]][j])
  results[[3]][[j]]=mle_net(data=dd,priori=GG)
  bb3=bicfunction(data=dd,G=results[[3]][[j]],nu=nu)
  bic3=c(bic3,bb3)

}

##estimate the network based on CO1
for (j in 1:length(rho[[4]])) {
  results[[4]][[j]]=selectFast(s,family="C01",K=rho[[4]][j])$C01$G
  comle=mle_net(data=dd,priori = results[[4]][[j]])
    bb4=bicfunction(data=dd,comle,nu=nu)
  bic4=c(bic4,bb4)
}



##estiamte the network based on LA
    for (j in 1:length(rho[[5]])) {
      results[[5]][[j]]=selectFast(s,family="LA",K=rho[[5]][j])$LA$G
      lamle=mle_net(data=dd,priori = results[[5]][[j]])
      bb5=bicfunction(data=dd,G=lamle,nu=nu)
      bic5=c(bic5,bb5)
    }


G1=results[[1]][[which.min(bic1)]]$omega
G2=results[[2]][[which.min(bic2)]]
G3=results[[3]][[which.min(bic3)]]
G4=results[[4]][[which.min(bic4)]]
G5=results[[5]][[which.min(bic5)]]

bb[[1]][h,]=as.numeric(comparison(graph,G1))
bb[[2]][h,]=as.numeric(comparison(graph,G2))
bb[[3]][h,]=as.numeric(comparison(graph,G3))
bb[[4]][h,]=as.numeric(comparison(graph,G4))
bb[[5]][h,]=as.numeric(comparison(graph,G5))

  }

aa=matrix(nrow=5,ncol=2)

  for (i in 1:5) {
      aa[i,]=apply(bb[[i]], 2, mean)
  }

bicc=cbind(bic1,bic2,bic3,bic4,bic5)
  return(aa)
}




