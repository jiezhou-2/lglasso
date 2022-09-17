

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

power_compare2=function(m,n,p,coe,l,rho1,rho2,prob,heter){
  results=vector("list",length = 5) # container for the final FPR and TPR
  results[[1]]=vector("list",length = length(rho1))
  results[[2]]=vector("list",length = length(rho1))
  results[[3]]=vector("list",length = length(rho1))
  results[[4]]=vector("list",length = length(rho2))
  results[[5]]=vector("list",length = length(rho2))
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
    for (j in 1:length(rho1)) {
print(paste0("rho1's ",j, " componts: ", rho1[j] ))
      ##estiamte the network based on lglasso
      if (heter==TRUE){
        if (length(which(coe==0))==2){
          results[[1]][[j]]=lglasso(data = dd, rho = 0.5*rho1[j],heter=T)
        }else{
          results[[1]][[j]]=lglasso(data = dd,x=covariates, rho = 0.5*rho1[j],heter=T)
        }
      }else{
        results[[1]][[j]]=lglasso(data = dd, rho = 0.5*rho1[j],heter=F)
      }
bb1=-2*results[[1]][[j]]$ll+  0.5*length(which(results[[1]][[j]]$omega!=0))*log(nrow(dd))
bic1=c(bic1,bb1)
      ##estiamte the network based on glasso
      s=cov(dd[,-c(1,2)])
      results[[2]][[j]]=glasso(s=s,rho=rho1[j])$wi
      #G2=glasso(s=s,rho=rho1[j])$wi
      #results[[2]][[j]]=mle_net(dd,priori=G2)
      bb2=bicfunction(data=dd,G=as.matrix(results[[2]][[j]]))
      bic2=c(bic2,bb2)
      ##estiamte the network based on nh
      results[[3]][[j]]=addition(data=dd,lambda=rho1[j])
     #G3=addition(data=dd,lambda=rho1[j])
     #results[[3]][[j]]=mle_net(data=dd,priori=G3)
      bb3=bicfunction(data=dd,G=results[[3]][[j]])
      bic3=c(bic3,bb3)
    }

    for (j in 1:length(rho2)) {
      ##estimate the network based on CO1
      results[[4]][[j]]=selectFast(s,family="C01",K=rho2[j])$C01$G
      #G4=selectFast(s,family="C01",K=rho2[j])$C01$G
      #results[[4]][[j]]=mle_net(data=dd,priori=G4)
      bb4=bicfunction(data=dd,G=results[[4]][[j]])
      bic4=c(bic4,bb4)
      ##estiamte the network based on EW
      results[[5]][[j]]=selectFast(s,family="LA",K=rho2[j])$LA$G
      #G5=selectFast(s,family="LA",K=rho2[j])$LA$G
      #results[[5]][[j]]=mle_net(data=dd,priori=G5)
      bb5=bicfunction(data=dd,G=results[[5]][[j]])
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
  return(bicc)
}




