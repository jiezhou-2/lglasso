

#' Main simulation function comparing five methods for network selection.
#'
#' @param m number of subjects
#' @param n  number of observations per subjects
#' @param p the dimension of the data to be generated
#' @param coe coefficients for covariates
#' @param l the simulation scale
#' @param rho1 tuning parameter for glasso
#' @param rho2 tuning parameter for glasso
#' @return list with length equal to 3.
#' @export

power_compare1=function(m,n,p,coe,l,rho,prob,heter,nu){
  results=vector("list",length = 5) # container for the final FPR and TPR
  results[[1]]=vector("list",length = length(rho[[1]]))
  results[[2]]=vector("list",length = length(rho[[2]]))
  results[[3]]=vector("list",length = length(rho[[3]]))
  results[[4]]=vector("list",length = length(rho[[4]]))
  results[[5]]=vector("list",length = length(rho[[5]]))
  RR=vector("list",length = 5)
  RR[[1]]=matrix(nrow=l,ncol = 2)
  RR[[2]]=matrix(nrow=l,ncol = 2)
  RR[[3]]=matrix(nrow=l,ncol = 2)
  RR[[4]]=matrix(nrow=l,ncol = 2)
  RR[[5]]=matrix(nrow=l,ncol = 2)

  age=vector("list",m) #generate the time points container
  for (i in 1:5) {
    for (j in 1:length(rho[[i]])) {
      results[[i]][[j]]= matrix(nrow = l,ncol = 2)
    }
  }


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
      ss=sim_homo(p = p,prob=prob,tau = 1/alpha,age = age)
    }
    simdata=ss$data
    tau=ss$tau
    lower=min(abs(tau))
    upper=max(abs(tau))
    graph=ss$precision
    dd=do.call(rbind,simdata)
    id=unique(dd[,1])
    covariate=cbind(id,x)


    bic1=c()
    bic2=c()
    bic3=c()
    bic4=c()
    bic5=c()
    aa=vector("list",length(rho[[1]]))
    aa1=vector("list",length(rho[[2]]))
    aa2=vector("list",length(rho[[3]]))
    aa3=vector("list",length(rho[[4]]))
    aa4=vector("list",length(rho[[5]]))

    ## lglasso
    for (j in 1:length(rho[[1]])) {
      if (heter==TRUE){
        if (length(which(coe==0))==2){
          aa[[j]]=lglasso(data = dd, rho = 0.5*rho[[1]][j],heter=T)
        }else{
          aa[[j]]=lglasso(data = dd,x=covariate, rho = 0.5*rho[[1]][j],heter=T)
        }
      }else{
        aa[[j]]=lglasso(data = dd, rho = 0.5*rho[[1]][j])
      }
      bb1=-2*aa[[j]]$ll +  0.5*length(which(aa[[j]]$omega!=0))*log(nrow(dd))*nu*3
      bic1=c(bic1,bb1)
      results[[1]][[j]][h,]=as.numeric(comparison(graph,aa[[j]]$omega))
    }

    print(paste0("lglasso's bic is ",bic1))

    ## glasso

    for (j in 1:length(rho[[2]])) {
      ##estiamte the network based on glasso
      s=cov(dd[,-c(1,2)])
      aa1[[j]]=glasso(s=s,rho=rho[[2]][j])$wi
      bb2=bicfunction(data=dd,G=as.matrix(aa1[[j]]),nu=3*nu)
      bic2=c(bic2,bb2)
      results[[2]][[j]][h,]=as.numeric(comparison(graph,aa1[[j]]))
    }


    ## nh

    for (j in 1:length(rho[[3]])) {
      bb2=addition(data=dd,lambda=rho[[3]][j])
      aa2[[j]]=mle_net(data=dd,priori=bb2)
      bb3=bicfunction(data=dd,G=aa2[[j]],nu=nu)
      bic3=c(bic3,bb3)
      results[[3]][[j]][h,]=as.numeric(comparison(graph,aa2[[j]]))
    }
    ## CO1
    for (j in 1:length(rho[[4]])) {
      aa3[[j]]=selectFast(s,family="C01",K=2*rho[[4]][j])$C01$G
      bb4=bicfunction(data=dd,G=aa3[[j]],nu=nu)
      bic4=c(bic4,bb4)
      results[[4]][[j]][h,]=as.numeric(comparison(graph,aa3[[j]]))

    }

    ## la
    for (j in 1:length(rho[[5]])) {
      aa4[[j]]=selectFast(s,family="LA",K=2*rho[[5]][j])$LA$G
      bb5=bicfunction(data=dd,G=aa4[[j]],nu=nu)
      bic5=c(bic5,bb5)
      results[[5]][[j]][h,]=as.numeric(comparison(graph,aa4[[j]]))
    }

    G1=aa[[which.min(bic1)]]$omega
    #print(bic2)
    G2=aa1[[which.min(bic2)]]
    #print(bic3)
    G3=aa2[[which.min(bic3)]]
    G4=aa3[[5]]
    G5=aa4[[5]]


    RR[[1]][h,]=as.numeric(comparison(graph,G1))
    RR[[2]][h,]=as.numeric(comparison(graph,G2))
    RR[[3]][h,]=as.numeric(comparison(graph,G3))
    RR[[4]][h,]=as.numeric(comparison(graph,G4))
    RR[[5]][h,]=as.numeric(comparison(graph,G5))
  }


  bb=vector("list",length = 5)
  bb[[1]]=matrix(nrow = length(rho[[1]]),ncol=2)
  bb[[2]]=matrix(nrow = length(rho[[2]]),ncol=2)
  bb[[3]]=matrix(nrow = length(rho[[3]]),ncol=2)
  bb[[4]]=matrix(nrow = length(rho[[4]]),ncol=2)
  bb[[5]]=matrix(nrow = length(rho[[5]]),ncol=2)
  names(bb)=c("lglasso","glasso","nh","co","la")
  for (i in 1:5) {
    for (j in 1:length(rho[[i]])) {
      iimatrix=results[[i]][[j]]
      index1=which( is.na(iimatrix[,1]))
      index2=which(is.na(iimatrix[,2]))
      index=union(index1,index2)
      if (length(index)==l){
        warning(paste0(" All results are NA for rho= ", i,"+",j,rho[[i]][j]))
        next
      }
      ff=if(length(index)==0){
        ff=iimatrix
      }else{
        ff=iimatrix[-index,]
      }
      bb[[i]][j,]=apply(ff, 2, mean)
    }
  }

  rr=matrix(nrow=5,ncol=2)

  for (i in 1:5) {
    rr[i,]=apply(RR[[i]], 2, mean)
  }


  return(list(whol=bb,selec=rr,resul=results))
}




