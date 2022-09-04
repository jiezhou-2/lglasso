


#' Title
#'
#' @param m number of subjects
#' @param n number of observations per subjects
#' @param p the dimension of the data to be generated
#' @param coe coefficients for covariates
#' @param l the simulation scale
#' @param rho tuning parameter for glasso
#' @param prob the density of the network
#' @param heter indicator
#'
#' @return list
#' @export
power_compare2=function(m,n,p,coe,l,rho,prob,heter){
  results=vector("list",length = 2) # container for the final FPR and TPR
  results[[1]]=vector("list",length = length(rho))
  results[[2]]=vector("list",length = length(rho))
  age=vector("list",m) #generate the time points container
  for (i in 1:2) {
    for (j in 1:length(rho)) {
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

    s=cov(dd[,-c(1,2)])


    for (j in 1:length(rho)) {
      print(rho[j])
      ##estimate the network based on CO1
      aa3=selectFast(s,family="C01",K=2*rho[j])$C01$G
      results[[1]][[j]][h,]=as.numeric(comparison(graph,aa3))
      ##estiamte the network based on EW
      aa4=selectFast(s,family="LA",K=2*rho[j])$LA$G
      results[[2]][[j]][h,]=as.numeric(comparison(graph,aa4))
    }
  }
  bb=vector("list",length = 2)
  names(bb)=c("co","la")
  bb[[1]]=matrix(nrow = length(rho),ncol=2)
  bb[[2]]=matrix(nrow = length(rho),ncol=2)

  for (i in 1:2) {
    for (j in 1:length(rho)) {
      iimatrix=results[[i]][[j]]
      index1=which( is.na(iimatrix[,1]))
      index2=which(is.na(iimatrix[,2]))
      index=union(index1,index2)
      if (length(index)==l){stop(paste0(" All results are NA for rho= ", rho[j]))}
      ff=if(length(index)==0){
        ff=iimatrix
      }else{
        ff=iimatrix[-index,]
      }

      bb[[i]][j,]=apply(ff, 2, mean)
    }

  }
  return(bb)
}


