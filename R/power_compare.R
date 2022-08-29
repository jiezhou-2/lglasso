

#' Main simulation function comparing five methods for network seletion.
#'
#' @param m number of subjects
#' @param n  number of observations per subjects
#' @param p the dimension of the data to be generated
#' @param coe coefficients for covariates
#' @param l the simulation scale
#' @param rho tuning parameter for glasso
#'
#' @return list with length equal to 5.
#' @export

power_compare=function(m,n,p,coe,l,rho,prob){
results=vector("list",length = 5) # container for the final FPR and TPR
age=vector("list",m) #generate the time points container
for (i in 1:5) {
  results[[i]]= matrix(nrow = l,ncol = 2)
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
  ss=sim_heter(p = p,prob=prob,alpha = alpha,age = age)
  simdata=ss$data
  tau=ss$tau
  lower=min(abs(tau))
  upper=max(abs(tau))
  graph=ss$precision
  dd=do.call(rbind,simdata)
  id=unique(dd[,1])
  covariate=cbind(id,x)
  ##estiamte the network based on lglasso
  aa=lglasso(data = dd,x=covariate, rho = rho,heter=T, type=1)
  results[[1]][h,]=as.numeric(comparison(graph,aa$omega))
  ##estiamte the network based on glasso
  s=cov(dd[,-c(1,2)])
  aa1=glasso(s=s,rho=rho)$wi
  results[[2]][h,]=as.numeric(comparison(graph,aa1))
  ##estiamte the network based on nh
  aa2=addition(data=dd,lambda=rho/2)
  results[[3]][h,]=as.numeric(comparison(graph,aa2))
  ##estiamte the network based on CO1
  aa3=selectFast(s,family="C01",K=0.00001)$C01$G
  results[[4]][h,]=as.numeric(comparison(graph,aa3))
  ##estiamte the network based on EW
  aa4=selectFast(s,family="LA",K=0.00001)$LA$G
  results[[5]][h,]=as.numeric(comparison(graph,aa4))
}
bb=matrix(nrow = 5,ncol = 2)
for (i in 1:5) {
  iimatrix=results[[i]]
  index1=which(iimatrix[,1]==1 | iimatrix[,1]==0)
  index2=which(iimatrix[,2]==1 | iimatrix[,2]==0)
  index=union(index1,index2)
  ff=if(length(index)==0){
    ff=iimatrix
    }else{
      ff=iimatrix[-index,]
}
  bb[i,]=apply(ff, 2, mean)
  }
return(bb)
}




