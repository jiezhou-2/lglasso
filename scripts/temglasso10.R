

library(glasso)
library(lglasso)
library(BDgraph)
## generate the data
m=20 ## number of subjects
ni=20 # number of observations per subjects
p=40 # the dimension of the data to be generated
alpha=10 #the parameter for exponential distribution
age=vector("list",m) #generate the time points container
l=30 # the simulation scale
rho=0.07 # tuning parameter for glasso
pool1=matrix(nrow = l,ncol = 2)
pool2=matrix(nrow = l,ncol = 2)
cograph=c(0.1,0.2,0.2)
for (h in 1:l){
print(paste0("The ",h,"th"," simulation:"))
  # subject level covariates
  x1=sample(x=c(0,1),size=m,prob=c(0.5,0.5),replace = T)
  x2=runif(m,min = 0,max = 1)
  alpha=exp(cograph[1]+cograph[2]*x1+cograph[3]*x2)
  x=cbind(x1,x2)
  ## generate the observation time for each subjects
for (k in 1:m) {
  a1=rpois(n=ni,lambda = 1) # the space between observations
  age[[k]][1]=max(a1[1],0.5)
  for (i in 2:ni) {
    age[[k]][i]=age[[k]][i-1]+max(a1[i-1],0.5)
  }
}

## generate the network data
prob=0.01 ## the edge density of the network
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
#alpha0=1
#omega0=diag(p)
a=lglasso(data = dd,x=covariate, rho = rho,heter=T, type=1)
mlenet_he=mle(data =as.matrix(dd), network = a$omega,type = 1,heter = F)
alpha_tau_mle=mlenet_he$alpha
tau_he=mlenet_he$tau
pool1[h,]=as.numeric(comparison(graph,a$omega))
rerror=mean(abs(tau-tau_he)/tau)
cotau=cor(tau,tau_he)
pool2[h,1]=rerror
pool2[h,2]=cotau
}

 a1=summary(pool1[,1])
 a2=summary(pool1[,2])
 b1=summary(pool2[,1])
 b2=summary(pool2[,2])
 a=rbind(a1,a2,b1,b2)
 write.csv(a,file="tem10.csv")




