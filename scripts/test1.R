
m=20
n=20
p=80
coe = c(2,0,0)
l=3
rho1 = 0.1
rho2= 0.1
prob=0.01
heter=T


x1=sample(x=c(0,1),size=m,prob=c(0.5,0.5),replace = T)
x2=runif(m,min = 0,max = 1)
#alpha=exp(coe[1]+coe[2]*x1+coe[3]*x2)
alpha=rep(0.916,m)
x=cbind(x1,x2)
age=vector("list",m)
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
lower=0.01
upper=5
graph=ss$precision
dd=do.call(rbind,simdata)
id=unique(dd[,1])
covariate=cbind(id,x)

##homo model
x=seq(lower,upper,length=50)
z=sapply(x, ll_homo,data=data,omega=graph,type=1)
plot(x,z,type = "l")

##heter
x=seq(lower,upper,length=50)
tauhat=c()
for (i in 1:length(id)) {
idata=dd[which(dd[,1]==id[i]),]
z=sapply(x, logdensity,idata=idata,omega=graph,alpha=alpha[1],type=1)
tauhat=c(tauhat,x[min(which(z==max(z)))])
}
plot(x,z,type = "l")
## covariate
idata=dd[which(dd[,1]==id[1]),]
icovariates=covariate[1,]
coe0=coe
x=seq(lower,upper,length=50)
z=sapply(x, covalogdensity,idata=idata,icovariates=icovariates,omega=graph,coe=coe0,type=type)
tau1[i]=x[min(which(z==max(z)))]
