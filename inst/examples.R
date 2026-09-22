library(lglasso)
library(fake)
# number of nodes
p=20
# number of edge in general network
m1=100
# the difference between the number of edges in individual networks and general network
m2=20
# number of subjects
n=100
set.seed(1)
## One-stage model
### Estimate the network based on homogeneous one-stage model
####simulate data
dd=lglasso:::Simulate(type="longihomo",n=n,p=p,m1=m1,m2=m2,tau=2,tt=10)
ddata=dd$data
dim(ddata)
ddata[1:2,1:5]
#### Estimation
aa=lglasso(data=ddata,lambda = 0.01,trace=TRUE)
#### estimated network
estimates=lapply(aa$wi,function(ll){ifelse(abs(ll)>10^(-5),1,0)})
#### estimated network
estimates
#### true network
dd$network
#### correlation parameter
aa$tau
#### likelihood
aa$ll



### Estimate the network based on heterogeneous one-stage model
####simulate data
dd=lglasso:::Simulate(type="longiheter",n=n,p=p,m1=m1,m2=m2,tt=10,alpha=5,group = 1)
ddata=dd$data$pre
dim(ddata)
ddata[1:2,1:5]
#### Estimation
aa=lglasso(data=ddata,lambda = 0.01,random=TRUE,trace=TRUE,N=100)
estimates=lapply(aa$wi,function(ll){ifelse(abs(ll)>10^(-5),1,0)})
#### estimated network
estimates
####  true network
dd$network
#### likelihood
aa$ll



## Two-stage model
### Estimate the networks based on homogeneous two-stage model
####simulte data
dd=lglasso:::Simulate(type="longihomo",n=n,p=p,m1=m1,m2=m2,tau=c(2,1),tt=10)
ddata=do.call(rbind,dd$data)
group=c(rep(0,nrow(ddata)/2),rep(1,nrow(ddata)/2))
dim(ddata)
ddata[1:2,1:5]
#### Estimation
aa=lglasso(data=ddata,lambda = c(0.01,0.01),group = group,trace=TRUE)
estimates=lapply(aa$wi,function(ll){ifelse(abs(ll)>10^(-5),1,0)})
#### estimated pre-treatment  network
estimates[[1]]
####  true pre-treatment network
dd$network$pre
#### estimated post-treatment  network
estimates[[2]]
####  true post-treatment network
dd$network$post
#### correlation parameters
aa$tau
#### likelihood
aa$ll


### Estimate the networks based on hetergeneous two-stage model
####simulate data
dd=lglasso:::Simulate(type="longiheter",n=n,p=p,m1=m1,m2=m2,tt=10,alpha=0.5,group=2)
ddata=do.call(rbind,dd$data)
group=c(rep(0,nrow(ddata)/2),rep(1,nrow(ddata)/2))
dim(ddata)
ddata[1:2,1:5]
#### Estimation
aa=lglasso(data=ddata,lambda = c(0.01,0.01),random=TRUE,group = group,trace=TRUE,N=100)
estimates=lapply(aa$wi,function(ll){ifelse(abs(ll)>10^(-5),1,0)})
#### estimated pre-treatment network
estimates[[1]]
#### true pre-treatment network
dd$network$pre
#### estimated post-treatment network
estimates[[2]]
####  true post-treatment network
dd$network$post
#### likelihood
aa$ll




