library(matrixcalc)
library(fake)
library(CVXR)
library(expm)
library(purrr)
library(stors)
library(MASS)
library(devtools)
#document()
#load_all()
library(lglasso)
p=20
m1=20
m2=0
m3=3
n=20
rho=0.8


## clustered data with no correlation structure

set.seed(5)
cc=matrix(rho,ncol = m3,nrow = m3)
diag(cc)=1
dd=Simulate(type="general",n=n,p=p,m1=m1,m2=m2,m3=m3,cc=cc)
ddata=dd$data
freq=table(ddata[,1])
group=c()
for (i in 1:length(freq)) {
  group=c(group,1:freq[i])
}
group=ifelse(group==1,1,2)
aa=lglasso(data=ddata,lambda = 0.1,type="general")
aa=lglasso(data=ddata,lambda = c(1,1),type="general",group = group)
lambda1=exp(seq(-2,0,length=3))
lambda2=lambda1
lambda=expand.grid(lambda1,lambda2)
system.time(cv.lglasso(type="general",data=ddata,lambda = lambda1,K=5, trace=TRUE))
system.time(cvp.lglasso(type="general",data=ddata,lambda = lambda1,K=5, trace=TRUE, cores=5))
system.time(cv.lglasso(type="general",group=group,data=ddata,lambda = lambda,K=5))
system.time(cvp.lglasso(type="general",group=group,data=ddata,lambda = lambda,K=5,cores=10))
plot.cvlglasso(bb)

## longitudinal data with structured homogeneous dampening rate
set.seed(5)
dd=Simulate(type="longihomo",n=n,p=p,m1=m1,tau=c(2,1))
ddata=dd$data
freq=table(ddata[,1])
group=c()
for (i in 1:length(freq)) {
  group=c(group,1:freq[i])
}
group=ifelse(group==1,1,2)
aa=lglasso(data=ddata,lambda = 1,type="expFixed",expFix=1,trace = TRUE)
aa=lglasso(data=ddata,lambda = c(1,1),type="expFixed",expFix=1,group = group, trace=TRUE)
lambda1=exp(seq(-2,0,length=3))
lambda2=lambda1
lambda=expand.grid(lambda1,lambda2)
system.time(cv.lglasso(type="expFixed",data=ddata,lambda = lambda1,K=5))
system.time(cvp.lglasso(type="expFixed",data=ddata,lambda = lambda1,K=5,cores=5))
system.time(cv.lglasso(type="expFixed",group=group,data=ddata,lambda = lambda,K=5))
system.time(cvp.lglasso(type="expFixed",group=group,data=ddata,lambda = lambda,K=5,cores=10))


## longitudinal data with structured homogeneous dampening rate

set.seed(5)
dd=Simulate(type="longiheter",n=n,p=p,m1=m1,alpha=2)
ddata=dd$data
freq=table(ddata[,1])
group=c()
for (i in 1:length(freq)) {
  group=c(group,1:freq[i])
}
group=ifelse(group==1,1,2)
aa=lglasso(data=ddata,lambda = 1,type="twoPara")
aa=lglasso(data=ddata,lambda = c(1,1),type="twoPara",group=group)
lambda1=exp(seq(-2,0,length=3))
lambda2=lambda1
lambda=expand.grid(lambda1,lambda2)
system.time(cv.lglasso(type="twoPara",data=ddata,lambda = lambda1,K=5))
system.time(cvp.lglasso(type="twoPara",data=ddata,lambda = lambda1,K=5,cores=5))
system.time(cv.lglasso(type="twoPara",group=group,data=ddata,lambda = lambda,K=5))
system.time(cvp.lglasso(type="twoPara",group=group,data=ddata,lambda = lambda,K=5,cores=10))




