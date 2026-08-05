





library(devtools)
library(forcats)
library(fake)
library(CVglasso)
load_all()
source("./scripts/simulations.R")
set.seed(3)

# number of nodes
p=50
# number of edge in general network
m1=122
# the difference between the number of edges in individual networks and general network
m2=12
#  number of clusters
m3=0
# number of subjects
n=40
# correlation between cluster
rho=0.8



## Model 1:  homogeneous individuals

set.seed(5)
dd=Simulate(type="longihomo",n=n,p=p,m1=m1,tt=10,tau=c(2,1))
data1=cbind(1,dd$data$pre)
colnames(data1)[1]="group"
data2=cbind(2,dd$data$post)
colnames(data2)[1]="group"
ddata=rbind(data1,data2)


lambda1=exp(seq(-5,-1,length=5))
lambda2=lambda1[1:2]
lambda=expand.grid(lambda1,lambda2)


### two-stage model

bb=CVlglasso(data=ddata[,-1],nlam = 10,group=ddata[,1],K=2,trace = T)
plot(bb)

## Model 2:  heterogeneous individuals


set.seed(5)
dd=Simulate(type="longiheter",n=n,p=p,m1=m1,alpha=2,group=2)
data1=cbind(1,dd$data$pre)
colnames(data1)[1]="group"
data2=cbind(2,dd$data$post)
colnames(data2)[1]="group"
ddata=rbind(data1,data2)

lambda1=exp(seq(-5,0,length=5))
lambda2=exp(seq(-5,-0,length=2))
lambda=expand.grid(lambda1,lambda2)


### one-stage model

bb=CVlglasso(data=data1[,-1],nlam = 50,K=5,random=TRUE)
plot(bb)




