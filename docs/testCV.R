





library(devtools)
library(forcats)
library(fake)
library(CVglasso)
load_all()
source("./scripts/simulations.R")
set.seed(3)




# number of nodes
p=10
# number of edge in general network
m1=10
# the difference between the number of edges in individual networks and general network
m2=3
#  number of clusters
m3=4
# number of subjects
n=10
# correlation between cluster
rho=0.8





set.seed(5)
dd=Simulate(type="longihomo",n=n,p=p,m1=m1,tt=10,tau=c(2,1))
data1=cbind(1,dd$data$pre)
colnames(data1)[1]="group"
data2=cbind(2,dd$data$post)
colnames(data2)[1]="group"
ddata=rbind(data1,data2)


lambda1=exp(seq(-5,-1,length=3))
lambda2=lambda1[1:2]
lambda=expand.grid(lambda1,lambda2)

system.time(CVlglasso(data=ddata[,-1],lambda = lambda,group=ddata[,1],K=3,cores=3,trace = T))


system.time(CVlglasso(data=ddata[,-1],lambda = lambda,group=ddata[,1],K=3,cores=5,trace = T))

system.time(CVlglasso(data=ddata[,-1],lambda = lambda,group=ddata[,1],K=3,cores=14,trace = T))

