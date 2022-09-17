

library(glasso)
library(lglasso)
library(BDgraph)
library(GGMselect)
rho=vector("list",5)
rho[[1]]=seq(0.001,3.5,length=10)
rho[[2]]=seq(0.001,0.3,length=10)
rho[[3]]=seq(0.001,0.3,length=10)
rho[[4]]=seq(0.00001,0.1,length=10)
rho[[5]]=seq(0.00001,0.1,length=10)

#set.seed(4)
results=power_compare1(m=20,n=5,p=20,coe = c(2,0,0),l=5,rho = rho,prob=0.1,heter=T)
plot(results[[2]][,2],results[[2]][,1],ylim = c(0.3,1),xlim=c(0,1),type="l")
lines(results[[3]][,2],results[[3]][,1],col=2)
lines(results[[1]][,2],results[[1]][,1],col=3)
lines(results[[4]][,2],results[[4]][,1],col=4)
lines(results[[5]][,2],results[[5]][,1],col=5)



m=20
n=5
p=20
coe = c(2,0.5,0.5)
l=5

prob=0.1
heter=T













