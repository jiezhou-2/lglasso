

library(glasso)
library(lglasso)
library(BDgraph)
library(GGMselect)
rho=vector("list",5)
rho[[1]]=seq(0.001,0.3,length=20)
rho[[2]]=seq(0.001,0.3,length=20)
rho[[3]]=seq(0.001,0.3,length=20)
rho[[4]]=seq(0.00001,0.1,length=10)
rho[[5]]=seq(0.00001,0.1,length=10)

set.seed(51)
results=power_compare1(m=20,n=5,p=80,coe = c(2,0,0),l=10,rho = rho,prob=0.01,heter=T,nu=0.5)
save(results)
plot(results[[1]][[2]][,2],results[[1]][[2]][,1],ylim = c(0.3,1),xlim=c(0,1),type="l",xlab = "FPR", ylab = "TPR")
lines(results[[1]][[3]][,2],results[[1]][[3]][,1],col=2)
lines(results[[1]][[1]][,2],results[[1]][[1]][,1],col=3)




m=20
n=5
p=20
coe = c(2,0,0)
l=10

prob=0.1
heter=T
nu=10












