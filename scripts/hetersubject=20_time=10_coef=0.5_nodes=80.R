library("lglasso")
library("glasso")
library("BDgraph")
library("GGMselect")
library("Matrix")
m=20
n=10
e=0
p=80
Nsim=50
rho=vector("list",5)
rho[[1]]=seq(0.001,0.3,length=20)
rho[[2]]=seq(0.001,0.3,length=20)
rho[[3]]=seq(0.001,0.3,length=20)
rho[[4]]=seq(0.00001,0.1,length=10)
rho[[5]]=seq(0.00001,0.1,length=10)
set.seed(m+n+e+p)
#debug(power_compare1)
simures=power_compare1(m=m,n=n,p=p,coe=c(2,e,e),l=Nsim,rho=rho,prob=0.01,heter=T,community2 = T,uu=c(5,10),zirate = c(0.2,0.6))
save(simures,file="subject=20_time=10_coef=0.5_nodes=80.Rd")
q("no")



coe=c(2,e,e)
l=Nsim
rho=rho
prob=0.01
heter=T
community2 = T
uu=c(5,10)
zirate=c(0.2,0.6)


idata=idata
omega=omega0
alpha=alpha0[i]
type=type
