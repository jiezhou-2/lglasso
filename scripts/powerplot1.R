library("lglasso")
 library("glasso")
 library("BDgraph")
 library("GGMselect")
 library("Matrix")
 uu1=0.1
 uu2=0.8
 Nsim=50
 rho=vector("list",5)
 rho[[1]]=seq(0.01,0.2,length=20)
 rho[[2]]=seq(0.01,0.2,length=20)
 rho[[3]]=seq(0.01,0.2,length=20)
 rho[[4]]=seq(0.001,0.1,length=10)
 rho[[5]]=seq(0.001,0.1,length=10)
 set.seed(uu2*30+uu1*10)
 debug(power_compare1)
 simures=power_compare1(m=20,n=20,p=80,coe=c(2.3,0,0),l=Nsim,rho=rho,prob=0.01,heter=T,community2=F,zirate=c(uu1, uu2))
 save(simures,file="uu1=0.1_uu2=0.8.Rd")
