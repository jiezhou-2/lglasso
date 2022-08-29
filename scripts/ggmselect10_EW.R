
source("joint1.R")
source("joint2.R")
library(glasso)
library(GGMselect)
source("simulation_bdgraph_2.R")
## generate the data
m=20
ni=20
p=80
alpha=14
age=vector("list",m)
l=10
rho=0.08
pool1=matrix(nrow = l,ncol = 2)
for (h in 1:l){
  print(h)
for (k in 1:m) {
  a1=rpois(n=ni,lambda = 1)
  age[[k]][1]=max(a1[1],0.5)
  for (i in 2:ni) {
    age[[k]][i]=age[[k]][i-1]+max(a1[i-1],0.5)
  }
}

ss=sim_heter(p = p,prob=0.01,alpha = alpha,age = age,type="abs")
data=ss$data
tau=ss$tau
graph=ss$precision
dd=data.frame()
for (i in 1:m) {
  dd=rbind(dd,data[[i]])
}
s=cov(dd[,-c(1,2)])
aa1=selectFast(s,family="EW",K=0.6)$EW$G
  pool1[h,]=as.numeric(comparison(graph,aa1)) 
}
a1=summary(pool1[,1])
a2=summary(pool1[,2])
a=rbind(a1,a2)
write.csv(a,file="a10_EW.csv")


