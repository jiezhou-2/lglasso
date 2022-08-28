
source("joint1.R")
source("joint2.R")
library(glasso)
source("simulation_bdgraph_2.R")
## generate the data
m=1
ni=100
p=80
alpha=10
age=vector("list",m)
l=500
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
##homogeneous for temglasso+ebic
  ##aa1=longraph_homo(data = dd,tau0 = tau,omega0 = pp,rho = rho,tole = 0.001, lower = 5,upper = 5, type = "abs")$omega
aa1=glasso(s=s,rho=rho)$wi
  pool1[h,]=as.numeric(comparison(graph,aa1)) 
}
a1=summary(pool1[,1])
a2=summary(pool1[,2])
a=rbind(a1,a2)
write.csv(a,file="a10.csv")


