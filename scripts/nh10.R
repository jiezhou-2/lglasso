source("joint1.R")
source("joint2.R")
source("candidateset.R")
source("comparison.R")
source("simulation_bdgraph_2.R")
library(glasso)
## generate the data
m=30
ni=7
p=80
alpha=10
age=vector("list",m)
l=500
rho=0.1
result=c()
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
  ##homogeneous for temglasso+ebic
  sglasso=cov(dd[,-c(1,2)])
    aa=addition(data=dd[,-c(1,2)],lambda=rho)
    result=rbind(result,as.numeric(comparison(graph,aa))) 
}
a1=summary(result[,1])
a2=summary(result[,2])
a=rbind(a1,a2)
write.csv(a,file="n10.csv")

