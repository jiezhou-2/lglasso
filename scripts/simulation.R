

library(glasso)
library(lglasso)
library(BDgraph)
library(GGMselect)
#set.seed(1)
results=power_compare(m=20,n=5,p=80,coe = c(2,1,1),l=10,rho = 0.3,prob=0.05)



m=10
n=20
p=80
coe = c(2,1,1)
l=10
rho = 0.05
prob=0.1

