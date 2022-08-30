

library(glasso)
library(lglasso)
library(BDgraph)
library(GGMselect)
#set.seed(1)
results=power_compare(m=20,n=10,p=80,coe = c(2,0.5,0.5),l=10,rho = 0.1,prob=0.01,heter=T)






