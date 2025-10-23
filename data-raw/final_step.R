source("./data-raw/roc_version2.R")
source("./data-raw/sim.R")
library(glasso)
library(EstimateGroupNetwork)
library(lglasso)
library(CVXR)
library(data.table)
rho=vector("list",3)
rho[[1]]=rbind(seq(0.1,0.9,length=10),seq(0.1,0.9,length=10))
rho[[2]]=rbind(seq(0.01,0.09,length=10),seq(0.01,0.09,length=10))
rho[[3]]=rbind(seq(0.01,0.09,length=10),seq(0.01,0.09,length=10))
results=power_compare1(n=20,tt=2,p=12,m1=10,m2=3,tau=c(1,1),l=10,rho=rho)

library(targets)
library(crew)
library(tarchetypes)                                                                                                                
tar_make()
aa=tar_read(simulated_data)
aa=tar_read(all_results)
aa=tar_read(glasso_fit)
### lglasso
k=0
aa=0
bb=0
for (i in 1:length(results_lglasso)) {
  for (j in 1:length(results_lglasso[i])) {
    aa=ifelse(abs(results_lglasso[[i]][[j]]$wi[[1]])<=10^(-5),0,1)+aa
    bb=ifelse(abs(results_lglasso[[i]][[j]]$wi[[2]])<=10^(-5),0,1)+bb
    k=k+1
  }
}
pre_prob_lglasso=aa/k
post_prob_lglasso=bb/k

### glasso
k=0
aa=0
bb=0
for (i in 1:length(results_glasso)) {
  for (j in 1:length(results_glasso[i])) {
    aa=ifelse(abs(results_glasso[[i]]$pre[[j]]$wi)<=10^(-5),0,1)+aa
    bb=ifelse(abs(results_glasso[[i]]$post[[j]]$wi)<=10^(-5),0,1)+bb
    k=k+1
  }
}
pre_prob_glasso=aa/k
post_prob_glasso=bb/k


### jgl
k=0
aa=0
bb=0
for (i in 1:length(results_jgl)) {
  for (j in 1:length(results_jgl[i])) {
    aa=ifelse(abs(results_jgl[[i]][[j]][[1]])<=10^(-5),0,1)+aa
    bb=ifelse(abs(results_jgl[[i]][[j]][[2]])<=10^(-5),0,1)+bb
    k=k+1
  }
}
pre_prob_jgl=aa/k
post_prob_jgl=bb/k




