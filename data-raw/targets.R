library(targets)
library(tarchetypes)
library(crew)
library(MASS)
library(glasso)
library(tibble)
source("./data-raw/sim.R")
source("./data-raw/roc_version2.R")


### globals
num_reps <- 10L
n <- 20L
tt=4
p <- 60L
m1 <- 200
m2 <- 20
rho_cov <- 0.5
tau <- c(1,1)
rho_seq <- rbind(seq(0.1, 0.5, length.out = 10),
                 seq(0.1, 0.5, length.out = 10))
rho_seq_glasso <- matrix(seq(0.1, 0.5, length.out = 10),nrow=1)
sigmaM=sim_stru(p=p,m1=m1,m2=m2)
simData=sim_data(n=n,tt=tt,tau=tau,covmat=sigmaM[3:4])
preData=simData$data$pre
postData=simData$data$post
fullData=rbind(rbind(preData,postData))
groupIndex=c(rep(1,nrow(preData)),rep(2,nrow(postData)))


