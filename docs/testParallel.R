


library(CVXR)
library(devtools)
source("/Users/f003r0s/Desktop/Anne_lab_projects/lglasso/docs/sim.R")
load_all("/Users/f003r0s/Desktop/Anne_lab_projects/lglasso")


## globals
n=20
p=10
m1=12
m2=3
tt1=5
tt2=5
tau1=1
tau2=1
nn=2





# Define the pipeline

    sigmaM=sim_stru(p=p,m1=m1,m2=m2)


    timepoints=sim_timepoints(n=n,tt=c(tt1,tt2))

    phiM = sim_phi(timepoints,tau=c(tau1,tau2))


      results = ebic_selection(phiM=phiM,sigmaM=sigmaM,timepoints=timepoints,
                               rho_lglasso = c(0.5,0.5),
                               rho_glasso =  c(0.5,0.5),
                               rho_jgl =  c(0.5,0.5),,tau=c(1,1))




