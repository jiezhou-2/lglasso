library(lglasso)
# number of nodes
p=15
# number of edge in general network
m1=20
# the difference between the number of edges in individual networks and general network
m2=5
# number of subjects
n=30
set.seed(1)
## One-stage model
### Estimate the network based on homogeneous one-stage model
####simulate data
dd=Simulate(type="homo",n=n,p=p,m1=m1,m2=m2,tau=2,tt=10)
ddata=dd$data
dim(ddata)
ddata[1:2,1:5]
#### Estimation
 aa=lglasso(data=ddata,lambda = 0.01,trace=TRUE)
# estimates=lapply(aa$wi,function(ll){ifelse(abs(ll)>10^(-5),1,0)})
# #### estimated network
# estimates
# #### true network
# dd$network
# #### correlation parameter
# aa$tau
# #### likelihood
# aa$ll


