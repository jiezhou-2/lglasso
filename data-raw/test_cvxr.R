library(matrixcalc)
library(fake)
library(CVXR)
library(expm)
library(purrr)
library(stors)
library(MASS)
library(devtools)
load_all()
p=50
m1=100
m2=0
m3=5
n=20
rho=0.8
set.seed(5)

## clustered data with no correlation structure
cc=matrix(rho,ncol = m3,nrow = m3)
diag(cc)=1
dd=simulate(n=n,p=p,m1=m1,m2=m2,m3=m3,cc=cc,homo=T)
ddata=dd$data
aa=lglasso(data=ddata,lambda = c(0.1,0.1),type="general",group = ddata[,2])
aa=lglasso(data=ddata,lambda = 0.1,type="general")
 cov2cor(aa$corMatrix)
#  homo=aa$preMatrix
# results=heternetwork(data=ddata,lambda = c(0.1,0.1),homo=homo)

## longitudinal data with structured homogeneous dampening rate
 set.seed(5)
dd=simulate_long(n=n,p=p,m1=m1,tau=c(2,1))
ddata=dd$data
group=rep(c(1,2),71)
aa=lglasso(data=ddata,lambda = c(0.1,0.1),type="expFixed",expFix=1,group = group)
aa=lglasso(data=ddata,lambda = 0.1,,type="expFixed",expFix=1)
aa$tauhat
cov2cor(aa$corMatrixList[[1]])


## longitudinal data with structured homogeneous dampening rate
set.seed(5)
dd=simulate_long(n=n,p=p,m1=m1,tau=c(2,1))
ddata=dd$data
group=rep(c(1,2),71)
aa=lglasso(data=ddata,lambda = c(0.001,0.1),type="twoPara",group=group)
aa=lglasso(data=ddata,lambda = 0.001,type="twoPara")
aa$tauhat
cov2cor(aa$corMatrixList[[1]])





