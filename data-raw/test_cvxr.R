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


## longitudinal data with structured homogeneous dampening rate
 set.seed(5)
dd=simulate_long(n=n,p=p,m1=m1,tau=c(2,1))
ddata=dd$data
group=rep(c(1,2),71)
aa=lglasso(data=ddata,lambda = c(0.1,0.1),type="expFixed",expFix=1,group = group, start="cold")
aa=lglasso(data=ddata,lambda = 0.1,,type="expFixed",expFix=1)
aa$tauhat
cov2cor(aa$vList[[1]])


## longitudinal data with structured homogeneous dampening rate
set.seed(5)
dd=simulate_long(n=n,p=p,m1=m1,tau=c(2,1))
ddata=dd$data
group=rep(c(1,2),71)
aa=lglasso(data=ddata,lambda = c(0.001,0.1),type="twoPara",group=group, trace = TRUE)
aa=lglasso(data=ddata,lambda = 0.001,type="twoPara")
aa$tauhat
cov2cor(aa$vList[[1]])





