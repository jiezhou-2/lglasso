

library(glasso)
library(lglasso)
library(BDgraph)
library(GGMselect)
#set.seed(100)
rho1=seq(0.01,0.3,length=10)
rho2=seq(0.1,1,length=5)
results=power_compare1(m=20,n=20,p=80,coe = c(0.916,0,0),l=5,rho1 = rho1,rho2= 0.1,prob=0.01,heter=T)
plot(results[[2]][,2],results[[1]][,1],ylim = c(0.25,1),type="l")
lines(results[[3]][,2],results[[1]][,1],col=2)
lines(results[[1]][,2],results[[1]][,1],col=3)

m=20
n=20
p=80
coe = c(0.916,0,0)
l=10
rho1 = 0.1
rho2= 0.1
prob=0.01
heter=F

load("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso/scripts/sample10_coe=0.5.Rd")
result10_0.5=results[[1]]
load("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso/scripts/sample20_coe=0.5.Rd")
result20_0.5=results[[1]]
load("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso/scripts/sample30_coe=0.5.Rd")
result30_0.5=results[[1]]


load("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso/scripts/sample10_coe=1.5.Rd")
result10_1.5=results[[1]]
load("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso/scripts/sample20_coe=1.5.Rd")
result20_1.5=results[[1]]
load("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso/scripts/sample30_coe=1.5.Rd")
result30_1.5=results[[1]]



load("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso/scripts/sample10_coe=2.5.Rd")
result10_2.5=results[[1]]
load("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso/scripts/sample20_coe=2.5.Rd")
result20_2.5=results[[1]]
load("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso/scripts/sample30_coe=2.5.Rd")
result30_2.5=results[[1]]



load("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso/scripts/lglasso210_coe=1.25.Rd")
result10_0.5=results[[1]]
load("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso/scripts/lglasso220_coe=1.25.Rd")
result20_0.5=results[[1]]
load("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso/scripts/lglasso230_coe=1.25.Rd")
result30_0.5=results[[1]]


