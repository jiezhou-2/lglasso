

library(glasso)
library(lglasso)
library(BDgraph)
library(GGMselect)
#set.seed(1)
results=power_compare(m=20,n=10,p=80,coe = c(1.25,0,0),l=10,rho = 0.1,prob=0.01,heter=T)


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


