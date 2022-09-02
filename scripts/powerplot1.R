result=vector("list",45)
for (i in c(0:1,3:44)) {
  load(paste0("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso/scripts/lglasso_20_coe=2.5_tune=_",i,".Rd"))
  result[[i+1]]=results[[1]]
}

index1=c(0,8,10,13,15,16,21)
index2=0:44
index=setdiff(index2,index1)

for (i in index) {
  load(paste0("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso/scripts/lglasso_20_coe=0.4_tune=_",i,".Rd"))
  result[[i+1]]=results[[1]]
}



for (i in c(0:44)) {
  load(paste0("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso/scripts/lglasso_20_coe=1.5_tune=_",i,".Rd"))
  result[[i+1]]=results[[1]]
}

lglasso1=c()
lglasso2=c()
glasso1=c()
glasso2=c()
nh1=c()
nh2=c()
co1=c()
co2=c()
ew1=c()
ew2=c()
for (i in c(1,2,3:45)) {
  lglasso1=c(lglasso1,result[[i]][1,1])
  lglasso2=c(lglasso2,result[[i]][1,2])

  glasso1=c(glasso1,result[[i]][2,1])
  glasso2=c(glasso2,result[[i]][2,2])

  nh1=c(nh1,result[[i]][3,1])
  nh2=c(nh2,result[[i]][3,2])

  co1=c(co1,result[[i]][4,1])
  co2=c(co2,result[[i]][4,2])

  ew1=c(ew1,result[[i]][5,1])
  ew2=c(ew2,result[[i]][5,2])

}
lglasso=cbind(lglasso1,lglasso2)
lglasso=lglasso[order(lglasso[,2]),]
glasso=cbind(glasso1,glasso2)
glasso=glasso[order(glasso[,2]),]
nh=cbind(nh1,nh2)
nh=nh[order(nh[,2]),]
co=cbind(co1,co2)
co=co[order(co[,2]),]
ew=cbind(ew1,ew2)
ew=ew[order(ew[,2]),]
plot(lglasso[,2],lglasso[,1],type = "l")
lines(glasso[,2],glasso[,1])
lines(nh[,2],nh[,1])
lines(co[,2],co[,1])
lines(ew[,2],ew[,1])

plot(glasso[,2],glasso[,1],type = "l")
plot(nh[,2],nh[,1],type = "l")
plot(co[,2],co[,1],type = "l")
plot(ew[,2],ew[,1],type = "l")

