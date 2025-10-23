rm(list=ls())
library(RColorBrewer)
library("Matrix")
library(pROC)
colors=c("#E04C5C","#7DAF4C","#23AECE", "#FB894B", "#E7DA36",  "#187A51",
         "#5EA4A2",  "#3D3C4E", "#4D1836", "#C51B7D",
         "#E9A3C9",  "#B35806", "#F1A340", "#FEE08B", "#D9EF8B",
         "#91CF60", "#C7EAE5", "#5AB4AC", "#01665E", "#E7D4E8",
         "#AF8DC3", "#762A83")


par(mfrow=c(1,2))
for (i in c(1,4)) {
  for (j in c(10)) {
k=100/j

      load(paste0("./output/figure1/simures",i, "_",j,"_",k,".Rd"))
      load(paste0("./output/figure1/K",i, "_",j,"_",k,".Rd"))

    lglassoPre=colSums(simures[[1]])/nrow(simures[[1]])
    glassoPre=colSums(simures[[2]])/nrow(simures[[2]])
    nhPre=colSums(simures[[3]])/nrow(simures[[3]])
    KK=K
    responce=KK[upper.tri(KK,diag = F)]
    responce_binary=ifelse(abs(responce)<=10^(-5),0,1)
    ##
    roc_nh=roc(responce_binary, nhPre)
    roc_glasso=roc(responce_binary,glassoPre)
    roc_lglasso=roc(responce_binary,lglassoPre)


    plot(roc_nh,col=colors[1])
    plot(roc_glasso,add=T,col=colors[2])
    plot(roc_lglasso,add=T,col=colors[3])
    legend(x=0.5,y=0.4, legend=c("NH", "GLASSO","LGLASSO"),
           col=c(colors[1], colors[2],colors[3]), lwd=2,cex = 0.8,box.lty = 0)
   # title(main = paste0("(",i,",",j,")"),adj = 0.5, line = 2.5)
   # title(main = paste0(i),adj = 0.5, line = 2.5)
#}
}
}



tar_load(comboIndex)
 tar_load(comboId)
 tar_load(preData)
 tar_load(postData)

 for (i in comboId) {
   print(i)
   index=comboIndex[i,]
aa=ff_lglasso(index=index,preData = preData,postData = postData,rho_lglasso = rho_lglasso)
 }
