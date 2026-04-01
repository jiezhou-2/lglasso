###
edgeListHeter21=read.csv("./docs/adjMatrixPre.CSV")[,-1]
edgeListHeter22=read.csv("./docs/adjMatrixPost.CSV")[,-1]
##centrality
centralityPre=apply(edgeListHeter21, 1, sum)
centralityPost=apply(edgeListHeter22, 1, sum)
cenDiff=centralityPost-centralityPre
names(cenDiff)=colnames(edgeListHeter22)
##diff
edgeDiff=edgeListHeter22-edgeListHeter21
onlyPre=ifelse(edgeDiff==-1,1,0)
onlyPost=ifelse(edgeDiff==1,1,0)
write.csv(onlyPre,file="onlyPre.CSV")
write.csv(onlyPost,file="onlyPost.CSV")
changes=apply(edgeDiff, 1, sum)
names(changes)=colnames(edgeListHeter22)
## cluster analysis
clusterPre=read.csv("./docs/edgeListHeter21_cluster.csv")
clusterPre=split(clusterPre,f = clusterPre$Cluster)
clusterPost=read.csv("./docs/edgeListHeter22_cluster.csv")
clusterPost=split(clusterPost,f = clusterPost$Cluster)
pp=matrix(nrow = length(clusterPre),ncol = length(clusterPost))
for (i in 1:length(clusterPre)) {
  a1=clusterPre[[i]]$name
  for (j in 1:length(clusterPost)) {
b1=clusterPost[[j]]$name
pp[i,j]=length(intersect(a1,b1))/length(a1)
  }
}


ppp=matrix(nrow = length(clusterPre),ncol = length(clusterPost))
for (i in 1:length(clusterPre)) {
  a1=clusterPre[[i]]$name
  for (j in 1:length(clusterPost)) {
    b1=clusterPost[[j]]$name
    ppp[i,j]=length(intersect(a1,b1))/length(b1)
  }
}


pppp=matrix(nrow = length(clusterPre),ncol = length(clusterPost))
for (i in 1:length(clusterPre)) {
  a1=clusterPre[[i]]$name
  for (j in 1:length(clusterPost)) {
    b1=clusterPost[[j]]$name
    pppp[i,j]=length(intersect(a1,b1))
  }
}
colnames(pppp)=paste0("clusterPost",1:6)
rownames(pppp)=paste0("clusterPre",1:7)


cluPre1=clusterPre$Cluster1
cluPre2=clusterPre$Cluster2
cluPost1=clusterPost$Cluster1
cluPost2=clusterPost$Cluster2
index=which(cluPost1$name%in% cluPre2$name)
comp11=cluPost1$name[index]
comp12=cluPost1$name[-index]


index=which(cluPost2$name%in% cluPre1$name)
comp21=cluPost2$name[index]
comp22=cluPost2$name[-index]
