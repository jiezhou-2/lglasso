library(matrixcalc)
library(RCy3)
library(igraph)
library(fake)
library(tidyr)
library(devtools)
load_all()
load("/Users/f003r0s/Desktop/Anne_lab_projects/lglassoProject/data/b_set.Rdata")
variables=c("subject_accession","visit_name","feature","value_reported","infant_arm","delivery_mode","arm_name")
data=data_extra[,variables]
set.seed(1)

## data


 wideData=reshape(data=data,idvar =  c("subject_accession","visit_name","infant_arm","delivery_mode","arm_name"),
                  timevar = "feature", direction = "wide")

wideData=wideData[,-c(3:5)]
timepoints=ifelse(wideData$visit_name=="prevaccinated",2,ifelse(wideData$visit_name=="mo4",4,ifelse(wideData$visit_name=="vaccinated",5,9)))
wideData$visit_name=timepoints
wideData=wideData[order(wideData$subject_accession,wideData$visit_name),]


## homogeneous subjects
### two-stage model (full data)

index=which(is.na(wideData),arr.ind = T)
aa=as.numeric(names(which(table(index[,2])>0.1*nrow(wideData))))
fulldata=wideData[,-aa]
index=which(is.na(fulldata),arr.ind = T)
for (i in 1:nrow(index)) {
  fulldata[index[i,1],index[i,2]]=0
}
fulldata[,-c(1,2)]=log10(1+fulldata[,-c(1,2)])
index=which(fulldata$subject_accession==23)
fulldata=fulldata[-index,]
group=ifelse(fulldata$visit_name==2,1,2)

lambda1=exp(seq(-10,0,length=2))
lambda2=exp(seq(-10,0,length=2))
lambda=expand.grid(lambda1,lambda2)
aa1=CVlglasso(data=fulldata,lambda = lambda,K=2,trace = FALSE,random = FALSE,group = group)
saveRDS(aa1,file="homo2.rds")



