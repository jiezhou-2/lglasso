result=vector("list",20)
for (i in 0:19) {
  load(paste0("C:/Users/Jie Zhou/Desktop/lglassoNew/lglasso/scripts/lglasso_20_coe=2.5_tune=_",i,".Rd"))
  result[[i+1]]=results[[1]]
}

