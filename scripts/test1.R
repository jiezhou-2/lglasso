id=unique(sample_data[,1])
x=matrix(rnorm(2*length(id)), nrow = length(id),ncol=2)
x=cbind(id,x)
a=lglasso(data = sample_data,x=x, rho = 0.7, heter=T, ty=1)
data = sample_data
x=x
rho = 0.7
heter=T
ty=1
