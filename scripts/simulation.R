
load("./scripts/subject=15_time=20_coef=2.302_nodes=80.Rd")
results1520=simures
load("./scripts/subject=15_time=40_coef=2.302_nodes=80.Rd")
results1540=simures
load("./scripts/subject=30_time=20_coef=2.302_nodes=80.Rd")
results3020=simures
load("./scripts/subject=30_time=40_coef=2.302_nodes=80.Rd")
results3040=simures

par(mfrow=c(2,2),mar=c(4,4,2,2),oma=c(0,0,2,0))
FPR=results1520[[1]][[1]][,2]
TPR=results1520[[1]][[1]][,1]
plot(FPR,TPR,ylim = c(0.2,1),xlim=c(0,1),xlab = "FPR",ylab = "TPR",type="l")
lines(results1520[[1]][[2]][,2],results1520[[1]][[2]][,1],type="l",lty=2)
lines(results1520[[1]][[3]][,2],results1520[[1]][[3]][,1],type="l",lty=3)
legend(0.38,0.7,legend=c("m=15, n=20"))


FPR=results1540[[1]][[1]][,2]
TPR=results1540[[1]][[1]][,1]
plot(FPR,TPR,ylim = c(0.2,1),xlim=c(0,1),xlab = "FPR",ylab = "TPR",type="l")
lines(results1540[[1]][[2]][,2],results1520[[1]][[2]][,1],type="l",lty=2)
lines(results1540[[1]][[3]][,2],results1520[[1]][[3]][,1],type="l",lty=3)
legend(0.38,0.7,legend=c("m=15, n=40"))





FPR=results3020[[1]][[1]][,2]
TPR=results3020[[1]][[1]][,1]
plot(FPR,TPR,ylim = c(0.2,1),xlim=c(0,1),xlab = "FPR",ylab = "TPR",type="l")
lines(results3020[[1]][[2]][,2],results3020[[1]][[2]][,1],type="l",lty=2)
lines(results3020[[1]][[3]][,2],results3020[[1]][[3]][,1],type="l",lty=3)
legend(0.38,0.7,legend=c("m=30, n=20"))

FPR=results3040[[1]][[1]][,2]
TPR=results3040[[1]][[1]][,1]
plot(FPR,TPR,ylim = c(0.2,1),xlim=c(0,1),xlab = "FPR",ylab = "TPR",type="l")
lines(results3040[[1]][[2]][,2],results3040[[1]][[2]][,1],type="l",lty=2)
lines(results3040[[1]][[3]][,2],results3040[[1]][[3]][,1],type="l",lty=3)
legend(0.38,0.7,legend=c("m=30, n=40"))


title("Comparison of LGLASSO, GLASSO  and NH", outer=TRUE, cex=1.5)




