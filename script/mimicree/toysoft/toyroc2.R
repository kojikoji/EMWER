library(ROCR)
spread<-c(0.0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010)


########################## BATCH 1 #####################################################

args<-commandArgs(TRUE)
alen=length(args)
filename=args[alen-1]
description=args[alen]

l1<-read.table(args[1])
f1<-read.table(args[2])
l2<-read.table(args[3])
f2<-read.table(args[4])


outfile1<-args[5]

prf1<-prediction(f1,l1)
prf2<-prediction(f2,l2)

pef1<-performance(prf1,"tpr","fpr")
pef2<-performance(prf2,"tpr","fpr")


pdf(file=paste(outfile1,".ps",sep=""),width=7,height=7)
par(mfrow=c(1,2))
plot(pef1,lwd=3,avg="vertical",spread.estimate="stderror",col="black")
plot(pef2,lwd=3,avg="vertical",spread.estimate="stderror",add=TRUE,col="red")

legend(0.0001,0.999,legend=c("s = 0.1","s = 0.025"),
       lty=c(1,1),lwd=c(2.5,2.5),col=c("black","red"))
dev.off()
