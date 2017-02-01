plhlist <- function(chil){
    pl <- NULL
    for(chiv in chil){
        if(!is.na(chiv)){
            if(chiv > 0){
                pl <- c(pl,1-pchisq(chiv,1))
            }
        }
    }
    return(pl)
}
dlfreqs <- function(dline){
    snum <- length(dline)/3
    freqs <- NULL
    for(i in c(1:snum)){
        freqs <- c(freqs,dline[3*i-2]/(dline[3*i-2]+dline[3*i-1]))
    }
    return(freqs)
}
dlgens <- function(dline){
    snum <- length(dline)/3
    gens <- NULL
    for(i in c(1:snum)){
        gens <- c(gens,dline[3*i])
    }
    return(gens)
}
matred <- function(mat,key,lw,up){
    mat <- mat[mat[,key]<up,]
    mat <- mat[mat[,key]>lw,]
    return(mat)
}
twolinecutter <- function(mat,lw2,up2,lw3,up3){
    mat <- matred(mat,2,lw2,up2)
    mat <- matred(mat,3,lw3,up3)
    return(mat[,1])
}
#main <- function(){
    range <- c(1:2000)
    dir <- "tmp/roc/sim_p1000_gen15,30,45,60_slc0.02_f100_dnum12_Dnum1000"
    lht <- read.table(paste(dir,"/rlt_wfe.dat",sep=""))[range,6]
    lhtnc <- read.table(paste(dir,"/rlt_wfenc.dat",sep=""))[range,6]
                                        #lhtsim <- read.table("tmp/fnum/rlt_s_p1000_gen60_slc0.00_f100_dnum6_Dnum100.dat")[,6]
    vlht <- read.table(paste(dir,"/rlt_wfe.dat",sep=""))[range,5]
    vlhtnc <- read.table(paste(dir,"/rlt_wfenc.dat",sep=""))[range,5]
bbg <- read.table(paste(dir,"/rlt_bbgp.dat",sep=""))[range,2]
bbgp <- bbg[1:1000]
bbg0 <- bbg[1001:2000]
    cmh <- read.table(paste(dir,"/rlt_cmh.dat",sep=""))[range,2]
    sim <- read.table(paste(dir,"/sim_s0.02.dat",sep=""))
    sim0 <- read.table(paste(dir,"/sim_s0.dat",sep=""))
                                        #lb <- read.table("tmp/mimicree/sim_p1000_s0.05_r6/lb_s0.05_r6.txt")[range,]
    lb <- c(rep(1,1000),rep(0,1000))
    plht <- plhlist(lht[lht<10])
ind <- c(1:2000)
pind <- c(1:1000)
nind <- c(1001:2000)
source("script/roc/rocr.r")
                                        #source("script/sync_process.r")
pdf("~/Documents/test.pdf")
                                        #roc
    prefl <- onerocr(lht,lb)
    prefln <- onerocr(lhtnc,lb)
    prefb <- onerocr(bbg,lb)
    prefc <- onerocr(cmh,lb)
    labels <- c("likelihood ratio","bbgp","cmh test")#,"wilcox test")
    xlim <- c(0,1)
    ylim <- c(0,1)
    ylab = "True positive rate"
    xlab = "False positive rate"
    cols <- c("black","red","blue")#,"green")
    par(col=cols[1])
    main <- "ROC curve"
    plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,main=main,cex=4,cex.lab=1.3,cex.main=2)
    #par(new="T",col=cols[1])
    #plot(prefl,ylab="",ann=F,lwd=1, xlim = xlim)
    par(new="T",col=cols[2])
    plot(prefb,ylab="",ann=F,lwd=1, xlim = xlim)
    par(new="T",col=cols[3])
    plot(prefc,ylab="",ann=F,lwd=1, xlim = xlim)
    par(new="T",col=cols[1])
    plot(prefln,ylab="",ann=F,lwd=1, xlim = xlim)
    par(new="T",col=cols[1])
    #lines(c(0,1),c(0,1))
    legend("bottomright", legend = labels, col = cols, lty = c(rep(1,4)),lwd=2,cex=1.5)                  
                                        #roc end
dev.off()
if(FALSE){
    #track each sample
brank <- rank(bbg[ind])
crank <- rank(cmh[ind])
lrank <- rank(lht[ind])
lncrank <- rank(lhtnc[ind])
lbr <- cbind(ind,lrank,brank)
hhind <- twolinecutter(lbr,1800,2000,1800,2000)
lhind <- twolinecutter(lbr,0,100,1400,1600)
hlind <- twolinecutter(lbr,1600,1800,600,900)
clr <- cbind(ind,crank,lncrank)
clhlind <- twolinecutter(clr,1400,1600,0,100)
cllhind <- twolinecutter(clr,1300,1600,1900,2000)
blr <- cbind(ind,brank,lncrank)
blhlind <- twolinecutter(blr,1500,1600,0,100)
bllhind <- twolinecutter(blr,0,100,1800,2000)
plind <- bllhind[1]
simfreqs <- dlfreqs(sim[plind,])
simgens <- dlgens(sim[plind,])
plot(simgens,simfreqs)
#track end
#rank 
brank <- rank(bbg[ind])
crank <- rank(cmh[ind])
lrank <- rank(lht[ind])
lncrank <- rank(lhtnc[ind])
plot(lncrank,brank)
points(lncrank[pind],brank[pind],col="red")
#rank end
    #bbg hist
lhtdt <- cbind(c(1:2000),lb,lhtnc)
cmhdt <- cbind(c(1:2000),lb,cmh)
bbgdt <- cbind(c(1:2000),lb,bbg)
sl <- sort(bbgdt[,3])
bbgdt <- bbgdt[sl,]
posar <- which(bbgdt[,2]==1)
hist(posar,breaks=20)
#bbg hist end
plot(lhtnc[1001:2000],lht[1001:2000],ylim=c(0,10),xlim=c(0,10))
        hist(plhtsim)
        hist(lht,xlim=c(0,5),breaks=1000)
        plot(lhnc,vlhnc,xlim=c(-10,40),ylim=c(-0.1,0.1))
    }    
#}
