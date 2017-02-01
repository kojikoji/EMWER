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
eachmin <- function(x,y){
    z <- NULL
    for(i in c(1:length(x))){
        z <- c(z,min(x[i],y[i]))
    }
    return(z)
}
mkdist <- function(bvec){
    ct <- 0
    pos <- 0
    pos2 <- 0
    dist <- NULL
    tmpdist <- NULL
    for(v in bvec){
        pos <- pos + 1
        pos2 <- pos2 +1
        tmpdist <- c(tmpdist,pos2)
        if(v == 1){
            tmpdist2 <- pos2 - tmpdist
            tmpdist <- eachmin(tmpdist,tmpdist2)
            dist <- c(dist,tmpdist)
            tmpdist <- NULL
            pos2 <- 0
        }
    }
    dist <- c(dist,tmpdist)
    return(dist)
}
distsummer <- function(dist,vals){
    mxdist <- 120
    divnum <- 6
    wd <- mxdist/divnum
    ms <- NULL
    sds <- NULL
    nums <- NULL
    srs <- NULL
    for(i in c(1:divnum)){
        svals <- vals[wd*(i-1)<dist]
        dist2 <- dist[wd*(i-1)<dist]
        svals <- svals[dist2<wd*i]
        ms <- c(ms,mean(svals))
        sds <- c(sd,sd(svals))
        nums <- c(nums,length(svals))
        srs <- c(srs,length(svals[svals>4])/length(svals))
    }
    return(cbind(wd*c(1:divnum),ms,sds,nums,srs))
}
distdiv <- function(dist,vals){
    mxdist <- 100
    mt <- NULL
    for(ds in c(1:mxdist)){
        svals <- vals[dist==ds]
        mt <- cbind(mt,svals)
    }
    return(mt)
}
                                        #range <- c(1:4000)
                                        range <- c(4001:8000)
#range <- c(8001:12000)
lht <- read.table("tmp/mimicree/sim_p1000_s0.05_r6/rlt_wfe.dat")[range,6]
lhtnc <- read.table("tmp/mimicree/sim_p1000_s0.05_r6/rlt_wfenc.dat")[range,6]
lhtsim <- read.table("tmp/fnum/rlt_s_p1000_gen60_slc0.00_f100_dnum6_Dnum100.dat")[,6]
vlhtsim <- read.table("tmp/fnum/rlt_s_p1000_gen60_slc0.00_f100_dnum6_Dnum100.dat")[,5]
vlht <- read.table("tmp/mimicree/sim_p1000_s0.05_r6/rlt_wfe.dat")[range,5]
vlhtnc <- read.table("tmp/mimicree/sim_p1000_s0.05_r6/rlt_wfenc.dat")[range,5]
bbg <- read.table("tmp/mimicree/sim_p1000_s0.05_r6/rlt_bbgp.dat")[range,2]
lb <- read.table("tmp/mimicree/sim_p1000_s0.05_r6/lb_s0.05_r6.txt")[range,]
lbdst <- mkdist(lb)
vlhset <- cbind(vlht,lht)
vlhset <- cbind(c(1:4000),vlhset)
plht <- plhlist(lht[lht<10])
plhtsim <- plhlist(lhtsim[lhtsim<10])
dsum <- distsummer(lbdst,lht)
dsumnc <- distsummer(lbdst,lhtnc)
dmt <- distdiv(lbdst,lht)
source("script/roc/rocr.r")
pdf("~/Documents/test.pdf")
hist(plhtsim)
dev.off()
if(FALSE){
plot(dsum[,1],dsum[,5],ylim=c(0,1))
boxplot(dmt,ylim=c(-1,30))
plot(lbdst,lht,xlim=c(-1,30),ylim=c(-1,10))
plot(lht,vlht,xlim=c(-10,40),ylim=c(-0.1,0.1))
hist(lht,xlim=c(0,5),breaks=1000)
plot(lhnc,vlhnc,xlim=c(-10,40),ylim=c(-0.1,0.1))
#roc
prefl <- onerocr(lht,lb)
prefln <- onerocr(lhtnc,lb)
prefb <- onerocr(bbg,lb)
lb <- c("likelihood ratio","no control","bbgp","t test")#,"wilcox test")
xlim <- c(0,1)
ylim <- c(0,1)
ylab = "True positive rate"
xlab = "False positive rate"
cols <- c("black","red","blue")
par(col=cols[1])
main <- "ROC curve"
plot(0, 0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,main=main,cex=4,cex.lab=1.5,cex.main=2)
par(new="T",col=cols[1])
plot(prefl,ylab="",ann=F,lwd=1)
par(new="T",col=cols[2])
plot(prefln,ylab="",ann=F,lwd=1)
par(new="T",col=cols[3])
plot(prefb,ylab="",ann=F,lwd=1)
lines(c(0,1),c(0,1))
#roc end
}    
