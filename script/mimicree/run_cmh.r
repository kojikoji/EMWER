runcmh <- function(inf,ifn0){
    cmhs <- NULL
    dtb <- read.table(ifn)
    dtb0 <- read.table(ifn0)
    for(i in c(1:nrow(dtb))){
        dl <- dtb[i,]
        dl0 <- dtb0[i,]
        l <- length(dl)
        ar <- NULL
        # nuber of table is 
        tbn <- as.integer(l/3)
        #array make
        for(j in c(1:tbn)){
            ar <- c(ar,as.numeric(dl[3*j-2]))
            ar <- c(ar,as.numeric(dl[3*j-1]))
            ar <- c(ar,as.numeric(dl0[3*j-2]))
            ar <- c(ar,as.numeric(dl0[3*j-1]))
        }
        ar <- array(ar,dim=c(2,2,tbn))
        rlt <- mantelhaen.test(ar,alternative=c("two.sided"))
        cmh <- matrix(rlt$statistic)[1,1]
        if(is.na(cmh)){
            cmh <- 0
        }
        cmhs <- c(cmhs,cmh)
    }
    #setwd("../..")
    return(cmhs)
}
print(getwd())
library(gptk)
setwd("script/BBGP")
source("loadBBGP.R")
loadBBGP()
setwd("../..")
args <- commandArgs(trailingOnly = T)
if(TRUE){
ifn <- args[1]
ofn <- args[2]
}else{
ifn <- "tmp/test/sim_test.dat"
ifn0 <- "tmp/test/simc_test.dat"
ofn <- "tmp/test/rlt_cmh.dat"
}     
cmhs <- runcmh(ifn)
print(cmhs)
write.table(cmhs,ofn, quote=F,col.names=F, append=F)
