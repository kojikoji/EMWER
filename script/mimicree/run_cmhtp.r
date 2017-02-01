runcmhtp <- function(inf){
    cmhs <- NULL
    dtb <- read.table(ifn)
    for(i in c(1:nrow(dtb))){
        dr <- dtb[i,]
        atbn <- 0
        ar <- NULL
        print(dr)
        for(dk in dr){
            print(dk)
            dl <- as.numeric(strsplit(as.character(dk),",")[[1]])
            print(dl)
            l <- length(dl)
        # number of row
            rwn <- as.integer(l/3)
                                        #make table (frq1 frq2 gen)
            frtb <- NULL
            for(j in c(1:rwn)){
                frtb <- rbind(frtb,c(as.numeric(dl[3*j-2]),as.numeric(dl[3*j-1]),as.numeric(dl[3*j])))
            }
                                        #print(frtb)
                                        #make array
            mngen <- min(frtb[,3])
            mxgen <- max(frtb[,3])
            mxfrtb <- as.matrix(frtb[frtb[,3]==mxgen,1:2])
            mnfrtb <- as.matrix(frtb[frtb[,3]==mngen,1:2])
            print(mxfrtb)
                                        #number of table
            tbn <- 1
            atbn <- atbn + tbn
            for(j in c(1:tbn)){
                ar <- c(ar,mxfrtb,mnfrtb)
            }
        }
        ar <- array(ar,dim=c(2,2,atbn))
        print(ar)
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
args <- commandArgs(trailingOnly = T)
if(TRUE){
ifn <- args[1]
ofn <- args[2]
}else{
ifn <- "tmp/test/sim_test.dat"
ofn <- "tmp/test/rlt_cmh.dat"
}     
cmhs <- runcmhtp(ifn)
print(cmhs)
write.table(cmhs,ofn, quote=F,col.names=F, append=F)
