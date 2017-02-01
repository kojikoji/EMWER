runcmhall <- function(inf){
    cmhs <- NULL
    dtb <- read.table(ifn)
    for(i in c(1:nrow(dtb))){
        dr <- dtb[i,]
        atbn <- 0
        ar <- NULL
        print(dr)
        for(dk in dr){
            dl <- as.numeric(strsplit(as.character(dk),",")[[1]])
            l <- length(dl)
        # number of row
            rwn <- as.integer(l/3)
                                        #make table (frq1 frq2 gen)
            frtb <- NULL
            for(j in c(1:as.integer((rwn-1)))){
                frtb <- c(frtb,as.numeric(dl[3*j-2]),as.numeric(dl[3*j-1]),as.numeric(dl[3*j+1]),as.numeric(dl[3*j+2]))
            }
                                        #print(frtb)
                                        #make array
                                        #number of table
            tbn <- rwn-1
            atbn <- atbn + tbn
            ar <- c(ar,frtb)
        }
        ar <- array(ar,dim=c(2,2,atbn))
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
cmhs <- runcmhall(ifn)
print(cmhs)
write.table(cmhs,ofn, quote=F,col.names=F, append=F)
