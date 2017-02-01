runcmhall <- function(inf){
    cmhs <- NULL
    dtb <- read.table(ifn)
    rwn <- 0
    for(i in c(1:nrow(dtb))){
        dr <- dtb[i,]
        atbn <- 0
        ar <- NULL
        print(dr)
        for(dk in dr){
            print(dk)
            dl <- as.numeric(strsplit(as.character(dk),",")[[1]])
            l <- length(dl)
        # number of row
            rwn <- as.integer(l/3)
            genl <- as.numeric(dl[seq(3,3*rwn,3)])
                                        #make table (frq1 frq2 gen)
            frtb <- NULL
            for(j in c(1:rwn)){
                frtb <- c(frtb,c(as.numeric(dl[3*j-2]),as.numeric(dl[3*j-1])))
            }
                                        #print(frtb)
                                        #make array
                                        #number of table
            tbn <- 1
            atbn <- atbn + tbn
            for(j in c(1:tbn)){
                ar <- c(ar,frtb)
            }
        }
        ar <- array(as.numeric(ar),dim=c(2,rwn,atbn)) + 1
        print(ar+1)
        rownames(ar) <- c("major","minor")
        colnames(ar) <- genl
        print(rownames(ar))
        print(colnames(ar))
        cmh <- CMHtest(ar, cscores="midrank",overall=TRUE)[["ALL"]][[1]][[1]]
        print("OK")
        if(is.na(cmh)){
            cmh <- 0
        }
        cmhs <- c(cmhs,cmh)
    }
    #setwd("../..")
    return(cmhs)
}
library(gptk)
library(vcdExtra)
args <- commandArgs(trailingOnly = T)
if(TRUE){
ifn <- args[1]
ofn <- args[2]
}else{
ifn <- "tmp/test/sim_test.dat"
ofn <- "tmp/test/rlt_cmh.dat"
}     
cmhs <- runcmhall(ifn)
write.table(cmhs,ofn, quote=F,col.names=F, append=F)
