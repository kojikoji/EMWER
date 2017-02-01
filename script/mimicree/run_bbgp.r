
runbbgp <- function(filename){
    BFs <- NULL
    dtb <- read.table(filename)
    for(i in c(1:nrow(dtb))){
        rdl <- dtb[i,]
        dl <- NULL
        for(onrec in rdl){
            dl <- c(dl,strsplit(rdl,",")[[1]])
        }
        l <- length(dl)
        frq <- NULL
        dep <- NULL
        tms <- NULL
        for(j in c(1:as.integer(l/3))){
            frq <- c(frq,as.numeric(dl[3*j-2]))
            dep <- c(dep,as.numeric(dl[3*j-1]+dl[3*j-2]))
            tms <- c(tms,as.numeric(dl[3*j]))
        }
        print(frq)
        print(dep)
        print(tms)
        bb_model=betabinomialModel(frq,dep,tms)
        x=bb_model$timeVector
        y=bb_model$posteriorMean
        v=bb_model$posteriorVariance
        indModelCovTypes=c("white","fixedvariance")
        depModelCovTypes=c("rbf","white","fixedvariance")
        rslt=bbgp_test(x,y,v,indModelCovTypes,depModelCovTypes)
        BFs <- c(BFs,as.numeric(rslt$BF))
    }
    #setwd("../..")
    return(BFs)
}
print(getwd())
library(gptk)
setwd("script/BBGP")
source("loadBBGP.R")
loadBBGP()
setwd("../..")
args <- commandArgs(trailingOnly = T)
ifn <- args[1]
ofn <- args[2]
bf <- runbbgp(ifn)
print(bf)
write.table(bf,ofn, quote=F,col.names=F, append=F)
