source("script/sync_process.r")
args <- commandArgs(trailingOnly = T)
ifn <- args[1]
genl <- as.numeric(strsplit(args[2],",")[[1]])
ofn <- args[3]
syncs <- read.table(ifn)
dt <- NULL
for(i in c(1:nrow(syncs))){
    if(length(synctodata(syncs[i,],genl)) > 30){
        print(syncs[i,])
        print(synctodata(syncs[i,],genl))
    }
    dt <- rbind(dt,synctodata(syncs[i,],genl))
}
write.table(dt,ofn,quote=FALSE,col.names=FALSE,row.names=FALSE)
    
