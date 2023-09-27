a=read.table("flps_ucsc_bp_overlap",header=F)
b=aggregate(a$V3 ~ a$V1 + a$V2, data = a, FUN = sum, na.rm = TRUE)
write.table(b,file="first_aggregate.txt",append=F,sep='\t', quote = F,row.names = F, col.names = F)
