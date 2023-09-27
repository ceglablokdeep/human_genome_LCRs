gns=read.table("ch17_cluster_genes",header=F)
exn=read.table("ch17window_exons",header=F)
rpt=read.table("ch17_clus_repeats",header=F)
gc=read.table("ch17_clus_gc",header=F)
jpeg('ch17_cluster_human.jpeg',width=48,height=18,units="in",res=300)
plot(1, type = "n", xlab = "Chromosomes coordinates",ylab = "", xlim = c(39990000,50010000),ylim = c(0,4),yaxt="n",main="Protein LCR hotspot of human genome chromosome 17",cex.main=3,cex.lab=2,cex.axis=2)
title(ylab="Chromosome 17",line=2,cex.lab=2)
rect(40000000,1,50000000,3)
gc$V5=(gc$V2+gc$V3)/2
gc$V6=gc$V4+3
gc=gc[order(gc$V2),]
lines(gc$V5,gc$V6,type="l",col="green",lwd=1)
segments(40000000,3+0.5639652,50000000,3+0.5639652,lty=3,lwd=0.7,col="purple")
for(exnrow in seq(1,dim(exn)[1])){
exnseg=exn[exnrow,]
rect(exnseg$V2,1,exnseg$V3,3,border="transparent",col=rgb(0,0.1,0.7,alpha=0.3))
}
for(gnrow in seq(1,dim(gns)[1])){
gnseg=gns[gnrow,]
rect(gnseg$V2,1,gnseg$V3,3,lty="dashed",border="orchid3",lwd=1)
text((gnseg$V2+gnseg$V3)/2,0.7,labels=as.character(gnseg$V4),cex=0.9,srt=90,adj=1)
}
for(rptrow in seq(1,dim(rpt)[1])){
rptseg=rpt[rptrow,]
rect(rptseg$V2,1,rptseg$V3,3,border="transparent",col="red")
}
dev.off()
