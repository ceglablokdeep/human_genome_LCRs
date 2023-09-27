gc=read.table("human_50kb_gc_per",header=F)
gc$V4=gc$V4*100
ec=read.table("exon_counts_50kb_window",header=F)
hc=read.table("hotspot_exon_count",header=F)
jpeg('gc_exon_count_human.jpeg',width=28,height=18,units="in",res=300)
plot(gc$V4,ec$V4,pch=20,col="rosybrown2",main="Distribution of exon count and GC% in 50kb window in human genome",xlab="GC percentage",ylab="",cex.main=3,cex.lab=2,cex.axis=2,cex=2)
title(ylab="Exon count",line=2.5,cex.lab=1.5)
for(hcrow in seq(1,dim(hc)[1])){
h1=hc[hcrow,]
h2=gc[gc$V1==as.character(h1$V1) & gc$V2==as.numeric(h1$V2) & gc$V3==as.numeric(h1$V3),]
h3=cbind(h2,h1$V4)
points(as.numeric(h3[1,4]),as.numeric(h3[1,5]),pch=18,col="red",cex=2)
}
legend("topright",legend="Hotspots of protein LCRs",col="red",horiz=F,fill="red",inset=c(0.042,-0.01),bty="n",cex=2)
dev.off()
