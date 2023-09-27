library(stringr)
a=read.table("repeats_coverage_50kb_windows",header=F)
b=read.table("exon_coverage_50kb_windows")
a$V5=a$V4/b$V4
a[is.na(a)] <- 0
###Selecting top 1% fraction
quantile(a$V5,na.rm = T,probs = 0.99)
c=a[a$V5>as.numeric(quantile(a$V5,na.rm = T,probs = 0.99)),]
write.table(c,file="hotspots_repeats_50kb.txt",append=T,sep='\t', quote = F,row.names = F, col.names = F)
##loading the GC% file
gc=read.table("human_50kb_gc_per",header=F)
quantile(gc$V4,prob=0.99)
#0.5639652
chsize=read.table("chromosomesizes")
pc_coor=read.table("human_protein_coding_genes_coordinates")
rpt_coor=read.table("flps_detected_single_chronly_exon_repeats_chr_rpts_combined_compact")
chrmsm=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")
ylimmax=(length(chrmsm)*4)+2
xlimmax=max(chsize$V2)+10000
##making a blank plot first
init=1
jpeg('karyoplot_human.jpeg',width=28,height=18,units="in",res=300)
plot(1, type = "n", xlab = "Chromosomes coordinates",ylab = "", xlim = c(0,xlimmax),ylim = c(0,ylimmax),yaxt="n",main="Hotspots of protein low-complexity regions (LCRs) in human genome",cex.main=3,cex.lab=2,cex.axis=2)
title(ylab="Chromosomes",line=2,cex.lab=2)
legend("topright",legend=c("Protein coding exons","Protein LCRs","Hotspots of protein LCRs","GC percentage"),col=c(rgb(0,0.1,0.7,alpha=0.3),"black","red","green2"),horiz=F,fill=c(rgb(0,0.1,0.7,alpha=0.3),"black","red","green"),inset=c(0.042,-0.01),bty="n",cex=2)
for(chr in chrmsm){
csz=as.numeric(chsize[chsize$V1==chr,2])
rect(0,init,csz,init+2)
gndf=pc_coor[pc_coor$V1==chr,]
rpdf=rpt_coor[rpt_coor$V1==chr,]
clusdf=c[c$V1==chr,]
gcdf=gc[gc$V1==chr,]
gcdf$V5=(gcdf$V2+gcdf$V3)/2
gcdf$V6=gcdf$V4+init+2
lines(gcdf$V5,gcdf$V6,type="l",col="green",lwd=1)
segments(0,init+2+0.5639652,csz,init+2+0.5639652,lty=3,lwd=0.7,col="purple")
for(pcrow in seq(1,dim(gndf)[1])){
pcseg=gndf[pcrow,]
rect(pcseg$V2,init,pcseg$V3,init+2,border="transparent",col=rgb(0,0.1,0.7,alpha=0.3))
}
if(dim(rpdf)[1]>0){
for(rprow in seq(1,dim(rpdf)[1])){
rpseg=rpdf[rprow,]
rect(rpseg$V2,init,rpseg$V3,init+2,border="transparent",col="black")
}
}
if(dim(clusdf)[1]>0){
for(clusrow in seq(1,dim(clusdf)[1])){
clusseg=clusdf[clusrow,]
rect(clusseg$V2-100,init-0.25,clusseg$V3+100,init+2.25,lty="dashed",col="red",lwd=0.5,border="red")
}
}
nch=str_remove(chr,"chr")
text(-5770000,init+1,labels=as.character(nch),cex=1.5)
init=init+4
}
dev.off()
