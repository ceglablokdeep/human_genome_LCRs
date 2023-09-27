library(igraph)
a=read.table("ch1_hotspot_genes")
b=read.table("paralog_in_ch1_cluster")
d=rbind(a,b)
g <- graph_from_data_frame(d)
d=d[order(d$V2), ]
d=d[!duplicated(d), ]
d=d[order(d$V2), ]
d$V3=c(rep(1,9),2,2,rep(3,32))
d$V4=c(rep(0.2,11),rep(0.7,32))
edgecolor=c("green4","olivedrab4","indianred1")
colrs=c("steelblue3","slateblue2","grey50","grey50","grey50","steelblue3","steelblue3","steelblue3","grey50","grey50","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","grey50","slateblue2","slateblue2","hotpink2","hotpink2","hotpink2","hotpink2","hotpink2","hotpink2","hotpink2")
##edge.color => arrows ka color
##vertex.color => circle ka color
##vertex.frame.color => circle k boundary ka color
jpeg("LCE_genes_in_ch1_hotspot_rpts.jpeg",width=28,height=18,units="in",res=300)
plot(g,edge.color=edgecolor[as.numeric(d$V3)],vertex.color=colrs, vertex.frame.color=colrs,edge.width=as.numeric(d$V4)*4,vertex.label.color="black",vertex.label.cex=1.4)
tt="Genes in protein repeat hotspot in Chromosome 1 in Homo sapiens"
title(main=tt ,line=1,cex.main=2.3)
dev.off()
