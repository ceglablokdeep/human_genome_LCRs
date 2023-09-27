library(igraph)
a=read.table("ch17_hotspot_genes_lps")
b=read.table("paralog_in_ch17_cluster")
d=rbind(a,b)
g <- graph_from_data_frame(d)
d=d[order(d$V2), ]
d=d[!duplicated(d), ]
d=d[order(d$V2), ]
d$V3=c(rep(1,11),rep(2,22))
d$V4=c(rep(0.2,11),rep(0.7,22))
edgecolor=c("green4","indianred1")
colrs=c("steelblue3","grey50","grey50","grey50","grey50","grey50","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","steelblue3","grey50","grey50","grey50","hotpink2","hotpink2","hotpink2","hotpink2","hotpink2","hotpink2","hotpink2")
g <- graph_from_data_frame(a)
##edge.color => arrows ka color
##vertex.color => circle ka color
##vertex.frame.color => circle k boundary ka color
jpeg("genes_in_ch17_hotspot_rpts.jpeg",width=28,height=18,units="in",res=300)
plot(g,edge.color=edgecolor[as.numeric(d$V3)],vertex.color=colrs, vertex.frame.color=colrs,edge.width=as.numeric(d$V4)*4,vertex.label.color="black",vertex.label.cex=1.4)
tt="Genes in protein repeat hotspot in Chromosome 17 in Homo sapiens"
title(main=tt ,line=1,cex.main=2.3)
dev.off()
