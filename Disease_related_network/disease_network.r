library(igraph)
a=read.table("new_filtered_interaction",fill=T)
b=read.table("disease_file_unique_paralog_numbered.txt")
rpts <- read.table("repeats.txt", header = FALSE, stringsAsFactors = FALSE)
combineListsAsOne=function(list1){
  n <- c()
  for(x in list1){
    n<-c(n, x)
  }
  return(n)
  }
diseases = unique(c(b$V3))
colors = rainbow(length(diseases))
a$V4 = ifelse(a$V3 == "", "grey", colors[match(a$V3,diseases)])
paras=unique(c(b$V4))
cols=rainbow(length(paras))
b$V5=ifelse(b$V4 == 0, "grey", cols[match(b$V4,paras)])
a2=a[a$V4!="grey" & a$V3!="INBORN_GENETIC_DISEASES",]
##I have to uniquely filter the dataframe
a2=unique(a2)
b2=b[b$V3!="INBORN_GENETIC_DISEASES",]
g2=graph_from_data_frame(a2,directed=F)

g2s <- simplify( g2, remove.multiple = F, remove.loops = T)
Isolated = which(degree(g2s)==0)
g2s1 = delete.vertices(g2s, Isolated)
g2s1df=as_data_frame(g2s1,what="vertices")
colv1=data.frame(v1=g2s1df$name, v2=b2[match(g2s1df$name, b2$V2), 5])
g2e1=as_data_frame(g2s1,what="edges")
g2e1=unique(g2e1)
dlab=unique(g2e1[,c(3,4)])
results <- list()
for (i in 1:length(g2s1df$name)) {
  gene = g2s1df$name[i]
  matches = b[b$V2 == gene, ]
  repeat_counts = lapply(rpts$V1, function(x) { sum(matches$V1 == x) })
  xyz=list(combineListsAsOne(repeat_counts))
  results[i] = xyz
}

colors_list <- c("red","green","blue","purple","orange","yellow","cyan","magenta","brown","pink","violet","turquoise","gold","sandybrown","maroon","plum","chartreuse","coral","slateblue1")
###finding the best seed number
for(s in 1:50){
set.seed(s)
naam=paste("LCR_genes_all_disorders_seed_",s,".jpeg",sep="")
jpeg(naam,width=48,height=30,units="in",res=300)
plot.igraph(g2s1,vertex.size=7,vertex.shape="pie", vertex.pie=results,vertex.pie.color=list(colors_list),vertex.label.cex=1.25,vertex.label.color="black",
    layout=layout.fruchterman.reingold(g2s1, niter=10000),edge.color=g2e1$V4,directed=F,edge.width =2, vertex.frame.color=colv1$v2)
legend(x=-1.2,y=1.1,legend=as.character(rpts$V1),col=colors_list,horiz=F,fill=colors_list,inset=c(0.042,-0.01),bty="n",cex=1.75)
legend(x=1.1,y=1.1,legend=as.character(dlab$V3),col=dlab$V4,horiz=F,fill=dlab$V4,inset=c(0.042,-0.01),bty="n",cex=1.2)
title("Interacting proteins with LCRs giving rise to genetic syndrome",cex.main=2.3,col.main="black")
dev.off()
}
####Take any one, they all are equally good and equally bad
######now without cancer also
a2=a[a$V4!="grey" & a$V3!="INBORN_GENETIC_DISEASES" & a$V3!="HEREDITARY_CANCER-PREDISPOSING_SYNDROME",]
##I have to uniquely filter the dataframe
a2=unique(a2)
b2=b[b$V3!="HEREDITARY_CANCER-PREDISPOSING_SYNDROME" & b$V3!="INBORN_GENETIC_DISEASES",]
g2=graph_from_data_frame(a2,directed=F)

g2s <- simplify( g2, remove.multiple = F, remove.loops = T)
Isolated = which(degree(g2s)==0)
g2s1 = delete.vertices(g2s, Isolated)
g2s1df=as_data_frame(g2s1,what="vertices")
colv1=data.frame(v1=g2s1df$name, v2=b2[match(g2s1df$name, b2$V2), 5])
g2e1=as_data_frame(g2s1,what="edges")
g2e1=unique(g2e1)
dlab=unique(g2e1[,c(3,4)])
results <- list()
for (i in 1:length(g2s1df$name)) {
  gene = g2s1df$name[i]
  matches = b[b$V2 == gene, ]
  repeat_counts = lapply(rpts$V1, function(x) { sum(matches$V1 == x) })
  xyz=list(combineListsAsOne(repeat_counts))
  results[i] = xyz
}
###finding the best seed number
for(s in 1:50){
set.seed(s)
naam=paste("LCR_genes_without_cancer_disorders_seed_",s,".jpeg",sep="")
jpeg(naam,width=30,height=18,units="in",res=300)
plot.igraph(g2s1,vertex.size=7,vertex.shape="pie", vertex.pie=results,vertex.pie.color=list(colors_list),
    layout=layout.fruchterman.reingold(g2s1, niter=10000),edge.color=g2e1$V4,directed=F,edge.width =2,vertex.label.cex=1.25, vertex.frame.color=colv1$v2,vertex.label.color="black")
legend(x=1.1,y=1.1,legend=as.character(dlab$V3),col=dlab$V4,horiz=F,fill=dlab$V4,inset=c(0.042,-0.01),bty="n",cex=0.85)
legend(x=-1.5,y=1.1,legend=as.character(rpts$V1),col=colors_list,horiz=F,fill=colors_list,inset=c(0.042,-0.01),bty="n",cex=1.5)
title("Interacting proteins with LCRs giving rise to genetic syndrome",cex.main=2.3,col.main="black")
dev.off()
}
#####now, only of cancer#########
a2=a[a$V4!="grey" & a$V3=="HEREDITARY_CANCER-PREDISPOSING_SYNDROME",]
##I have to uniquely filter the dataframe
a2=unique(a2)
b2=b[b$V3=="HEREDITARY_CANCER-PREDISPOSING_SYNDROME",]
g2=graph_from_data_frame(a2,directed=F)

g2s <- simplify( g2, remove.multiple = F, remove.loops = T)
Isolated = which(degree(g2s)==0)
g2s1 = delete.vertices(g2s, Isolated)
g2s1df=as_data_frame(g2s1,what="vertices")
colv1=data.frame(v1=g2s1df$name, v2=b2[match(g2s1df$name, b2$V2), 5])
g2e1=as_data_frame(g2s1,what="edges")
g2e1=unique(g2e1)
dlab=unique(g2e1[,c(3,4)])
results <- list()
for (i in 1:length(g2s1df$name)) {
  gene = g2s1df$name[i]
  matches = b[b$V2 == gene, ]
  repeat_counts = lapply(rpts$V1, function(x) { sum(matches$V1 == x) })
  xyz=list(combineListsAsOne(repeat_counts))
  results[i] = xyz
}

###finding the best seed number
for(s in 1:50){
set.seed(s)
naam=paste("LCR_genes_cancer_disorders_seed_",s,".jpeg",sep="")
jpeg(naam,width=28,height=18,units="in",res=300)
plot.igraph(g2s1,vertex.size=7,vertex.shape="pie", vertex.pie=results,vertex.pie.color=list(colors_list),
    layout=layout.fruchterman.reingold(g2s1, niter=10000),edge.color=g2e1$V4,directed=F,edge.width =2,vertex.label.cex=1.25, vertex.frame.color=colv1$v2,vertex.label.color="black")
legend(x=1.1,y=1.1,legend=as.character(dlab$V3),col=dlab$V4,horiz=F,fill=dlab$V4,inset=c(0.042,-0.01),bty="n",cex=1.5)
legend(x=-2,y=1.1,legend=as.character(rpts$V1),col=colors_list,horiz=F,fill=colors_list,inset=c(0.042,-0.01),bty="n",cex=2.25)
title("Interacting proteins with LCRs giving rise to Hereditary cancer-predisposing syndrome",cex.main=2.3,col.main="black")
dev.off()
}
