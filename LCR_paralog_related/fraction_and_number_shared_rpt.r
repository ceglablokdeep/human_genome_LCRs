a=read.table("hotspot_genes_lps")
a$V7="hotspot"
b=read.table("lps_paralog_without_hostpot_genes")
b$V7="no_hotspot"
a1=a[,c(6,7)]
b1=b[,c(6,7)]
data=rbind(a1,b1)
data$obs="Fraction"
data=data[,c(3,2,1)]
colnames(data)=c("obs","status","frc")
new_order <- with(data, reorder(obs , frc, mean , na.rm=T))
jpeg("fraction_shared_repeats_hotspot_nohotspot_comparison.jpeg",width=28,height=18,units="in",res=300)
par(mar=c(5,6,4,1))
myplot <- boxplot(frc ~ status*new_order , data=data  , outline=F,ylim=c(-0.1,1.1),
        boxwex=0.4,
        ann=F,
        col=c("seagreen4" , "coral3") , 
        border=c("seagreen3","coral2"),cex.axis=1.25,
        xaxt="n")
tt="Fraction of shared LCRs between paralogs in hotspot vs. no hotspot"
title(main=tt ,line=1,cex.main=2.3)
title(ylab="Fraction of shared LCRs", line=2, cex.lab=2.5)        
title(xlab="Paralogs", line=3, cex.lab=2.5)
my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
my_names <- my_names[seq(1 , length(my_names) , 2)]
axis(1, 
     at = 1.5, 
     labels = my_names , 
     tick=FALSE , cex.axis=1.75)
     
# Add a legend ##Adjust the position of legend by changing x and y below
legend(x=2,y=1, legend = c("Fraction of shared LCRs in hotspot", "Fraction of shared LCRs in no hotspot"), 
       col=c("seagreen4" , "coral3"),
       pch = 15, pt.cex = 6, cex = 2,  horiz = F, inset = c(0.1, 0.1),box.lwd = 0,box.col = "white",bg = "white")

###adding the significant values for each clade paired comparison
xpos=1
for(pop in unique(my_names)){
a1=data[data$obs==pop & data$status=="hotspot",]
b1=data[data$obs==pop & data$status=="no_hotspot",]
yax=max(quantile(a1$frc,probs=0.75),quantile(b1$frc,probs=0.75))
pval=signif(as.numeric(wilcox.test(a1$frc,b1$frc,alternative="greater")[3]),digits=3)
segments(xpos,yax+0.01,xpos+1,yax+0.01,lty=3,col="purple")
text(xpos+0.5,yax+0.05,labels=pval,cex=2.5)
yminax=min(quantile(a1$frc,probs=0.25),quantile(b1$frc,probs=0.25))
text(xpos,-0.075,bquote("n" == .(dim(a1)[1])),cex=2.1)
text(xpos+1,-0.075,bquote("n" == .(dim(b1)[1])),cex=2.1)
xpos=xpos+2
}
dev.off()
#####now for number of paralogs####
a1=a[,c(3,7)]
b1=b[,c(3,7)]
data=rbind(a1,b1)
data$obs="Number"
data=data[,c(3,2,1)]
colnames(data)=c("obs","status","num")
new_order <- with(data, reorder(obs , num, mean , na.rm=T))
jpeg("number_paralogs_hotspot_nohotspot_comparison.jpeg",width=28,height=18,units="in",res=300)
par(mar=c(5,6,4,1))
myplot <- boxplot(num ~ status*new_order , data=data  , outline=F,
        boxwex=0.4,
        ann=F,
        col=c("seagreen4" , "coral3") , 
        border=c("seagreen3","coral2"),cex.axis=1.25,
        xaxt="n")
tt="Number of paralogs between hotspot vs. no hotspot"
title(main=tt ,line=1,cex.main=2.3)
title(ylab="Number of paralogs", line=2, cex.lab=2.5)        
title(xlab="Paralogs", line=3, cex.lab=2.5)
my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
my_names <- my_names[seq(1 , length(my_names) , 2)]
axis(1, 
     at = 1.5, 
     labels = my_names , 
     tick=FALSE , cex.axis=1.75)
     
# Add a legend ##Adjust the position of legend by changing x and y below
legend(x=2,y=25, legend = c("Number of paralogs in hotspot", "Number of paralogs in no hotspot"), 
       col=c("seagreen4" , "coral3"),
       pch = 15, pt.cex = 6, cex = 2,  horiz = F, inset = c(0.1, 0.1),box.lwd = 0,box.col = "white",bg = "white")

###adding the significant values for each clade paired comparison
xpos=1
for(pop in unique(my_names)){
a1=data[data$obs==pop & data$status=="hotspot",]
b1=data[data$obs==pop & data$status=="no_hotspot",]
yax=max(quantile(a1$num,probs=0.75),quantile(b1$num,probs=0.75))
pval=signif(as.numeric(wilcox.test(a1$num,b1$num,alternative="greater")[3]),digits=3)
segments(xpos,yax+0.01,xpos+1,yax+0.01,lty=3,col="purple")
text(xpos+0.5,yax+1,labels=pval,cex=2.5)
yminax=min(quantile(a1$frc,probs=0.25),quantile(b1$frc,probs=0.25))
text(xpos,-0.5,bquote("n" == .(dim(a1)[1])),cex=2.1)
text(xpos+1,-0.5,bquote("n" == .(dim(b1)[1])),cex=2.1)
xpos=xpos+2
}
dev.off()
