###In R###
data=read.table("all_fst_combined")
##variety treatment note
colnames(data)=c("pops","status","fst")
new_order <- with(data, reorder(pops , fst, mean , na.rm=T))
jpeg("fst_popwise_repeats_nonrpeeats_comparison.jpeg",width=28,height=18,units="in",res=300)
par(mar=c(5,6,4,1))
myplot <- boxplot(fst ~ status*new_order , data=data  , outline=F,
        boxwex=0.4,
        ann=F,
        col=c("seagreen4" , "coral3") , 
        border=c("seagreen3","coral2"),cex.axis=1.25,
        xaxt="n")
tt=(as.expression(bquote('F'['ST']*' comparison of LCR and non-LCR')))
title(main=tt ,line=1,cex.main=2.3)
title(ylab=as.expression(bquote('F'['ST'])), line=2, cex.lab=2.5)        
title(xlab="Populations", line=3, cex.lab=2.5)
my_names <- sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
my_names <- my_names[seq(1 , length(my_names) , 2)]
##change the number in seq below according to your number of populations. IT should be 2*(number of populations)
axis(1, 
     at = seq(1.5 , 20 , 2), 
     labels = my_names , 
     tick=FALSE , cex.axis=1.75)
     
# Add the grey vertical lines
##change the seq range (in our case it is 30) below according to your number of populations 
for(i in seq(0.5 , 30 , 2)){ 
  abline(v=i,lty=3, col="grey",lwd=1.75)
  }
 
# Add a legend ##Adjust the position of legend by changing x and y below
legend(x=0.4,y=0.04, legend = c("Fst of LCRs", "Fst of non-LCRs"), 
       col=c("seagreen4" , "coral3"),
       pch = 15, pt.cex = 6, cex = 2,  horiz = F, inset = c(0.1, 0.1),box.lwd = 0,box.col = "white",bg = "white")

###adding the significant values for each clade paired comparison
xpos=1
for(pop in unique(my_names)){
a1=data[data$pops==pop & data$status=="LCR",]
b1=data[data$pops==pop & data$status=="Non-LCR",]
yax=max(quantile(a1$fst,probs=0.75),quantile(b1$fst,probs=0.75))
pval=signif(as.numeric(wilcox.test(a1$fst,b1$fst,alternative="less")[3]),digits=3)
segments(xpos,yax+0.005,xpos+1,yax+0.005,lty=3,col="purple")
text(xpos+0.5,yax+0.0075,labels=pval,cex=1.5)
yminax=min(quantile(a1$fst,probs=0.25),quantile(b1$fst,probs=0.25))
text(xpos,-0.00075,bquote("n" == .(dim(a1)[1])),cex=1.1)
text(xpos+1,-0.00075,bquote("n" == .(dim(b1)[1])),cex=1.1)
xpos=xpos+2
}
dev.off()
