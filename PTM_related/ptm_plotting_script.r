##Plotting in R
a=read.table("genes_length")
b=read.table("lcrs_info")
b$st=(b$V2/b$V4)*100
b$en=(b$V3/b$V4)*100
c=read.table("ptms_info")
c$col="NA"
c$col[c$V4=="Serine_Phosphorylation"]="orchid3"
c$col[c$V4=="Threonine_Phosphorylation"]="brown2"
c$col[c$V4=="Tyrosine_Phosphorylation"]="dodgerblue4"
c$rp=(c$V3/c$V1)*100
jpeg("lcr_ptm_overlap.jpeg",width=20,height=11,units='in',res=300)
plot(0,xlim=c(0,102),ylim=c(0,28),type='n',xlab='',ylab='',axes=F,ann=F)
y=0
for (g in a$V1){
g1=as.character(g)
rect(0,y,100,y+1)
l=b[b$V1==g,]
for(z in 1:dim(l)[1]){
l1=l[z,]
rect(l1$st,y,l1$en,y+1,border="khaki2",lty=3,lwd=1,col=rgb(0.1,0,0.3,0.2))
legend(100.1,y+1.25,legend=g1,text.font=3,bty='n',adj=0,xpd=T)
ptm=c[c$V2==g,]
for(p in 1:dim(ptm)[1]){
p1=ptm[p,]
segments(p1$rp,y,p1$rp,y+1,lwd=2,col=as.character(p1$col))
}
}
y=y+2
}
legend(50,28,legend="Overlap of LCR and disease-susceptible PTM sites",adj=0.5,bty='n')
legend("bottomright",legend=c("LCR","Serine Phosphorylation","Threonine Phosphorylation","Tyrosine Phosphorylation"),fill=c(rgb(0.1,0,0.3,0.2),"orchid3","brown2","dodgerblue4"),col=c(rgb(0.1,0,0.3,0.2),"orchid3","brown2","dodgerblue4"),border=c(rgb(0.1,0,0.3,0.2),"orchid3","brown2","dodgerblue4"),bty='n',horiz=T)
dev.off()
