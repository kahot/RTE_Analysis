library(vegan)
library(fields)
library(MASS)
library(fossil)
library(mgcv)
library(DescTools)
library(ggplot2)
library(plotrix)
cols<-c("#ED8636","#0433FF")

#####  Using 12 biological replicates (tehc rep pooled)
rte_bio<-read.csv('Data/RTE.csv', header=T)
rte_bio1=t(rte_bio)
rte_bio1<-apply(rte_bio1,2,as.numeric)
rtebio1<-metaMDS(rte_bio1)
ordiplot(rtebio1,type= "text", display="sites" )

treatment<-c(rep("N->N",times=9), rep("N->O",times=9), rep("O->N",times=8), rep("O->O",times=9))
#rtedata<-cbind(treatment,rte_bio2)

NN<-c(1:9)
CN<-c(10:18)
NC<-c(19:26)
CC<-c(27:35)

nn=rtebio1$points[NN,]	
cn=rtebio1$points[CN,]	
nc=rtebio1$points[NC,]
cc=rtebio1$points[CC,]

pdf(file="Output/NMDSplots_ellipse.pdf", width=5.5, height=5)
ordiplot(rtebio1, type= "n", display="sites", xlim=c(-0.35, 0.35), ylim=c(-0.25, 0.3),xlab="NMDS Axis1", ylab="NMDS Axix2")  
points(nn, pch=21, col=cols[1], bg=cols[1], cex=2)
points(cn, pch=21, col=cols[2], bg=cols[1],cex=2,lwd=2)
points(nc, pch=21, col=cols[1],bg=cols[2],cex=2,lwd=2)
points(cc, pch=21, col=cols[2], bg=cols[2],cex=2)
draw.ellipse(x=-0.06, y=-0.11, a=.12, b=.12, angle=10, border=cols[1], lwd=1.5)
draw.ellipse(x=-0.015, y=0.096, a=.32, b=.08, angle=157, border=cols[2], lwd=1.5)

points(0.23,-0.15, pch=21,col=cols[1], bg=cols[1],cex=1.5)
text(0.29, -0.15, labels=expression("N" %->% "N"),cex=.8)
points(0.23,-0.175, pch=21,col=cols[2], bg=cols[1],cex=1.5)
text(0.29, -0.175, labels=expression("N" %->% "O"),cex=.8)
points(0.23,-0.20, pch=21,col=cols[2], bg=cols[2],cex=1.5)
text(0.29, -0.20, labels=expression("O" %->% "O"),cex=.8)
points(0.23,-0.225, pch=21,col=cols[1], bg=cols[2],cex=1.7)
text(0.29, -0.225, labels=expression("O" %->% "N"),cex=.8)

dev.off()

library(ggforce)
##ggplot
rte1 = as.data.frame(scores(rtebio1))
rte1$site<-c(rep("N", times=9),rep("O", times=9),rep("N", times=8),rep("O", times=9))
rte1$Origin<-c(rep("Nearshore", times=18),rep("Offshore", times=17))



ggplot(rte1, aes(x=NMDS1,y=NMDS2,color=Origin,shape=site)) + 
    geom_point(size=3.5) + theme_bw() +
    theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
    scale_shape_manual(values=c(19,17), name="Transplante site",labels=c("Nearshore", "Offshore"))+
    xlab("Axis 1") + ylab("Axis 2")+
    theme(axis.text=element_text(size=10))+
    scale_color_manual(values=paste0(cols,"CC"))+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))+
    geom_ellipse(aes(x0=-0.06, y0=-0.11, a=.12, b=.12, angle=10), color=cols[1])+
    geom_ellipse(aes(x0=0, y0=0.12, a=.40, b=0.06, angle=70), color=cols[2])+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=13))

ggsave("./Output/RTE_nmds.ellipse2.pdf", height = 5, width = 6.8)

