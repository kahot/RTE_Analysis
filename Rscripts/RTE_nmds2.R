library(vegan)
library(fields)
library(MASS)
library(fossil)
library(mgcv)
library(DescTools)
library(ggplot2)
library(plotrix)
library(grid)
cols<-c("#ED8636","#0433FF")

#####  Using 12 biological replicates (tehc rep pooled)
rte_bio<-read.csv('Data/RTE_BioRep12.csv', header=T)
rte_bio1=t(rte_bio)
rte_bio2<-rte_bio1[-c(1:2),]
rte_bio3<-apply(rte_bio2,2,as.numeric)
rtebio1<-metaMDS(rte_bio3)
ordiplot(rtebio1,type= "text", display="sites" )

treatment<-c("NN","NN","NN","ON","ON","ON","NO","NO","NO","OO","OO","OO")
#rtedata<-cbind(treatment,rte_bio2)

NN<-c(1,2,3)
CN<-c(4,5,6)
NC<-c(7,8,9)
CC<-c(10,11,12)

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
treatment<-c("NN","NN","NN","ON","ON","ON","NO","NO","NO","OO","OO","OO")
##ggplot
rte1 = as.data.frame(scores(rtebio1))
rte1$site<-c(rep("N", times=3),rep("O", times=3),rep("N", times=3),rep("O", times=3))
rte1$Origin<-c(rep("Nearshore", times=6),rep("Offshore", times=6))



ggplot(rte1, aes(x=NMDS1,y=NMDS2,color=Origin,shape=site)) + 
    geom_point(size=3.5) + theme_bw() +
    theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
    scale_shape_manual(values=c(19,17), name="Transplante site",labels=c("Nearshore", "Offshore"))+
    xlab("Axis 1") + ylab("Axis 2")+
    theme(axis.text=element_text(size=10))+
    scale_color_manual(values=paste0(cols,"CC"))+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))+
    geom_ellipse(aes(x0=-0.06, y0=-0.11, a=.12, b=.12, angle=10), color=cols[1])+
    geom_ellipse(aes(x0=0, y0=0.095, a=.33, b=0.06, angle=207), color=cols[2])+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=13))

ggsave("./Output/RTE_nmds.ellipse2.pdf", height = 5, width = 6.8)


    
######################################################## 
### Find Variables that expalin the ordination axis 
##use rtebio1 ordination results


##%%%%  1. Using Biological Process -NSAF values
###
#data formatting
biopro1<-read.csv("Data/RTE_GO_BP_nsaf.csv", header=F)
biopro1<-t(biopro1)
biopro1.1<-biopro1
colnames(biopro1.1) <- biopro1.1[1, ]
biopro1.1<-biopro1.1[-1,]
rownames(biopro1.1)<-biopro1.1[,1]
biopro1.1<-biopro1.1[,-1]

biopro1.1<-apply(biopro1.1,2,as.numeric)
biopro1.1<-data.frame(biopro1.1)


#run ENVFIT  
en1<-envfit(rtebio1,biopro1.1, choices=1:2, w=weights(rtebio1), na.rm = TRUE)
en1
#                                                  NMDS1    NMDS2     r2 Pr(>r)    
#Intracellular.protein.transport                 0.85792  0.51378 0.3743  0.120    
#Fatty.acid.beta.oxidation.using.ACAD            0.00432  0.99999 0.5855  0.028 *  
#Tricarboxylic.acid.cycle                       -0.69883  0.71529 0.1755  0.424    
#Purine.ribonuc..monophosphate.biosynth.        -0.93532  0.35382 0.1477  0.481    
#ATP.metabolic.process                          -0.95768  0.28784 0.4190  0.112    
#Arp2.3.complex.mediated.actin.nucleation       -0.30175  0.95339 0.0381  0.813    
#Nicotinamide.nucleotide.metabolic.process      -0.99307  0.11750 0.5552  0.018 *  
#tRNA.aminoacylation.for.protein.translation     0.76077  0.64902 0.8870  0.001 ***
#Oxidation.reduction.process                     0.28870  0.95742 0.8381  0.002 ** 
#Translation                                     0.99703  0.07700 0.5943  0.019 *  
#cellular.amino.acid.metabolic.process           0.58957  0.80772 0.9078  0.001 ***
#Cellular.oxidant.detoxification                -0.21997 -0.97551 0.5918  0.018 *  
#Coenzyme.biosynthetic.process                   0.25307  0.96745 0.1996  0.367    
#Oxidoreduction.coenzyme.metabolic.process      -0.98657  0.16335 0.5443  0.022 *  
#Chitin.catabolic.process                       -0.73726 -0.67561 0.1241  0.561    
#Hexose.metabolic.process                       -0.96760  0.25249 0.6697  0.005 ** 
#Lipid.homeostasis                               0.02993  0.99955 0.5868  0.026 *  
#Protein.folding                                -0.95267  0.30400 0.4730  0.066 .  
#Ribonucleoprotein.complex.assembly              0.96830  0.24978 0.7383  0.004 ** 
#Dicarboxylic.acid.metabolic.process            -0.98566  0.16875 0.1382  0.515    
#Alpha.amino.acid.biosynthetic.process           0.25653  0.96654 0.8833  0.001 ***
#Small.molecule.biosynthetic.process             0.23527  0.97193 0.4855  0.055 .  
#Carboxylic.acid.catabolic.process               0.56053  0.82813 0.3276  0.180    
#Vesicle.organization                            0.89128  0.45345 0.6772  0.005 ** 
#Translational.initiation                        0.05981 -0.99821 0.2563  0.264    
#Formation.of.translation.preinitiation.complex  0.32336  0.94627 0.2607  0.237   


NN<-c(1,2,3)
CN<-c(4,5,6)
NC<-c(7,8,9)
CC<-c(10,11,12)

nn=rtebio1$points[NN,]	
cn=rtebio1$points[CN,]	
nc=rtebio1$points[NC,]
cc=rtebio1$points[CC,]

ordiplot(rtebio1, type= "n", display="sites")  

ordiplot(rtebio1, type= "n", display="sites", xlim=c(-0.35, 0.35), ylim=c(-0.25, 0.3),xlab="NMDS Axis1", ylab="NMDS Axix2")  
points(nn, pch=21, col=cols[1], bg=cols[1], cex=2)
points(cn, pch=21, col=cols[2],bg=cols[2],cex=2,lwd=2)
points(nc, pch=21, col=cols[1],bg=cols[2],cex=2,lwd=2)
points(cc, pch=21, col=cols[2], bg=cols[1],cex=2)

plot(en1, cex=0.4, col='steelblue')


#plot with ggplot
Arrows<-data.frame(scores(en1, display="vectors"))
Arrows$bp<-rownames(Arrows)
Arrows2<-Arrows
Arrows2[,1:2]<-Arrows[,1:2]/3.5


########
ggplot() + 
    geom_point(data=rte1, aes(x=NMDS1,y=NMDS2,color=Origin,shape=site),size=3.5) + theme_bw() +
    theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
    scale_shape_manual(values=c(19,17), name="Transplante site",labels=c("Nearshore", "Offshore"))+
    xlab("Axis 1") + ylab("Axis 2")+
    theme(axis.text=element_text(size=10))+
    scale_color_manual(values=paste0(cols,"CC"))+
    geom_segment(data=Arrows2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow= arrow(length=unit(0.2,"cm")), 
                 color="gray50", inherit_aes=FALSE)+
    geom_text(data=Arrows2,aes(x=NMDS1, y=NMDS2, label=bp), size=2, color="gray50" )+
    #coord_fixed()+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))+
    geom_ellipse(aes(x0=-0.06, y0=-0.11, a=.12, b=.12, angle=10), color=cols[1])+
    geom_ellipse(aes(x0=0, y0=0.095, a=.33, b=0.06, angle=207), color=cols[2])+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=13))

ggsave("./Output/RTE_nmds.envfit.bp.pdf", height = 5, width = 6.8)

##plot with numbers instead of BP names
Arrows<-data.frame(scores(en1, display="vectors"))
Arrows$bp<-1:nrow(Arrows)
Arrows2<-Arrows
Arrows2[,1:2]<-Arrows[,1:2]/3.5



########
ggplot() + 
    geom_point(data=rte1, aes(x=NMDS1,y=NMDS2,color=Origin,shape=site),size=3.5) + theme_bw() +
    theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
    scale_shape_manual(values=c(19,17), name="Transplante site",labels=c("Nearshore", "Offshore"))+
    xlab("Axis 1") + ylab("Axis 2")+
    theme(axis.text=element_text(size=10))+
    scale_color_manual(values=paste0(cols,"CC"))+
    geom_segment(data=Arrows2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow= arrow(length=unit(0.2,"cm")), 
                 color="gray50", inherit_aes=FALSE)+
    geom_text(data=Arrows2,aes(x=NMDS1, y=NMDS2, label=bp), size=2, color="gray50" )+
    #coord_fixed()+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))+
    geom_ellipse(aes(x0=-0.06, y0=-0.11, a=.12, b=.12, angle=10), color=cols[1])+
    geom_ellipse(aes(x0=0, y0=0.095, a=.33, b=0.06, angle=207), color=cols[2])+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=13))

ggsave("./Output/RTE_nmds.envfit.bp_numbers1.pdf", height = 5, width = 6.8)

#no numbers/labels for the arrows
ggplot() + 
    geom_point(data=rte1, aes(x=NMDS1,y=NMDS2,color=Origin,shape=site),size=3.5) + theme_bw() +
    theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
    scale_shape_manual(values=c(19,17), name="Transplante site",labels=c("Nearshore", "Offshore"))+
    xlab("Axis 1") + ylab("Axis 2")+
    theme(axis.text=element_text(size=10))+
    scale_color_manual(values=paste0(cols,"CC"))+
    geom_segment(data=Arrows2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow= arrow(length=unit(0.2,"cm")), 
                 color="gray50", inherit_aes=FALSE)+
    geom_text(data=Arrows2,aes(x=NMDS1, y=NMDS2, label=""), size=2, color="gray50" )+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))+
    geom_ellipse(aes(x0=-0.06, y0=-0.11, a=.12, b=.12, angle=10), color=cols[1])+
    geom_ellipse(aes(x0=0, y0=0.095, a=.33, b=0.06, angle=207), color=cols[2])+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=13))

ggsave("./Output/RTE_nmds.envfit.bp_no-numbers.pdf", height = 5, width = 6.8)

###
##%%%%  Using Biological Process -separate detoxification into 2-Peroxiredoxins & the rest
###
#data formatting
biopro2<-read.csv("Data/RTE_GO_BP2.csv", header=F)
biopro21<-t(biopro1)
biopro21.1<-biopro1
colnames(biopro1.1) <- biopro1.1[1, ]
biopro1.1<-biopro1.1[-1,]
rownames(biopro1.1)<-biopro1.1[,1]
biopro1.1<-biopro1.1[,-1]

biopro1.1<-apply(biopro1.1,2,as.numeric)
biopro1.1<-data.frame(biopro1.1)


#run ENVFIT  
en1<-envfit(rtebio1,biopro1.1, choices=1:2, w=weights(rtebio1), na.rm = TRUE)
en1





##
#remove non-significant arrows
env1<-data.frame(scores(en1, display="vectors"))
Pvalue<-en1[[1]][[4]] #pvalues
env1<-cbind(env1, Pvalue)
BP<-rownames(env1)
env1<-data.frame(cbind(env1,BP) )
#for (i in 1:3) env1[,i]<-as.numeric(as.character(env1[,i]))
arrow2<-env1[env1$Pvalue<=0.05,]            
arrow2[,1:2]<-arrow2[,1:2]/3.5
arrow2$no<-1:nrow(arrow2)

ggplot() + geom_point(data=rte1, aes(x=NMDS1,y=NMDS2,color=Origin,shape=site),size=3.5) + theme_bw() +
    theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
    scale_shape_manual(values=c(19,17), name="Transplante site",labels=c("Nearshore", "Offshore"))+
    xlab("Axis 1") + ylab("Axis 2")+
    theme(axis.text=element_text(size=10))+
    scale_color_manual(values=paste0(cols,"CC"))+
    geom_segment(data=arrow2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow= arrow(length=unit(0.2,"cm")), 
                 color="gray50", inherit_aes=FALSE)+
    geom_text(data=arrow2,aes(x=NMDS1, y=NMDS2, label=BP), size=2, color="gray50" )+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))+
    geom_ellipse(aes(x0=-0.06, y0=-0.11, a=.12, b=.12, angle=10), color=cols[1])+
    geom_ellipse(aes(x0=0, y0=0.095, a=.33, b=0.06, angle=207), color=cols[2])+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=13))

ggsave("./Output/RTE_nmds.envfit.bp_sig.only.pdf", height = 5, width = 6.8)

ggplot() + geom_point(data=rte1, aes(x=NMDS1,y=NMDS2,color=Origin,shape=site),size=3.5) + theme_bw() +
    theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
    scale_shape_manual(values=c(19,17), name="Transplante site",labels=c("Nearshore", "Offshore"))+
    xlab("Axis 1") + ylab("Axis 2")+
    theme(axis.text=element_text(size=10))+
    scale_color_manual(values=paste0(cols,"CC"))+
    geom_segment(data=arrow2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow= arrow(length=unit(0.2,"cm")), 
                 color="gray50", inherit_aes=FALSE)+
    geom_text(data=arrow2,aes(x=NMDS1, y=NMDS2, label=""), size=2, color="gray50" )+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))+
    geom_ellipse(aes(x0=-0.06, y0=-0.11, a=.12, b=.12, angle=10), color=cols[1])+
    geom_ellipse(aes(x0=0, y0=0.095, a=.33, b=0.06, angle=207), color=cols[2])+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=13))

ggsave("./Output/RTE_nmds.envfit.bp_sig.only_no_text.pdf", height = 5, width = 6.8)


            
##%%%%  1. Using MF/CC -NSAF values
#data formatting
mfcc<-read.csv("Data/RTE_GO_MF_CC_nsaf.csv", header=F)
mfcc<-t(mfcc)
colnames(mfcc) = mfcc[1, ]
mfcc<-mfcc[-1,]
mfcc1.1<-mfcc[,-1]

mfcc1.1<-apply(mfcc1.1,2,as.numeric)
mfcc1.1<-data.frame(mfcc1.1)


#run ENVFIT  
en2<-envfit(rtebio1,mfcc1.1, choices=1:2, w=weights(rtebio1), na.rm = TRUE)
en2
 
#                                           NMDS1    NMDS2     r2 Pr(>r)   
#Serine.type.endopeptidase.inhib..activity  0.88119  0.47276 0.7190  0.008 **
#GTP.binding                                0.54660  0.83739 0.4715  0.058 . 
#Hydro.lyase.activity                      -0.82388 -0.56676 0.2703  0.238   
#Phosphotransferases                        0.99405  0.10891 0.0853  0.649   
#NAD.binding                               -0.61697  0.78698 0.5825  0.045 * 
#Oxidoreductase.activity..flavin.           0.00436  0.99999 0.5855  0.017 * 
#Pyridoxal.phosphate.binding                0.23828  0.97120 0.7017  0.004 **
#Transaminase.activity                      0.30186  0.95335 0.6536  0.005 **
#Translation.initiation.factor.activity     0.83173  0.55518 0.5193  0.037 * 
#Catalytic.step.2.spliceosome               0.78079  0.62480 0.5864  0.014 * 
#Clathrin.coat                              0.82207  0.56939 0.3896  0.087 . 
#Cytosolic.large.ribosomal.subunit          0.99472  0.10258 0.8172  0.003 **
#Endosome                                   0.74313  0.66915 0.5851  0.016 * 
#Eukaryotic.transl..ini..F3.complex         0.99999  0.00345 0.0034  0.980   
#Extracellular.exosome                      0.80955  0.58705 0.0400  0.815   
#Extracellular.region.part                 -0.72513 -0.68861 0.1476  0.482   
#Mitochondrion                              0.33388  0.94261 0.6205  0.006 **
#Oxidoreductase.complex                     0.81010  0.58629 0.0803  0.666   
#Proteasome.regulatory.particle..b.s..      0.41972  0.90765 0.6123  0.011 * 
#Transport.vesicle.membrane                -0.51713  0.85591 0.1473  0.500 


ordiplot(rtebio1, type= "n", display="sites", xlim=c(-0.35, 0.35), ylim=c(-0.25, 0.3),xlab="NMDS Axis1", ylab="NMDS Axix2")  
points(nn, pch=21, col='green', bg='green', cex=2)
points(cn, pch=21, col='#0000FF', bg='#00FF00',cex=2,lwd=2)
points(nc, pch=21, col='green',bg='blue',cex=2,lwd=2.5)
points(cc, pch=21, col='blue', bg='blue',cex=2)

plot(en2, cex=0.7, col='steelblue')


#remove non-significant arrows
env2<-en2[[1]][[1]]
Pvalue<-en2[[1]][[4]] #pvalues
env2<-cbind(env2, Pvalue)
GO<-rownames(env2)
env2<-data.frame(cbind(env2,GO) )
for (i in 1:3) env2[,i]<-as.numeric(as.character(env2[,i]))
arrow3<-env2[env2$Pvalue<=0.05,]            
arrow3[,1:2]<--arrow3[,1:2]/3.5
arrow3$no<-1:nrow(arrow3)

ggplot() + geom_point(data=rte1, aes(x=NMDS1,y=NMDS2,color=Origin,shape=site),size=3.5) + theme_bw() +
    theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
    scale_shape_manual(values=c(19,17), name="Transplante site",labels=c("Nearshore", "Offshore"))+
    xlab("Axis 1") + ylab("Axis 2")+
    theme(axis.text=element_text(size=10))+
    scale_color_manual(values=paste0(cols,"CC"))+
    geom_segment(data=arrow3, aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow= arrow(length=unit(0.2,"cm")), 
                 color="gray50", inherit_aes=FALSE)+
    geom_text(data=arrow3,aes(x=NMDS1, y=NMDS2, label=GO), size=2, color="gray50" )+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))+
    geom_ellipse(aes(x0=-0.06, y0=-0.11, a=.12, b=.12, angle=10), color=cols[1])+
    geom_ellipse(aes(x0=0, y0=0.095, a=.33, b=0.06, angle=207), color=cols[2])+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=13))

ggsave("./Output/RTE_nmds.envfit.bp_sig.only.pdf", height = 5, width = 6.8)

ggplot() + geom_point(data=rte1, aes(x=NMDS1,y=NMDS2,color=Origin,shape=site),size=3.5) + theme_bw() +
    theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
    scale_shape_manual(values=c(19,17), name="Transplante site",labels=c("Nearshore", "Offshore"))+
    xlab("Axis 1") + ylab("Axis 2")+
    theme(axis.text=element_text(size=10))+
    scale_color_manual(values=paste0(cols,"CC"))+
    geom_segment(data=arrow2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow= arrow(length=unit(0.2,"cm")), 
                 color="gray50", inherit_aes=FALSE)+
    geom_text(data=arrow2,aes(x=NMDS1, y=NMDS2, label=""), size=2, color="gray50" )+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))+
    geom_ellipse(aes(x0=-0.06, y0=-0.11, a=.12, b=.12, angle=10), color=cols[1])+
    geom_ellipse(aes(x0=0, y0=0.095, a=.33, b=0.06, angle=207), color=cols[2])+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=13))

ggsave("./Output/RTE_nmds.envfit.bp_sig.only_no_text.pdf", height = 5, width = 6.8)



#################################################

go<-read.csv("Data/RTE_GO_all.csv", header=F)
go<-t(go)
colnames(go) = go[1, ]
go<-go[-1,]
go1<-go[,-1]

go1<-apply(go1,2,as.numeric)
go1<-data.frame(go1)


#run ENVFIT  
en3<-envfit(rtebio1,go1, choices=1:2, w=weights(rtebio1), na.rm = TRUE)
en3
##plot with numbers instead of BP names
Arrows3<-data.frame(scores(en3, display="vectors"))
Arrows3$go<-rownames(Arrows3)
Arrows3.1<-Arrows3
Arrows3.1[,1:2]<-Arrows3[,1:2]/3.5



########
ggplot() + 
    geom_point(data=rte1, aes(x=NMDS1,y=NMDS2,color=Origin,shape=site),size=3.5) + theme_bw() +
    theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
    scale_shape_manual(values=c(19,17), name="Transplante site",labels=c("Nearshore", "Offshore"))+
    xlab("Axis 1") + ylab("Axis 2")+
    theme(axis.text=element_text(size=10))+
    scale_color_manual(values=paste0(cols,"CC"))+
    geom_segment(data=Arrows3.1, aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow= arrow(length=unit(0.2,"cm")), 
                 color="gray50", inherit_aes=FALSE)+
    geom_text(data=Arrows3.1,aes(x=NMDS1, y=NMDS2, label=go), size=2, color="gray50" )+
    #coord_fixed()+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))+
    geom_ellipse(aes(x0=-0.06, y0=-0.11, a=.12, b=.12, angle=10), color=cols[1])+
    geom_ellipse(aes(x0=0, y0=0.095, a=.33, b=0.06, angle=207), color=cols[2])+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=13))

ggsave("./Output/RTE_nmds.envfit.all.sig.GO.pdf", height = 5, width = 6.8)
