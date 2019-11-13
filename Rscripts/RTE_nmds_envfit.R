library(vegan)
library(fields)
library(MASS)
library(fossil)
library(mgcv)
library(DescTools)
library(ggplot2)
library(plotrix)
library(grid)
library(ggforce)
cols<-c("#ED8636","#0433FF")

#####  Using 12 biological replicates (tehc rep pooled)
rte_bio<-read.csv('Data/RTE_BioRep12.csv', header=T)
rte_bio1=t(rte_bio)
rte_bio2<-rte_bio1[-c(1:2),]
rte_bio3<-apply(rte_bio2,2,as.numeric)
rtebio1<-metaMDS(rte_bio3)
ordiplot(rtebio1,type= "text", display="sites" )


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

##%%%%  Using Biological Process -separate detoxification into 2-Peroxiredoxins & the rest
#data formatting
biopro2<-read.csv("Data/RTE_GO_BP2.csv", header=F)
biopro2<-t(biopro2)
biopro2.1<-biopro2
colnames(biopro2.1) <- biopro2.1[1, ]
biopro2.1<-biopro2.1[-1,]
rownames(biopro2.1)<-biopro2.1[,1]
biopro2.1<-biopro2.1[,-1]

biopro2.1<-apply(biopro2.1,2,as.numeric)
biopro2.1<-data.frame(biopro2.1)


#run ENVFIT  
en1.2<-envfit(rtebio1,biopro2.1, choices=1:2, w=weights(rtebio1), na.rm = TRUE)
en1.2

#***VECTORS
#                                                  NMDS1    NMDS2     r2 Pr(>r)    
#Intracellular.protein.transport                 0.85755  0.51441 0.3742  0.117    
#Fatty.acid.beta.oxidation.using.ACAD            0.00386  0.99999 0.5857  0.020 *  
#Tricarboxylic.acid.cycle                       -0.69883  0.71528 0.1753  0.443    
#Purine.ribonuc..monophosphate.biosynth.        -0.93518  0.35417 0.1477  0.497    
#ATP.metabolic.process                          -0.95772  0.28772 0.4190  0.110    
#Arp2.3.complex.mediated.actin.nucleation       -0.30221  0.95324 0.0380  0.831    
#Nicotinamide.nucleotide.metabolic.process      -0.99312  0.11712 0.5553  0.027 *  
#tRNA.aminoacylation.for.protein.translation     0.76050  0.64934 0.8870  0.001 ***
#Oxidation.reduction.process                     0.28827  0.95755 0.8383  0.002 ** 
#Translation                                     0.99701  0.07724 0.5942  0.033 *  
#cellular.amino.acid.metabolic.process           0.58924  0.80796 0.9078  0.001 ***
#Coenzyme.biosynthetic.process                   0.25249  0.96760 0.1996  0.386    
#Oxidoreduction.coenzyme.metabolic.process      -0.98663  0.16297 0.5444  0.033 *  
#Chitin.catabolic.process                       -0.73732 -0.67554 0.1241  0.569    
#Hexose.metabolic.process                       -0.96773  0.25200 0.6698  0.003 ** 
#Lipid.homeostasis                               0.02946  0.99957 0.5870  0.021 *  
#Protein.folding                                -0.95276  0.30371 0.4731  0.064 .  
#Ribonucleoprotein.complex.assembly              0.96815  0.25035 0.7385  0.002 ** 
#Dicarboxylic.acid.metabolic.process            -0.98578  0.16804 0.1380  0.548    
#Alpha.amino.acid.biosynthetic.process           0.25618  0.96663 0.8833  0.001 ***
#Small.molecule.biosynthetic.process             0.23501  0.97199 0.4853  0.043 *  
#Carboxylic.acid.catabolic.process               0.56025  0.82833 0.3277  0.194    
#esicle.organization                            0.89112  0.45376 0.6769  0.003 ** 
#Translational.initiation                        0.06039 -0.99817 0.2562  0.273    
#Formation.of.translation.preinitiation.complex  0.32299  0.94640 0.2606  0.253    
#Cellular.oxidant.detoxification1                0.61843  0.78584 0.3263  0.166    
#Peroxiredoxins                                 -0.45663 -0.88966 0.6139  0.016 *  

#plot with ggplot
Arrows1.2<-data.frame(scores(en1.2, display="vectors"))
Arrows1.2$bp<-rownames(Arrows1.2)
Arrows1.3<-Arrows1.2
Arrows1.3[,1:2]<-Arrows1.2[,1:2]/3.5

ggplot() + geom_point(data=rte1, aes(x=NMDS1,y=NMDS2,color=Origin,shape=site),size=3.5) + theme_bw() +
    theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
    scale_shape_manual(values=c(19,17), name="Transplante site",labels=c("Nearshore", "Offshore"))+
    xlab("Axis 1") + ylab("Axis 2")+
    theme(axis.text=element_text(size=10))+
    scale_color_manual(values=paste0(cols,"CC"))+
    geom_segment(data=Arrows1.3, aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow= arrow(length=unit(0.2,"cm")), color="gray60", size=0.3)+
    geom_text(data=Arrows1.3,aes(x=NMDS1, y=NMDS2, label=bp), size=2, color="gray50" )+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))+
    geom_ellipse(aes(x0=-0.06, y0=-0.11, a=.12, b=.12, angle=10), color=cols[1])+
    geom_ellipse(aes(x0=0, y0=0.095, a=.33, b=0.06, angle=207), color=cols[2])+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=13))

ggsave("./Output/RTE_NMDS.envfit.bp.text.pdf", height = 5, width = 6.8)


##plot with numbers instead of BP names
Arrows2<-data.frame(scores(en1.2, display="vectors"))
Arrows2$bp<-1:nrow(Arrows2)
Arrows2[,1:2]<-Arrows2[,1:2]/3.5

ggplot() + geom_point(data=rte1, aes(x=NMDS1,y=NMDS2,color=Origin,shape=site),size=3.5) + theme_bw() +
    theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
    scale_shape_manual(values=c(19,17), name="Transplante site",labels=c("Nearshore", "Offshore"))+
    xlab("Axis 1") + ylab("Axis 2")+
    theme(axis.text=element_text(size=10))+
    scale_color_manual(values=paste0(cols,"CC"))+
    geom_segment(data=Arrows2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow= arrow(length=unit(0.2,"cm")), 
                 color="gray50", inherit_aes=FALSE)+
    geom_text(data=Arrows2,aes(x=NMDS1, y=NMDS2, label=bp), size=2, color="gray50" )+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))+
    geom_ellipse(aes(x0=-0.06, y0=-0.11, a=.12, b=.12, angle=10), color=cols[1])+
    geom_ellipse(aes(x0=0, y0=0.095, a=.33, b=0.06, angle=207), color=cols[2])+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=13))

ggsave("./Output/RTE_envfit.bp_numbers.pdf", height = 5, width = 6.8)

#no numbers/labels for the arrows
ggplot() + 
    geom_point(data=rte1, aes(x=NMDS1,y=NMDS2,color=Origin,shape=site),size=3.5) + theme_bw() +
    theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
    scale_shape_manual(values=c(19,17), name="Transplante site",labels=c("Nearshore", "Offshore"))+
    xlab("Axis 1") + ylab("Axis 2")+
    theme(axis.text=element_text(size=10))+
    scale_color_manual(values=paste0(cols,"CC"))+
    geom_segment(data=Arrows2, aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow= arrow(length=unit(0.2,"cm")), 
                 color="gray60", size=0.3)+
    geom_text(data=Arrows2,aes(x=NMDS1, y=NMDS2, label=""), size=2, color="gray60" )+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))+
    geom_ellipse(aes(x0=-0.06, y0=-0.11, a=.12, b=.12, angle=10), color=cols[1])+
    geom_ellipse(aes(x0=0, y0=0.095, a=.33, b=0.06, angle=207), color=cols[2])+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=13))

ggsave("./Output/RTE_envfit.bp_no-numbers.pdf", height = 5, width = 6.8)



##
#remove non-significant arrows
env1<-data.frame(scores(en1.2, display="vectors"))
Pvalue<-en1.2[[1]][[4]] #pvalues
env1<-cbind(env1, Pvalue)
BP<-rownames(env1)
env1<-data.frame(cbind(env1,BP) )
#for (i in 1:3) env1[,i]<-as.numeric(as.character(env1[,i]))
arrow3<-env1[env1$Pvalue<0.1,]            
arrow3[,1:2]<-arrow3[,1:2]/3.5
arrow3$no<-1:nrow(arrow3)

ggplot() + geom_point(data=rte1, aes(x=NMDS1,y=NMDS2,color=Origin,shape=site),size=3.5) + theme_bw() +
    theme(legend.text = element_text(size=12), legend.title = element_text(size=13))+
    scale_shape_manual(values=c(19,17), name="Transplante site",labels=c("Nearshore", "Offshore"))+
    xlab("Axis 1") + ylab("Axis 2")+
    theme(axis.text=element_text(size=10))+
    scale_color_manual(values=paste0(cols,"CC"))+
    geom_segment(data=arrow3, aes(x=0, xend=NMDS1, y=0, yend=NMDS2), arrow= arrow(length=unit(0.2,"cm")), 
                 color="gray50", inherit_aes=FALSE)+
    geom_text(data=arrow3,aes(x=NMDS1, y=NMDS2, label=BP), size=2, color="gray50" )+
    guides(shape = guide_legend(order = 2),color=guide_legend(order=1))+
    geom_ellipse(aes(x0=-0.06, y0=-0.11, a=.12, b=.12, angle=10), color=cols[1])+
    geom_ellipse(aes(x0=0, y0=0.095, a=.33, b=0.06, angle=207), color=cols[2])+
    theme(axis.text = element_text(size=12), axis.title = element_text(size=13))

ggsave("./Output/RTE_envfit.BP_sig.only.pdf", height = 5, width = 6.8)

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

