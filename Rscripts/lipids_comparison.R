#Lipid Content: 10.5.15
lipids<-read.csv("Data/lipids_after6.5.csv", header=T) 

library(Rcmdr)
leveneTest(Lipids ~ Site, data=lipids)
#      Df F value Pr(>F)
#group  1  2.2327 0.1524

leveneTest(Lipids ~ Treatment, data=lipids)
#      Df F value Pr(>F)
#group  3  1.0455 0.3994

a1 <- aov(Lipids~Treatment,data=lipids)
summary(a1)
#            Df  Sum Sq  Mean Sq F value Pr(>F)  
#Treatment    3 0.02502 0.008339   3.151 0.0539 .
#Residuals   16 0.04234 0.002646                 
plot(a1)
shapiro.test(lipids$Lipids) 
#data:  lipids$Lipids
#W = 0.92709, p-value = 0.1357

TukeyHSD(x=a1,conf.level=0.95)
#              diff         lwr         upr     p adj
#NO-NN  0.016651917 -0.07643267 0.109736503 0.9550729
#ON-NN  0.081475239 -0.01160935 0.174559824 0.0973067 
#OO-NN -0.008889362 -0.10197395 0.084195224 0.9925862
#ON-NO  0.064823321 -0.02826126 0.157907906 0.2315865
#OO-NO -0.025541279 -0.11862586 0.067543306 0.8601419
#OO-ON -0.090364600 -0.18344919 0.002719985 0.0586391  .

pairwise.t.test(lipids$Lipids,lipids$Treatment, p.adjust.method="none")
#   NN    NO    ON   
#NO 0.616 -     -    
#ON 0.023 0.064 -    
#OO 0.788 0.444 0.013

pairwise.t.test(lipids$Lipids,lipids$Treatment, p.adjust.method="fdr")

#   NN   NO   ON  
#NO 0.74 -    -   
#ON 0.07 0.13 -   
#OO 0.79 0.67 0.07

####
#### Site Effect on Genotype N ###
Ng<-lipids[1:10,]
a2 <- aov(Lipids~Site,data=Ng)
summary(a2)
#            Df   Sum Sq  Mean Sq F value Pr(>F)
#Site         1 0.000693 0.000693   0.214  0.656
#Residuals    8 0.025930 0.003241

#### Site Effect on Genotype C ###
Cg<-lipids[11:20,]
a2 <- aov(Lipids~Site,data=Cg)
summary(a2)
#            Df   Sum Sq  Mean Sq F value Pr(>F)
#Site         1 0.02041 0.020414   9.951 0.0135 *
#Residuals    8 0.01641 0.002051



#######################################################################
##################################
#Lipid Content: 10.5.15
lipids2<-read.csv("/Users/kahotisthammer/Documents/R/lipids.csv", header=T) 
str(lipids2)

library(Rcmdr)
leveneTest(Lipids ~ Site, data=lipids2)
#      Df F value Pr(>F)
#group  1  1.6757 0.2033
leveneTest(Lipids ~ Treatment, data=lipids2)
#       Df F value Pr(>F)
#group  1  0.0343  0.854
#      38
# % data but transformation not necessary based on Leven Test:


### N->N  comparing before and after###
lipid_N<-lipids2[1:10,]
#1-way anova
lipN1<-aov(Lipids~ Treatment, data=lipid_N )
summary(lipN1)
op<-par(mfrow=c(2,2)); plot(lipN1); par(op)
#            Df  Sum Sq  Mean Sq F value Pr(>F)
#Treatment    1 0.00277 0.002765    0.16    0.7
#Residuals    8 0.13868 0.017335

### N->O before and after ###
lipid_N2<-lipids2[11:20,]
#1-way anova
lipN2<-aov(Lipids~ Treatment, data=lipid_N2 )
summary(lipN2)
#            Df  Sum Sq  Mean Sq F value Pr(>F)
#Treatment    1 0.00788 0.007882   0.424  0.533
#Residuals    8 0.14862 0.018578
op<-par(mfrow=c(2,2)); plot(lipN2); par(op)


### O->N  before and after ###
lipid_O<-lipids2[21:30,]
#1-way anova
lipO1<-aov(Lipids~ Treatment, data=lipid_O )
summary(lipO1)
#            Df  Sum Sq Mean Sq F value  Pr(>F)   
#Treatment    1 0.14075 0.14075   22.22 0.00151 **
#Residuals    8 0.05068 0.00633

leveneTest(Lipids ~ Treatment, data=lipid_O)
#       Df F value  Pr(>F)  
#group  1  5.1205 0.05348 .
#       8 

### O->O  before and after ###
lipid_O2<-lipids2[31:40,]
#1-way anova
lipO2<-aov(Lipids~ Treatment, data=lipid_O2 )
summary(lipO2)
#            Df   Sum Sq Mean Sq F value   Pr(>F)    
#Treatment    1 0.023793 0.02379   29.75 0.000605 ***
#Residuals    8 0.006398 0.00080 

op<-par(mfrow=c(2,2)); plot(lipO2); par(op)

leveneTest(Lipids ~ Treatment, data=lipid_O2)
#      Df F value Pr(>F)
#group  1  1.1948 0.3062


# Extracting before exp values
nb<-lipids2[1:5,]
nc<-lipids2[21:25,]
before<-rbind(nb,nc)

lip.b<-aov(Lipids~ Origin, data=before )
summary(lip.b)
#             Df  Sum Sq Mean Sq F value Pr(>F)
#Origin       1 0.01128 0.01128   0.781  0.403
#Residuals    8 0.11554 0.01444

#######################################################################
library(plotrix)
library(ggplot2)
library(gridExtra)
library(ggthemes)
cols<-c("#ED8636","#0433FF")
lipids<-read.csv("Data/lipids_after6.5.csv", header=T) 

Ncoral<-subset(lipids, Origin=="N")
Ocorals<-subset(lipids, Origin=="O")

Nmean<-tapply(Ncoral$Lipids, Ncoral$Treatment,mean)
Omean<-tapply(Ocorals$Lipids, Ocorals$Treatment,mean)

nn<-subset(lipids, Treatment=="NN")
no<-subset(lipids, Treatment=="NO")
on<-subset(lipids, Treatment=="ON")
oo<-subset(lipids, Treatment=="OO")
Origin<-c("Nearshore","Nearshore","Offshore","Offshore")

stats<-data.frame(Site=c("Nearshore","Offshore","Nearshore","Offshore"),
                  Origin=c("Nearshore","Nearshore","Offshore","Offshore"),
                  Type=c("NN","NO","ON","OO"), 
                  Meanlipids=c(mean(nn$Lipids)*100, mean(no$Lipids)*100,mean(on$Lipids)*100,mean(oo$Lipids)*100), 
                  SE=c(std.error(nn$Lipids*100) ,std.error(no$Lipids*100),std.error(on$Lipids*100),std.error(oo$Lipids*100)))

lipid1<- ggplot(stats, aes(x=Site,y=Meanlipids,color=factor(Origin)))+
  geom_errorbar(mapping=aes(ymin=Meanlipids-SE, ymax=Meanlipids+SE,color=factor(Origin)), width=.1)+
  scale_color_manual("Origin", breaks=c(1,2), values=cols)+
  geom_point(aes(fill=factor(Origin)),size=5)+
  scale_fill_manual("Origin", breaks=c(1,2), values=cols)+
  labs(x="Transplant site",y="Tissue lipid content (w/w%)") +theme_base()+theme(legend.position="none")+
geom_segment(data=stats, mapping=aes(x=1, y=stats[1,4], xend=1.9, yend=stats[2,4]), size=.8, arrow=arrow(length = unit(0.3, "cm")),color=cols[1])+
geom_segment(data=stats, mapping=aes(x=2, y=stats[4,4], xend=1.1, yend=stats[3,4]), size=.8, arrow=arrow(length = unit(0.3, "cm")), color=cols[2])
lipid1
ggsave(filename="~/Dropbox/R/RTE_Analysis/Lipids1.pdf",width=4, height=3.8, units='in',device='pdf', plot=lipid1 )
library(cowplot)
#no arrows 
lipid2<-ggplot(stats, aes(x=Site,y=Meanlipids,color=factor(Origin)))+
        geom_errorbar(mapping=aes(ymin=Meanlipids-SE, ymax=Meanlipids+SE,color=factor(Origin)), width=.17)+
        scale_color_manual("Origin", breaks=c(1,2), values=cols)+
        geom_point(aes(fill=factor(Origin)),size=5)+
        scale_fill_manual("Origin", breaks=c(1,2), values=cols)+
        labs(x="Transplant site",y="Tissue lipid content (w/w%)") +theme_bw()+theme(legend.position="none")+
        geom_segment(data=stats, mapping=aes(x=1, y=stats[1,4], xend=2, yend=stats[2,4]), size=.7, color=cols[1])+
        geom_segment(data=stats, mapping=aes(x=2, y=stats[4,4], xend=1, yend=stats[3,4]), size=.7,  color=cols[2])+
        theme(axis.text.x =element_text(size=12, color="black"), panel.grid.major.x = element_blank(),
              axis.title = element_text(size=14))


lipid2<-lipid2+annotate(geom="text", x=1.1, y=28.5, label=paste0("a^{.}~(italic(P) == 0.059)"),color ='black', size=4,parse=TRUE,hjust=0)+
        annotate(geom="text", x=1.1, y=20.7, label="a",color ='black', size=4)+
        annotate(geom="text", x=2.1, y=22.5, label="a",color ='black', size=4)+
        annotate(geom="text", x=2.1, y=19.5, label="a",color ='black', size=4)

lipid2        
ggsave(filename="Output/Lipids1.pdf",width=4, height=3.8, units='in',device='pdf', plot=lipid2 )


## From Tissue_layer.R 
pdf("Output/Tissue_Lipids.pdf",width=8, height=3.8)
grid.arrange(plot2,lipid2,nrow=1)
dev.off()

#
p<-arrangeGrob(plot2,lipid2,nrow=1)
p
ggsave(filename="Output/Tissue_Lipids3.pdf",width=8, height=3.8, units='in',device='pdf', plot=p )
               


####### bwfore and after figures########
lipids2<-read.csv("Data//lipids.csv", header=T) 
Ncoral<-subset(lipids2, Origin=="N")
Ocoral<-subset(lipids2, Origin=="O")


Nstats<-data.frame(name=c("N-before","N->N","N->O"),treatment=c("1.Before","2.After","2.After"),
                   Meanlipids=c(mean(Ncoral$Lipids[Ncoral$Treatment=="Before"&Ncoral$Site=="N"])*100,
                                mean(Ncoral$Lipids[Ncoral$Treatment=="After"& Ncoral$Site=="N"])*100,
                                mean(Ncoral$Lipids[Ncoral$Treatment=="After"& Ncoral$Site=="O"])*100),
                   SE=c(std.error(Ncoral$Lipids[Ncoral$Treatment=="Before"&Ncoral$Site=="N"])*100,
                        std.error(Ncoral$Lipids[Ncoral$Treatment=="After"& Ncoral$Site=="N"])*100,
                        std.error(Ncoral$Lipids[Ncoral$Treatment=="After"& Ncoral$Site=="O"])*100)
                   )

library(grid)
text_before <- textGrob("Before", gp=gpar(fontsize=11))
text_after <- textGrob("After", gp=gpar(fontsize=11))


N_plot<-ggplot(Nstats, aes(x=treatment,y=Meanlipids,color=factor(name)))+
        scale_x_discrete(breaks=1:2,labels=c("Before", "After"))+
        geom_segment(aes(x=1, y=Nstats[1,3]-Nstats[1,4],       xend=1, yend=Nstats[1,3]+Nstats[1,4]), size=.4,color=cols[1])+
        geom_segment(aes(x=1.95, y=Nstats[2,3]-Nstats[2,4], xend=1.95, yend=Nstats[2,3]+Nstats[2,4]), size=.4,color=cols[1])+
        geom_segment(aes(x=2   , y=Nstats[3,3]-Nstats[3,4],    xend=2, yend=Nstats[3,3]+Nstats[3,4]), size=.4,color=cols[1])+
        geom_segment(aes(x=0.95, y=Nstats[1,3]-Nstats[1,4], xend=1.05, yend=Nstats[1,3]-Nstats[1,4]), size=.4,color=cols[1])+
        geom_segment(aes(x=0.95, y=Nstats[1,3]+Nstats[1,4], xend=1.05, yend=Nstats[1,3]+Nstats[1,4]), size=.4,color=cols[1])+
        geom_segment(aes(x=1.9, y=Nstats[2,3]-Nstats[2,4],     xend=2, yend=Nstats[2,3]-Nstats[2,4]), size=.4,color=cols[1])+
        geom_segment(aes(x=1.9, y=Nstats[2,3]+Nstats[2,4],     xend=2, yend=Nstats[2,3]+Nstats[2,4]), size=.4,color=cols[1])+
        geom_segment(aes(x=1.95, y=Nstats[3,3]-Nstats[3,4], xend=2.05, yend=Nstats[3,3]-Nstats[3,4]), size=.4,color=cols[1])+
        geom_segment(aes(x=1.95, y=Nstats[3,3]+Nstats[3,4], xend=2.05, yend=Nstats[3,3]+Nstats[3,4]), size=.4,color=cols[1])+
        labs(x="",y="Tissue lipid content (w/w %)") +
        theme_bw()+theme(legend.position="none")+
        geom_segment(data=Nstats, mapping=aes(x=1, y=Nstats[1,3], xend=1.9, yend=Nstats[2,3]), size=.5, arrow=arrow(length = unit(0.3, "cm")),color=cols[1])+
        geom_segment(data=Nstats, mapping=aes(x=1, y=Nstats[1,3], xend=1.9, yend=Nstats[3,3]), size=.5, arrow=arrow(length = unit(0.3, "cm")), color=cols[2])+
        geom_point(aes(x=1, y=Nstats[1,3]),col=cols[1], bg=cols[1], pch=21, size=3)+
        geom_point(aes(x=1.95, y=Nstats[2,3]),col=cols[1],bg=cols[1], pch=21, size=3)+
        geom_point(aes(x=2   , y=Nstats[3,3]),col="#0433FF",bg=cols[1], pch=21, size=3)+
        theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
        ylim(0,33)+
        theme(plot.margin = unit(c(1,1,2,1), "lines"))+
        annotation_custom(text_before, xmin=1, xmax=1, ymin=-5, ymax=-5)+
        annotation_custom(text_after, xmin=2, xmax=2, ymin=-5, ymax=-5)+
        annotate("point", x = 2, y =5 , colour = cols[1], size=3)+annotate(geom="text", x=2.3, y=5, label=expression("N" %->% "N"),size=3)+
        annotate("point", x = 2, y =2.5 , col="#0433FF",bg=cols[1], pch=21, size=3)+annotate(geom="text", x=2.3, y=2.5, label=expression("N" %->% "O"),size=3)+
        coord_cartesian(clip="off")
N_plot

## 
Ostats<-data.frame(name=c("O-before","O->O","O->N"),treatment=c("1.Before","2.After","2.After"),
                   Meanlipids=c(mean(Ocoral$Lipids[Ocoral$Treatment=="Before"&Ocoral$Site=="O"])*100,
                                mean(Ocoral$Lipids[Ocoral$Treatment=="After"& Ocoral$Site=="O"])*100,
                                mean(Ocoral$Lipids[Ocoral$Treatment=="After"& Ocoral$Site=="N"])*100),
                   SE=c(std.error(Ocoral$Lipids[Ocoral$Treatment=="Before"&Ocoral$Site=="O"])*100,
                        std.error(Ocoral$Lipids[Ocoral$Treatment=="After"& Ocoral$Site=="O"])*100,
                        std.error(Ocoral$Lipids[Ocoral$Treatment=="After"& Ocoral$Site=="N"])*100))

O_plot<-ggplot(Ostats, aes(x=treatment,y=Meanlipids,color=factor(name)))+
        scale_x_discrete(breaks=1:2,labels=c("Before", "After"))+xlab(NULL)+
        geom_segment(aes(x=1,    y=Ostats[1,3]-Ostats[1,4], xend=1,    yend=Ostats[1,3]+Ostats[1,4]), size=.4,color=cols[2])+
        geom_segment(aes(x=1.95, y=Ostats[2,3]-Ostats[2,4], xend=1.95, yend=Ostats[2,3]+Ostats[2,4]), size=.4,color=cols[2])+
        geom_segment(aes(x=2   , y=Ostats[3,3]-Ostats[3,4], xend=2,    yend=Ostats[3,3]+Ostats[3,4]), size=.4,color=cols[2])+
        geom_segment(aes(x=0.95, y=Ostats[1,3]-Ostats[1,4], xend=1.05, yend=Ostats[1,3]-Ostats[1,4]), size=.4,color=cols[2])+
        geom_segment(aes(x=0.95, y=Ostats[1,3]+Ostats[1,4], xend=1.05, yend=Ostats[1,3]+Ostats[1,4]), size=.4,color=cols[2])+
        geom_segment(aes(x=1.9,  y=Ostats[2,3]-Ostats[2,4], xend=2,    yend=Ostats[2,3]-Ostats[2,4]), size=.4,color=cols[2])+
        geom_segment(aes(x=1.9,  y=Ostats[2,3]+Ostats[2,4], xend=2,    yend=Ostats[2,3]+Ostats[2,4]), size=.4,color=cols[2])+
        geom_segment(aes(x=1.95, y=Ostats[3,3]-Ostats[3,4], xend=2.05, yend=Ostats[3,3]-Ostats[3,4]), size=.4,color=cols[2])+
        geom_segment(aes(x=1.95, y=Ostats[3,3]+Ostats[3,4], xend=2.05, yend=Ostats[3,3]+Ostats[3,4]), size=.4,color=cols[2])+
        labs(x="",y="") +
        theme_light()+theme(legend.position="none")+
        theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
        geom_segment(data=Ostats, mapping=aes(x=1, y=Ostats[1,3], xend=1.9, yend=Ostats[2,3]), size=.5, arrow=arrow(length = unit(0.3, "cm")), color="#0433FF")+
        geom_segment(data=Ostats, mapping=aes(x=1, y=Ostats[1,3], xend=1.9, yend=Ostats[3,3]), size=.5, arrow=arrow(length = unit(0.3, "cm")), color=cols[1])+
        geom_point(aes(x=1,    y=Ostats[1,3]),col="#0433FF", bg="#0433FF", pch=21, size=3)+
        geom_point(aes(x=1.95, y=Ostats[2,3]),col="#0433FF",bg="#0433FF", pch=21, size=3)+
        geom_point(aes(x=2   , y=Ostats[3,3]),col=cols[1],bg="#155FF2", pch=21, size=3.5)+
        theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
        ylim(0,33)+
        theme(plot.margin = unit(c(1,1,2,1), "lines"))+
        annotate(geom="text", x=1.1, y=11, label="a",color ='black', size=3)+
        annotate(geom="text", x=2.1, y=19, label="b",color ='black', size=3)+
                annotate(geom="text", x=2.12, y=27.5, label="c",color ='black', size=3)+
        annotation_custom(text_before, xmin=1, xmax=1, ymin=-5, ymax=-5)+
        annotation_custom(text_after, xmin=2, xmax=2, ymin=-5, ymax=-5)+
        annotate("point", x = 2, y =5 , colour = "#0433FF", size=3)+annotate(geom="text", x=2.3, y=5, label=expression("O" %->% "O"),size=3)+
        annotate("point", x = 2, y =2.5 , bg="#155FF2",col=cols[1], pch=21, size=3.5)+
        annotate(geom="text", x=2.3, y=2.5, label=expression("O" %->% "N"),size=3)+
        coord_cartesian(clip="off")

O_plot 

pdf("Output/FigS9_Lipids.pdf",width=6, height=3.8)
grid.arrange(N_plot,O_plot,nrow=1)
dev.off()



