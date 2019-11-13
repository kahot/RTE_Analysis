library(plotrix)
library(ggplot2)
library(ggthemes)
cols<-c("#ED8636","#0433FF")

after<-read.csv("Data/tissue_thickness_A.csv", header=T) 
nn<-subset(after, Treatment=="NN")
no<-subset(after, Treatment=="NO")
on<-subset(after, Treatment=="ON")
oo<-subset(after, Treatment=="OO")

Origin<-c("Nearshore","Nearshore","Offshore","Offshore")
stats2<-data.frame(Site=c("Nearshore","Offshore","Nearshore","Offshore"),
                   Origin=c("Nearshore","Nearshore","Offshore","Offshore"),
                   Type=c("NN","ON","NO","OO"), 
                   Tissue=c(mean(nn$Thickness), mean(no$Thickness),mean(on$Thickness),mean(oo$Thickness)), 
                   SE=c(std.error(nn$Thickness), std.error(no$Thickness), std.error(on$Thickness),std.error(oo$Thickness)))



ggplot(stats2, aes(x=Site,y=Tissue,color=factor(Origin)))+geom_point()+
    geom_errorbar(mapping=aes(ymin=Tissue-SE, ymax=Tissue+SE,color=factor(Origin)), width=.1)+
    scale_color_manual("Origin", breaks=c(1,2), values=cols)+
    geom_point(aes(fill=factor(Origin)),size=5)+
    scale_fill_manual("Origin", breaks=c(1,2), values=cols)+
    labs(x="Transplant site",y="Tissue layer thickness (mm)") +theme_base()+theme(legend.position="none", text = element_text(size=15))+scale_y_continuous(limits = c(2.4, 3.2))+
    geom_segment(data=stats2, mapping=aes(x=1, y=stats2[1,4], xend=1.90, yend=stats2[2,4]), size=.8, arrow=arrow(length = unit(0.3, "cm")),color='#ED8636')+
    geom_segment(data=stats2, mapping=aes(x=2, y=stats2[4,4], xend=1.1, yend=stats2[3,4]), size=.8, arrow=arrow(length = unit(0.3, "cm")), color='#0433FF')
ggsave(filename="~/Dropbox/R/RTE_Analysis/TissueLayer.pdf",width=4, height=3.8, units='in',device='pdf', plot=plot1 )

#no arrows
plot2<-ggplot(stats2, aes(x=Site,y=Tissue,color=factor(Origin)))+geom_point()+
    geom_errorbar(mapping=aes(ymin=Tissue-SE, ymax=Tissue+SE,color=factor(Origin)), width=.17)+
    scale_color_manual("Origin", breaks=c(1,2), values=cols)+
    geom_point(aes(fill=factor(Origin)),size=5)+
    scale_fill_manual("Origin", breaks=c(1,2), values=cols)+
    theme(plot.margin = unit(c(0,0,0,0), "lines"))+
    labs(x="Transplant site",y="Tissue layer thickness (mm)") +theme_bw()+theme(legend.position="none")+scale_y_continuous(limits = c(2.4, 3.2))+
    geom_segment(data=stats2, mapping=aes(x=1, y=stats2[1,4], xend=2, yend=stats2[2,4]), size=.8, color=cols[1])+
    geom_segment(data=stats2, mapping=aes(x=2, y=stats2[4,4], xend=1, yend=stats2[3,4]), size=.8,  color='#0433FF')+
    theme(axis.text.x =element_text(size=12, color="black"), panel.grid.major.x = element_blank(), axis.title = element_text(size=14))

plot2<-plot2+annotate("point", x = 2, y =2.52 , colour = cols[1], size=4)+annotate(geom="text", x=2.3, y=2.52, label="Nearshore",color =cols[1], size=4)+
    annotate("point", x = 2, y =2.45 , colour = cols[2], size=4)+annotate(geom="text", x=2.26, y=2.45, label="Offshore",color =cols[2], size=4)+
    annotate(geom="text", x=1.1, y=3.05, label="a",color ='black', size=4)+
    annotate(geom="text", x=1.05, y=2.65, label="b",color ='black', size=4)+
    annotate(geom="text", x=2.1, y=3.1, label="a",color ='black', size=4)+
    annotate(geom="text", x=2.1, y=2.95, label="a",color ='black', size=4)
plot2
ggsave(filename="Output/TissueLayer1.pdf",width=4, height=3.8, units='in',device='pdf', plot=plot2 )

