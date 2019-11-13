library(ggplot2)
library(reshape2)
#Growth RATE
cols<-c("#ED8636","#0433FF")

gr<-read.table('Data/growthrate.txt',header=T)

####### Levene's Test  #####
#check assumption
library(Rcmdr)
leveneTest(GrowthRate ~ Sample, data=gr)
#      Df F value Pr(>F)
#group  1  0.4826 0.5069

#"Shapiro-Wilk normality test".
shapiro.test(gr$GrowthRate) 
# = 0.88411, p-value = 0.1454

## ANOVA
ano1<-aov(GrowthRate~Sample,data=gr) 
summary(ano1)     
#            Df   Sum Sq   Mean Sq F value Pr(>F)  
#Sample       1 0.002242 0.0022420   8.767 0.0181 *
#Residuals    8 0.002046 0.000255
print(model.tables(ano1,"means"),digits=3)    
##Grand mean   0.0406796 

##Sample
#N      O 
#0.0557 0.0257 

op<-par(mfrow=c(2,2)); 
plot(ano1); 


ano2<-aov(asin(GrowthRate)~Sample,data=gr) 
summary(ano2)  
#            Df   Sum Sq   Mean Sq F value Pr(>F)  
#Sample       1 0.002247 0.0022474   8.752 0.0182 *
#Residuals    8 0.002054 0.0002568                 
leveneTest(asin(GrowthRate) ~ Sample, data=gr)
#Levene's Test for Homogeneity of Variance (center = median)
#      Df F value Pr(>F)
#group  1  0.4351  0.528
#8 

plot(ano2); 

ano3<-aov(log(GrowthRate)~Sample,data=gr) 
summary(ano3)  
#            Df Sum Sq Mean Sq F value  Pr(>F)   
#Sample       1 1.4904  1.4904   13.09 0.00681 **
#Residuals    8 0.9111  0.1139                   

plot(ano3); 
leveneTest(log(GrowthRate) ~ Sample, data=gr)
#      Df F value Pr(>F)
#group  1  0.0952 0.7656
#8 
plot(ano1)


logitTransform <- function(p) { log(p/(1-p)) }
ano4<-aov(logitTransform(GrowthRate)~Sample,data=gr) 
summary(ano4)  
#            Df Sum Sq Mean Sq F value Pr(>F)  
#Sample       1  1.6140  1.6140   12.96 0.00698 **
#Residuals    8 0.9961  0.1245    
plot(ano4)

leveneTest(logitTransform(GrowthRate) ~ Sample, data=gr)
#      Df F value Pr(>F)
#group  1  0.0647 0.8057
#8 

############################################################

library("Hmisc")
gr<-read.csv('Data/growth_rate2.csv')

pdf("Output/GrowthRate_Plot.pdf", ,width=6.5, height=5)
plot(gr$Nearshore,type="o",pch=19,col=cols[1],yaxt="n",xaxt="n",xlab='',ylab="Growth rate (%)", ylim=c(0,0.085))
axis(1,at=1:12,labels=gr$X, las=2,cex=0.7 )
lines(gr$Offshore,type="o",pch=19,col=cols[2])
ylabel=c('0','','2.0','','4.0','','6.0','','8.0')
axis(2,at=c(0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08),labels=ylabel,tick=1:7)
with(data=gr, expr=errbar(x=1:12,Nearshore, Nearshore+SE_N,Nearshore-SE_N,add=T,col=cols[1],pch=1,errbar.col=cols[1]))
with(data=gr, expr=errbar(x=1:12,Offshore, Offshore+SE_O,Offshore+-SE_O,add=T,pch=1,col="#0433EE",errbar.col='#0433EE'))
legend('topleft',c('Nearshore corals','Offshore corals'), bty='n', col=cols, pch=c(19,19), cex=0.9)
mtext(expression(italic(P) == 0.0068),side=3, adj=1,cex=0.9)

dev.off()
mtext(expression(paste0('ANOVA ', italic('P'), '= 0.0068')),side=3, adj=1,cex=0.9)

##
gr2<-gr[,1:2]
gr3<-gr[,c(1,3)]
gr2$Site<-"Nearshore"
colnames(gr2)<-c("Week","Rate","Site")
gr3$Site<-"Offshore"
colnames(gr3)<-c("Week","Rate","Site")
gr2<-rbind(gr2,gr3)

gr4<-gr[,c(1,4)]
gr5<-gr[,c(1,5)]
colnames(gr4)<-c("Week","SE")
colnames(gr5)<-c("Week","SE")
gr4<-rbind(gr4,gr5)
gr2$SE<-gr4$SE
gr2$Week<-factor(gr2$Week, c("Week 0","Week 1","Week 2","Week 3","Week 4","Week 5","Week 6","Week 7","Week 8","Week 9","Week 10","Week 11"))
gr2$Rate<-100*gr2$Rate
gr2$SE<-100*gr2$SE


scaleFUN <- function(x) sprintf("%.1f", x)
ggplot(data=gr2, aes(x=Week, y=Rate, color=Site))+
    geom_point()+
    geom_errorbar(aes(ymin=Rate-SE, ymax=Rate+SE), width=.2, size=.2)+ 
    geom_path(aes(x=Week, y=Rate, group=Site),linetype = 1, size=0.2)+
    scale_color_manual(values=cols)+
    theme_bw()+
    theme(axis.text.x = element_text(size=10, angle=90, color="black"))+
    xlab("")+
    ylab("Growth rate (%)")+
    scale_y_continuous(labels=scaleFUN)+
    theme(legend.title = element_blank())+
    theme(axis.text.y = element_text(size=10,color="black"), panel.grid.major.x =element_blank())+
    annotate("text", x=10.5, y=8.4,label= expression(~italic(P)~" = 0.0068"), size=3.5)
ggsave("./Output/GrowthRate.plot.pdf", height = 5, width = 7.4)

    




