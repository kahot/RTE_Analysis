
#7.15.16
#####################################################################
### Nested Anova test on after 30day samples ##
after<-read.csv("/Users/kahotisthammer/Documents/R/tissue_thickness_A.csv", header=T) 
str(after)
library(Rcmdr)
leveneTest(Thickness ~ Treatment, data=after)
#        Df F value Pr(>F)
#group   3  0.8642 0.4606
#     196   


##### comparing 4 treatments  #####
mod<-aov(Thickness~Treatment/Indiv,data=after)
summary(mod)
#                 Df Sum Sq Mean Sq F value   Pr(>F)    
#Treatment         3   6.87   2.289   30.81  4e-16 ***
#Treatment:Indiv  16  71.91   4.494   60.49 <2e-16 ***
#Residuals       180  13.37   0.074                   
 

TukeyHSD(x=mod,'Treatment', conf.level=0.99)

#$Treatment
#        diff         lwr         upr     p adj
#NO-NN  0.0392 -0.132922  0.21132199 0.8893995
#ON-NN -0.4256 -0.597722 -0.25347801 0.0000000
#OO-NN -0.0518 -0.223922  0.12032199 0.7777140
#ON-NO -0.4648 -0.636922 -0.29267801 0.0000000
#OO-NO -0.0910 -0.263122  0.08112199 0.3429679
#OO-ON  0.3738  0.201678  0.54592199 0.0000000






####  nested 2-way anova   
mod1<-aov(Thickness~Site*Origin/Indiv,data=after)
summary(mod1)
#                   Df Sum Sq Mean Sq F value   Pr(>F)    
#Site                1   2.13   2.132   28.70 2.56e-07 ***
#Origin              1   3.34   3.336   44.90 2.57e-10 ***
#Site:Origin         1   1.40   1.399   18.84 2.37e-05 ***
#Site:Origin:Indiv  16  71.91   4.494   60.49  < 2e-16 ***
#Residuals         180  13.37   0.074

TukeyHSD(x=mod1,'Site:Origin', conf.level=0.95)
#                           diff         lwr        upr     p adj
#SiteN:SiteC-SiteC:SiteC -0.3738 -0.51516368 -0.2324363 0.0000000
#SiteC:SiteN-SiteC:SiteC  0.0910 -0.05036368  0.2323637 0.3429679
#SiteN:SiteN-SiteC:SiteC  0.0518 -0.08956368  0.1931637 0.7777140
#SiteC:SiteN-SiteN:SiteC  0.4648  0.32343632  0.6061637 0.0000000
#SiteN:SiteN-SiteN:SiteC  0.4256  0.28423632  0.5669637 0.0000000
#SiteN:SiteN-SiteC:SiteN -0.0392 -0.18056368  0.1021637 0.8893995

#Rewriting to Genotyp:Site style, which is the same as 1-way anaova results
#ON-OO -0.3738 -0.51516368 -0.2324363 0.0000000
#NO-OO  0.0910 -0.05036368  0.2323637 0.3429679
#NN-OO  0.0518 -0.08956368  0.1931637 0.7777140
#NO-ON  0.4648  0.32343632  0.6061637 0.0000000
#NN-ON  0.4256  0.28423632  0.5669637 0.0000000
#NN-NO -0.0392 -0.18056368  0.1021637 0.8893995






mod2<-aov(Thickness~Site+Origin+ Origin/Indiv + Site:Origin/Indiv,data=after)
summary(mod2)
#                   Df Sum Sq Mean Sq F value   Pr(>F)    
#Site                1   2.13   2.132   28.70 2.56e-07 ***
#Origin              1   3.34   3.336   44.90 2.57e-10 ***
#Origin:Indiv        8  64.56   8.069  108.62  < 2e-16 ***
#Site:Origin         1   1.40   1.399   18.84 2.37e-05 ***
#Site:Origin:Indiv   8   7.35   0.919   12.37 4.66e-14 ***
#Residuals         180  13.37   0.074 4                     

TukeyHSD(x=mod2, conf.level=0.99,"Origin")

#              diff       lwr       upr p adj
#SiteN-SiteC 0.2583 0.1579466 0.3586534     0

mod3<-aov(Thickness~Origin/Indiv,data=after)
summary(mod3)



boxplot(Thickness ~ Treatment, data=after)
mean(after$Thickness)  #2.90765

NN<-after[1:50,]
mean(NN$Thickness)  #3.0172
mean(NN$Thickness, trim = 0.5) # 2.865

NO<-after[101:151,]
mean(NO$Thickness)  #3.06
mean(NO$Thickness, trim = 0.5) #3.04

OO<-after[151:200,]
mean(OO$Thickness)  #2.9654
mean(OO$Thickness, trim = 0.5) #2.975

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################


Nafter<-which(after$Origin=='SiteN')
N.after<-after[Nafter,]

Oafter<-which(after$Origin=='SiteC')
O.after<-after[Oafter,]

## Nested ANOVA
mod<-aov(Thickness~Site/Indiv,trim=0, data=N.after)
summary(mod)
#same as above
mod<-aov(Thickness~Site + Indiv %in% Site,data=N.after)
summary(mod)
#      Df Sum Sq Mean Sq F value  Pr(>F)    
#Site         1   0.038   0.038   0.471  0.494    
#Site:Indiv   8 30.633   3.829  46.942 <2e-16 ***
#Residuals   90  7.342   0.082

#ignore the nesting: 
mod2<-aov(Thickness~Site,data=N.after)
summary(mod2)
#           Df Sum Sq Mean Sq F value Pr(>F)
#Site         1   0.04  0.0384   0.099  0.754
#Residuals   98  37.97  0.3875 


## Offshore Samples#3##
mod3<-aov(Thickness~Site/Indiv,data=O.after)
summary(mod3)
#            Df Sum Sq Mean Sq F value  Pr(>F)    
#Site         1   3.49   3.493   52.13 1.6e-10 ***
#Site:Indiv   8  41.27   5.159   76.99 < 2e-16 ***
#Residuals   90   6.03   0.067                    
---
#ignore the nesting: 
mod4<-aov(Thickness~Site,data=O.after)
summary(mod4)
#            Df Sum Sq Mean Sq F value Pr(>F)   
#Site         1   3.49   3.493   7.237 0.0084 **
#Residuals   98  47.30   0.48



#####################################################################

## Nested ANOVA on before samples ####
before<-read.csv("/Users/kahotisthammer/Documents/R/tissue_thickness_B.csv", header=T) 
str(before)

mod5<-aov(Thickness~Origin/Indiv,data=before)
summary(mod5)
#            Df Sum Sq Mean Sq F value  Pr(>F)    
#Origin        1  4.541   4.541   27.19 1.17e-06 ***
#Origin:Indiv  8 21.748   2.718   16.28 1.22e-14 ***
#Residuals    90 15.030   0.167                    

#ignore the nesting: 
mod6<-aov(Thickness~Origin,data=before)
summary(mod6)
#            Df Sum Sq Mean Sq F value  Pr(>F)    
#Origin       1   4.54   4.541    12.1 s ***
#Residuals   98  36.78   0.375                     
#####################################################################

## Before After Comparison ###
# 1. Nearshore
N<-read.csv("/Users/kahotisthammer/Dropbox/ResearchResults/RTE/N.csv", header=T) 
str(N)

tNN<-N[1:100,]
tNO1<-N[1:50,]
tNO2<-N[101:150,]
tNO<-rbind(tNO1,tNO2)

#Nearshore->Nearshore
mod7<-aov(Thickness~Time/Indiv,data=tNN)
summary(mod7)
#            Df Sum Sq Mean Sq F value  Pr(>F)    
#Time         1  3.382   3.382   31.53 2.16e-07 ***
#Time:Indiv   8 28.545   3.568   33.27  < 2e-16 ***
#Residuals   90  9.652   0.107                    

#ignore the nesting: 
mod8<-aov(Thickness~Time,data=tNN)
summary(mod8)
#            Df Sum Sq Mean Sq F value  Pr(>F)    
#Time         1   3.38   3.382   8.677 0.00403 **
#Residuals   98  38.20   0.390                   

#Nearshore->offshore
mod9<-aov(Thickness~Time/Indiv,data=tNO)
summary(mod9)
#            Df Sum Sq Mean Sq F value  Pr(>F)    
#Time         1   6.76   6.760   50.91 2.38e-10 ***
#Time:Indiv   8  29.92   3.740   28.17  < 2e-16 ***
#Residuals   90  11.95   0.133                     
mod<-aov(Thickness~Time,data=tNO)
summary(mod)
#            Df Sum Sq Mean Sq F value   Pr(>F)    
#Time         1   6.76   6.760   15.82 0.000133 ***
#Residuals   98  41.87   0.4

## comparing the 3 treatments

mod<-aov(Thickness~Treatment/Indiv, data=N.after)
summary(mod)

####################
#2. Offshore
O<-read.csv("/Users/kahotisthammer/Dropbox/ResearchResults/RTE/O.csv", header=T) 
str(O)

tON<-O[1:100,]
tOO1<-O[1:50,]
tOO2<-O[101:150,]
tOO<-rbind(tOO1,tOO2)
mo1<-aov(Thickness~Time/Indiv,data=tOO)
summary(mo1)
#            Df Sum Sq Mean Sq F value   Pr(>F)    
#Time         1   0.00  0.0002   0.002  0.967    
#Time:Indiv   8  24.80  3.1006  26.596 <2e-16 ***
#Residuals   90  10.49  0.1166  

mo2<-aov(Thickness~Time,data=tOO)
summary(mo2)
#            Df Sum Sq Mean Sq F value   Pr(>F)    
#Time         1    0.0  0.0002   0.001  0.981
#Residuals   98   35.3  0.3602  


mo3<-aov(Thickness~Time/Indiv,data=tON)
summary(mo3)
#            Df Sum Sq Mean Sq F value   Pr(>F)    
#Time         1   3.44   3.441   28.18 7.91e-07 ***
#Time:Indiv   8  34.71   4.339   35.52  < 2e-16 ***
#Residuals   90  10.99   0.122  

mo4<-aov(Thickness~Time,data=tON)
summary(mo4)
#            Df Sum Sq Mean Sq F value   Pr(>F)    
#Time         1   3.44   3.441   7.379 0.0078 **
#Residuals   98  45.70   0.466      


#############################################################
## NN vs. ON Pairwise t test showed variable results based on the correction menthods

NN<-N[51:100,]
ON<-O[51:100,]
SiteN<-rbind(NN,ON)

mo5<-aov(Thickness~Origin/Indiv,data=SiteN)
summary(mo5)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
#Origin        1   2.44   2.443   39.16 1.29e-08 ***
#Origin:Indiv  8  41.51   5.188   83.17  < 2e-16 ***
#Residuals    90   5.61   0.062

## NN-OO
OO<-O[101:150,]
NNOO<-rbind(NN,OO)
mo<-aov(Thickness~Origin/Indiv,data=NNOO)
summary(mo)
#             Df Sum Sq Mean Sq F value Pr(>F)    
#Origin        1  0.094   0.094   1.648  0.203    
#Origin:Indiv  8 31.602   3.950  69.515 <2e-16 ***
#Residuals    90  5.114   0.057 

NO<-N[101:150,]
OONO<-rbind(OO,NO)
mo<-aov(Thickness~Origin/Indiv,data=OONO)
summary(mo)
#             Df Sum Sq Mean Sq F value Pr(>F)    
#Origin        1   0.21   0.207   2.513  0.116    
#Origin:Indiv  8  32.98   4.122  50.048 <2e-16 ***
#Residuals    90   7.41   0.082  

NOON<-rbind(NO,ON)
mo<-aov(Thickness~Origin/Indiv,data=NOON)
summary(mo)
#             Df Sum Sq Mean Sq F value Pr(>F)    
#Origin        1   5.40   5.401   61.43 8.8e-12 ***
#Origin:Indiv  8  42.88   5.360   60.97 < 2e-16 ***
#Residuals    90   7.91   0.088 




#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#BEfore after cahnge
mm<-read.csv("/Users/kahotisthammer/Dropbox/ResearchResults/RTE/TissueLayer_change.csv", header=T) 
str(mm)
mo<-aov(mm_change~Treatment,data=mm)
summary(mo)


leveneTest(mm_change ~ Treatment, data=mm)
#        Df F value Pr(>F)
#group   3  0.3075 0.8196

pairwise.t.test(mm$mm_change, mm$Treatment , p.adj = "none",pool.sd= FALSE)
#   NN      NO      ON     
#NO 0.00022 -       -      
#ON 1.8e-05 0.03536 -      
#OO 2.3e-05 0.18441 0.1665

pairwise.t.test(mm$mm_change, mm$Treatment, p.adj = "fdr")
#   NN      NO      ON     
#NO 1.9e-06 -     -    
#ON 9.4e-08 0.022 -    
#OO 2.4e-07 0.158 0.268


#################################################################################################################################
############################################################
# Pairwise -not appropriate #
tissue<-read.csv("/Users/kahotisthammer/Documents/R/tissue_thickness.csv", header=T) 
str(tissue)


pairwise.t.test(tissue$Thickness, tissue$Treatment , p.adj = "none",pool.sd= FALSE)
#   N       NN      NO      O       ON     
#NN 0.00403 -       -       -       -      
#NO 0.00013 0.23770 -       -       -      
#O  0.00076 0.62676 0.45780 -       -      
#ON 0.69654 0.02653 0.00171 0.00789 -      
#OO 0.00089 0.61825 0.48019 0.98144 0.00846

pairwise.t.test(tissue$Thickness, tissue$Treatment, p.adj = "fdr")
#   N       NN      NO      O       ON     
#NN 0.0107 -      -      -      -     
#NO 0.0012 0.4046 -      -      -     
#O  0.0044 0.7195 0.6608 -      -     
#ON 0.7195 0.0316 0.0031 0.0107 -     
#OO 0.0044 0.7195 0.6608 0.9828 0.0107


pairwise.t.test(tissue$Thickness, tissue$Treatment, p.adj = "holm")
#   N       NN      NO      O       ON    
#NN 0.0478 -      -      -      -     
#NO 0.0012 1.0000 -      -      -     
#O  0.0142 1.0000 1.0000 -      -     
#ON 1.0000 0.1347 0.0057 0.0478 -     
#OO 0.0142 1.0000 1.0000 1.0000 0.0478

pairwise.t.test(tissue$Thickness, tissue$Treatment, p.adj = "bonf")
#   N       NN      NO      O       ON    
#NN 0.0750 -      -      -      -     
#NO 0.0012 1.0000 -      -      -     
#O  0.0176 1.0000 1.0000 -      -     
#ON 1.0000 0.2526 0.0062 0.0696 -     
#OO 0.0163 1.0000 1.0000 1.0000 0.0651




