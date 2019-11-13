##1. 7.15.16 SOD1

#read data
SOD<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_SOD.txt", header=T)  #C->C normalized to O->N samples
#SOD<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_SOD2.txt", header=T) #C->C normalized to N->O samples

str(SOD)

#check assumption
library(Rcmdr)
leveneTest(Intensity ~ Site, data=SOD)
#Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
#group 1   0.1873 0.6703
leveneTest(Intensity ~ Origin, data=SOD)


#ANOVA

anova_sod2<-aov(Intensity~Site*Origin,data=SOD)
summary(anova_sod2)
#           Df    Sum Sq   Mean Sq F value  Pr(>F)   
#Site         1 7.450e+09 7.450e+09   8.592 0.00979 **
#Origin       1 1.341e+09 1.341e+09   1.547 0.23155   
#Site:Origin  1 1.613e+09 1.613e+09   1.860 0.19149   
#Residuals   16 1.387e+10 8.671e+08 

#SOD2 Results:
#            Df    Sum Sq   Mean Sq F value   Pr(>F)    
#Site         1 1.244e+10 1.244e+10  16.555 0.000893 ***
#Origin       1 3.822e+09 3.822e+09   5.088 0.038447 *  
#Site:Origin  1 2.238e+08 2.238e+08   0.298 0.592710    
#Residuals   16 1.202e+10 7.511e+08

op<-par(mfrow=c(2,2)); 
plot(anova_sod2); 
par(op)

TukeyHSD(x=anova_sod2,'Site', conf.level=0.95)
#             diff      lwr      upr     p adj
#SiteN-CW 38599.94 10683.84 66516.03 0.0097854

TukeyHSD(x=anova_sod2,'Site:Origin', conf.level=0.95)
#$`Site:Origin`
#                      diff        lwr       upr     p adj
#SiteN:C-CW:C     20639.873 -32641.340  73921.09 0.6897617
#CW:N-CW:C        -1583.318 -54864.531  51697.90 0.9997699
#SiteN:N-CW:C     54976.682   1695.469 108257.90 0.0419512
#CW:N-SiteN:C    -22223.191 -75504.404  31058.02 0.6395669
#SiteN:N-SiteN:C  34336.809 -18944.404  87618.02 0.2900912
#SiteN:N-CW:N     56560.000   3278.787 109841.21 0.0355508


sod2<-aov(Intensity~Treatment,data=SOD)
summary(sod2)

#            Df    Sum Sq   Mean Sq F value Pr(>F)  
#Treatment    3 1.040e+10 3.468e+09       4 0.0266 *
#Residuals   16 1.387e+10 8.671e+08                 

#SOD2 Result
#            Df    Sum Sq   Mean Sq F value   Pr(>F)    
#Treatment    3 1.648e+10 5.494e+09   7.314 0.00264 **
#Residuals   16 1.202e+10 7.511e+08




TukeyHSD(x=sod2,'Treatment', conf.level=0.95)
#            diff        lwr       upr     p adj
#CN-CC  20639.873 -32641.340  73921.09 0.6897617
#NC-CC  -1583.318 -54864.531  51697.90 0.9997699
#NN-CC  54976.682   1695.469 108257.90 0.0419512*
#NC-CN -22223.191 -75504.404  31058.02 0.6395669
#NN-CN  34336.809 -18944.404  87618.02 0.2900912
#NN-NC  56560.000   3278.787 109841.21 0.0355508*

#SOD2 Results
#           diff        lwr       upr     p adj
#CN-CC  43179.79  -6411.830  92771.40 0.0996328 
#NC-CC  20956.60 -28635.021  70548.21 0.6302132
#NN-CC  77516.60  27924.979 127108.21 0.0019635 **
#NC-CN -22223.19 -71814.807  27368.43 0.5866509
#NN-CN  34336.81 -15254.807  83928.43 0.2357535
#NN-NC  56560.00   6968.384 106151.62 0.0227514 *

## **SOD and SOD2 results are consistent-> use the data SOD for the paper 

###################################################################
###################################################################
#2.Ferrochelatase (compare across the membranes adjusting to C-> N or N->C )

f<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_Ferrochelatase.txt", header=T) 
str(f)

#check assumption
leveneTest(Intensity ~ Site*Origin, data=f)
#       Df F value Pr(>F)
#group  3  1.1574 0.3566

#2-Way ANOVA
anova_f<-aov(Intensity~Site*Origin, data=f)
summary(anova_f)
#            Df    Sum Sq   Mean Sq F value   Pr(>F)    
#Site         1 135933361 135933361   4.655 0.0465 * 
#Origin       1 273229908 273229908   9.357 0.0075 **
#Site:Origin  1 154897103 154897103   5.305 0.0350 * 
#Residuals   16 467196589  29199787  

TukeyHSD(x=anova_f,'Site:Origin', conf.level=0.95)
#                      diff        lwr      upr     p adj
#SiteN:C-CW:C       351.832  -9425.958 10129.62 0.9995920
#CW:N-CW:C        12958.208   3180.418 22736.00 0.0078300
#SiteN:N-CW:C      2178.208  -7599.582 11956.00 0.9183804
#CW:N-SiteN:C     12606.376   2828.586 22384.17 0.0096520
#SiteN:N-SiteN:C   1826.376  -7951.414 11604.17 0.9493697
#SiteN:N-CW:N    -10780.000 -20557.790 -1002.21 0.0282335





## 1-way ANOVA
a1 <- aov(Intensity~Treatment,data=f)
summary(a1)
#            Df    Sum Sq  Mean Sq F value  Pr(>F)   
#Treatment    3 564060373 1.88e+08   6.439 0.00458 **
#Residuals   16 467196589 2.92e+07

TukeyHSD(x=a1,conf.level=0.95)

#           diff        lwr       upr     p adj
#CN-CC    351.832  -9425.958 10129.62 0.9995920
#NC-CC  12958.208   3180.418 22736.00 0.0078300 *
#NN-CC   2178.208  -7599.582 11956.00 0.9183804
#NC-CN  12606.376   2828.586 22384.17 0.0096520 *
#NN-CN   1826.376  -7951.414 11604.17 0.9493697
#NN-NC -10780.000 -20557.790 -1002.21 0.0282335 *

a2 <- aov(Intensity~Site*Origin,data=f)
TukeyHSD(x=a2,conf.level=0.95)
#$Site
#              diff       lwr       upr     p adj
#SiteN-CW -5214.084 -10337.05 -91.12021 0.0465054

#$Origin
#               diff      lwr      upr     p adj
#N-C 7392.292 2269.328 12515.26 0.0074969



## Separate by transplant site ###
fN<-f[1:10,]
fC<-f[11:20,]

anova_fN<-aov(Intensity~Site,data=fN)
summary(anova_fN)
#            Df   Sum Sq  Mean Sq F value Pr(>F)  
#Site         1 290521000 290521000   6.487 0.0343 *
#Residuals    8 358308000  44788500    

#genotype effects at Site CW.
anova_fC<-aov(Intensity~Site,data=fC)
summary(anova_fC)
#            Df    Sum Sq  Mean Sq F value Pr(>F)
#Site         1    309464   309464   0.023  0.884
#Residuals    8 108888589 13611074   

fT1<-f[1:5,]
fT2<-f[11:15,]
fTT<-rbind(fT1,fT2)
anova_fTT<-aov(Intensity~Origin,data=fTT)
summary(anova_fTT)
#            Df    Sum Sq   Mean Sq F value Pr(>F)  
#Origin       1 397301794 397301794   10.77 0.0112 *
#Residuals    8 295116158  36889520

###################################################################
###################################################################
#PGK

# 1.PGK Reanalyzed on 7.15.16
p<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_PGK3-2.txt", header=T) (adjusted to N->C )

a <- aov(Intensity~Treatment,data=p)
summary(a)
#            Df   Sum Sq Mean Sq F value Pr(>F)
#Treatment    3 161550569 53850190   2.431  0.103
#Residuals   16 354489740 221556093               
TukeyHSD(x=a,conf.level=0.95)
#           diff         lwr       upr     p adj
#CN-CC  7668.436   -848.6831 16185.555 0.0853934.
#NC-CC  1769.813  -6747.3057 10286.932 0.9322729
#NN-CC  3443.813  -5073.3057 11960.932 0.6612178
#NC-CN -5898.623 -14415.7415  2618.496 0.2355759
#NN-CN -4224.623 -12741.7415  4292.496 0.5061275
#NN-NC  1674.000  -6843.1188 10191.119 0.9417772



#3. PGK #3 (adjusted to N->C )
p<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_PGK2.txt", header=T) 

leveneTest(Intensity ~ Treatment, data=p)

a1 <- aov(Intensity~Site*Origin,data=p)
summary(a1)
#            Df   Sum Sq  Mean Sq F value  Pr(>F)   
#Site         1 42177687 42177687   8.806 0.00907 **
#Origin       1  2434263  2434263   0.508 0.48617   
#Site:Origin  1  9280664  9280664   1.938 0.18297   
#Residuals   16 76633391  4789587            

a <- aov(Intensity~Treatment,data=p)
summary(a)
#            Df   Sum Sq  Mean Sq F value Pr(>F)
#Treatment    3 53892615 17964205   3.751 0.0325 *
#Residuals   16 76633391  4789587                 

TukeyHSD(x=a,conf.level=0.95)
#            diff        lwr      upr     p adj
#CN-CC  4266.7993   306.7564 8226.842 0.0325135*
#NC-CC   664.6514 -3295.3915 4624.694 0.9623925
#NN-CC  2206.6514 -1753.3915 6166.694 0.4092381
#NC-CN -3602.1480 -7562.1909  357.895 0.0813320
#NN-CN -2060.1480 -6020.1909 1899.895 0.4667498
#NN-NC  1542.0000 -2418.0429 5502.043 0.6863939



###################################################################
###################################################################
# Hsp60
h60<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_hsp60.txt", header=T)  

a2 <- aov(Intensity~Treatment,data=h60)
summary(a2)
#             Df    Sum Sq  Mean Sq F value Pr(>F)  
#Treatment    3 122611952 40870651   5.028 0.0131 *
#Residuals   15 121919302  8127953

TukeyHSD(x=a2,conf.level=0.95)
#            diff        lwr       upr     p adj
#CN-CC -2426.5411 -7938.5922  3085.510 0.5953074
#NC-CC  4751.3296  -445.4821  9948.141 0.0789190 .
#NN-CC   725.4589 -4471.3527  5922.271 0.9771630
#NC-CN  7177.8706  1665.8195 12689.922 0.0092503 **
#NN-CN  3152.0000 -2360.0511  8664.051 0.3832458
#NN-NC -4025.8706 -9222.6823  1170.941 0.1592030

a1<- aov(Intensity~Origin* Site,data=h60)
summary(a1)
#            Df    Sum Sq  Mean Sq F value  Pr(>F)   
#Origin       1  44096112 44096112   5.648 0.03030 * 
#Site         1  75506599 75506599   9.670 0.00674 **
#Origin:Site  1   3009241  3009241   0.370 0.55198   
#Residuals   16 124928543  7808034

###################################################################
###################################################################
##################################
## 5. Catalse anlysis3 (3.8.16 data, C->N adjusted)###
cat<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_CatalaseC-N.txt", header=T) 

#ANOVA
anova_cat1<-aov(Intensity~Site*Origin,data=cat)
summary(anova_cat1)
#            Df    Sum Sq   Mean Sq F value  Pr(>F)   
#Site         1 5.645e+08 564519532     9.4 0.00700 **
#Origin       1 6.906e+08 690597624    11.5 0.00347 **
#Site:Origin  1  46562075  46562075   0.765 0.39483
#Residuals   17 1.021e+09  60053401

#### pairwise comparison of treatmnet ###
a1 <- aov(Intensity~Treatment,data=cat)
TukeyHSD(x=a1,conf.level=0.95)
#           diff        lwr       upr     p adj
#CN-CC -13677.250 -27797.66   443.1636 0.0593309 .
#NC-CC -14804.051 -28924.46  -683.6377 0.0382540 *
#NN-CC -22378.051 -36498.46 -8257.6377 0.0017320 *
#NC-CN  -1126.801 -15247.21 12993.6123 0.9956311
#NN-CN  -8700.801 -22821.21  5419.6123 0.3259108
#NN-NC  -7574.000 -21694.41  6546.4136 0.4412050




###################################################################
###################################################################
# CaM binding protein
cam<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_CaM55.txt", header=T) 

## 2-way anova
a1 <- aov(Intensity~Origin* Site,data=cam)
summary(a1)
#            Df    Sum Sq   Mean Sq F value  Pr(>F)   
#Origin       1 1.419e+09 1.419e+09   1.302 0.26970   
#Site         1 9.556e+09 9.556e+09   8.767 0.00875 **
#Origin:Site  1 1.218e+09 1.218e+09   1.126 0.30444  
#Residuals   17 1.853e+10 1.090e+09

#### pairwise comparison of treatmnet ###
a2 <- aov(Intensity~Treatment,data=cam)
TukeyHSD(x=a2,conf.level=0.95)
#$Treatment
#           diff        lwr       upr     p adj
#CN-CC  1238.354 -58279.6621  60756.37 0.9999209
#NC-CC 28108.862 -31409.1539  87626.88 0.5457189 
#NN-CC 60561.815   1043.7989 120079.83 0.0454008 *
#NC-CN 26870.508 -32647.5079  86388.52 0.5809572
#NN-CN 59323.461   -194.5551 118841.48 0.0509040 .
#NN-NC 32452.953 -27065.0633  91970.97 0.4274471


###################################################################
###################################################################
# CYP1A 

cyp<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_CYP2.txt", header=T) 
str(cyp)

anova_cyp1<-aov(Intensity~Site*Origin,data=cyp)
summary(anova_cyp1)
#            Df    Sum Sq   Mean Sq F value  Pr(>F)   
#Site         1    291829   291829   0.023  0.881
#Origin       1   2022945  2022945   0.161  0.694
#Site:Origin  1   6839703  6839703   0.543  0.472
#Residuals   16 201538030 12596127               
---
  op<-par(mfrow=c(2,2))
plot(anova_cyp1)

anova_cyp2<-aov(Intensity~Treatment,data=cyp)
summary(anova_cyp2)
#            Df    Sum Sq  Mean Sq F value Pr(>F)
#Treatment    3   9154476  3051492   0.242  0.866
#Residuals   16 201538030 12596127 

TukeyHSD(x=anova_cyp2,conf.level=0.95)
#           diff       lwr      upr     p adj
#CN-CC 1411.1799 -5010.806 7833.166 0.9213035
#NC-CC 1805.6630 -4616.323 8227.649 0.8513829
#NN-CC  877.6630 -5544.323 7299.649 0.9790026
#NC-CN  394.4831 -6027.503 6816.469 0.9979884
#NN-CN -533.5169 -6955.503 5888.469 0.9950799
#NN-NC -928.0000 -7349.986 5493.986 0.9753725


#################################################
###8. Actin #2 (adjusted to N->C )
a<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_actin.txt", header=T) 

## Separate by Genotype(Origin) ###
#aN<-a[1:10,]
#aC<-a[11:20,]

#ANOVA
anova_a1<-aov(Intensity~Site*Origin,data=a)
summary(anova_a1)
#           Df    Sum Sq   Mean Sq F value  Pr(>F)   
#Site         1    246968   246968   0.013  0.909
#Origin       1  46601239 46601239   2.538  0.131
#Site:Origin  1  31758213 31758213   1.730  0.207
#Residuals   16 293730643 18358165               
op<-par(mfrow=c(2,2)); 
plot(anova_p1); 
par(op)

anova_a2<-aov(Intensity~Treatment,data=a)
summary(anova_a2)
#            Df    Sum Sq  Mean Sq F value Pr(>F)
#Treatment    3  78606420 26202140   1.427  0.272
#Residuals   16 293730643 18358165               
TukeyHSD(x=anova_a2,conf.level=0.95)
#            diff        lwr       upr     p adj
#CN-CC  2742.4931  -5010.430 10495.417 0.7448627
#NC-CC  -532.6615  -8285.585  7220.262 0.9971951
#NN-CC -2830.6615 -10583.585  4922.262 0.7265174
#NC-CN -3275.1546 -11028.078  4477.769 0.6304556
#NN-CN -5573.1546 -13326.078  2179.769 0.2092273
#NN-NC -2298.0000 -10050.924  5454.924 0.8308623

> 
