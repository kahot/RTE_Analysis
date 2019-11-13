rm(list=ls(all=TRUE))
#######################################################################

## 1. Catalse anlysis1 (3.8.16 cata)###
cat<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_Catalase.txt", header=T) 
str(cat)

#check assumption
library(Rcmdr)
leveneTest(Intensity ~ Origin, data=cat)
#Levene's Test for Homogeneity of Variance (center = median)
#      Df    0.5468 0.4692

leveneTest(Intensity ~ Origin*Site, data=cat)
#      Df F value Pr(>F)
#group  3  0.6863 0.5734

leveneTest(Intensity ~ Site, data=cat)
#      Df F value Pr(>F)
#group  1  0.6683 0.4243


## Separate by Genotype ###
catN<-cat[1:10,]
catC<-cat[11:20,]

#ANOVA
anova_cat1<-aov(Intensity~Site+Origin,data=cat)
summary(anova_cat1)
#            Df    Sum Sq   Mean Sq F value  Pr(>F)   
#Site         1 381510293 381510293   8.537 0.00951 **
#Origin       1 645206441 645206441  14.438 0.00143 **
#Residuals   17 759677638  44686920


op<-par(mfrow=c(2,2))
plot(anova_cat1)


anova_cat2<-aov(Intensity~Site*Origin,data=cat)
summary(anova_cat2)
#           Df    Sum Sq   Mean Sq F value  Pr(>F)   
#Site         1 381510293 381510293   8.107 0.01165 * 
#Origin       1 645206441 645206441  13.711 0.00193 **
#Site:Origin  1   6740824   6740824   0.143 0.71005   
#Residuals   16 752936815  4705855                


op<-par(mfrow=c(2,2)); 
plot(anova_cat2); 
par(op)

### any transplant effects on GenotypeN?
anova_catN<-aov(Intensity~Site,data=catN)
summary(anova_catN)
#            Df   Sum Sq  Mean Sq F value Pr(>F)  
#Origin       1 143413690 143413690   3.228   0.11
#Residuals    8 355410400  44426300  

#transplant effects on Genotype CW
anova_catC<-aov(Intensity~Site,data=catC)
summary(anova_catC)
#            Df   Sum Sq  Mean Sq F value Pr(>F)  
#Site         1 213911428 213911428   3.004  0.121
#Residuals    8 569589572  71198697


#### pairwise comparison of treatmnet ###
a1 <- aov(Intensity~Treatment,data=cat)
TukeyHSD(x=a1,conf.level=0.95)
#           diff        lwr       upr     p adj
#CN-CC  -9250.112 -23008.32  4508.092 0.2574952
#NC-CC -10376.913 -24135.12  3381.290 0.1774531
#NN-CC -17950.913 -31709.12 -4192.710 0.0088229 **
#NC-CN  -1126.801 -14885.00 12631.402 0.9952827
#NN-CN  -8700.801 -22459.00  5057.402 0.3050078
#NN-NC  -7574.000 -21332.20  6184.204 0.4194165


## pairwise t test with correction 
#Correction Methods:
# Bonferroni correction ("bonferroni") 
# Holm (1979) ("holm"), 
# Hochberg (1988) ("hochberg"),
# Hommel (1988) ("hommel"),
# Benjamini & Hochberg (1995) ("BH" or its alias "fdr"), and 
# Benjamini & Yekutieli (2001) ("BY"),

pairwise.t.test(cat$Intensity,cat$Treatment, p.adjust.method="fdr")
#   CC     CN     NC    
#0.134 -     -    
#NC 0.134 0.818 -    
#NN 0.011 0.134 0.162

pairwise.t.test(cat$Intensity,cat$Treatment, p.adjust.method="BY")

pairwise.t.test(cat$Intensity,cat$Treatment, p.adjust.method="none")
#   CC      CN      NC     
#CN 0.0724 -      -     
#NC 0.0465 0.8177 -     
#NN 0.0018 0.0892 0.1348

##################################
## 1.2. Catalse anlysis2 (3.21.16 data)###
cat<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_Catalase1.txt", header=T) 

#ANOVA
anova_cat1<-aov(Intensity~Site+Origin,data=cat)
summary(anova_cat1)
#            Df    Sum Sq   Mean Sq F value  Pr(>F)   
#Site         1 353813428 353813428   6.478 0.0209 *
#Origin       1 454948996 454948996   8.330 0.0103 *
#Residuals   17 928511661  54618333

anova_cat2<-aov(Intensity~Site*Origin,data=cat)
summary(anova_cat2)
#           Df    Sum Sq   Mean Sq F value  Pr(>F)   
#Site         1 353813428 353813428   6.120 0.0250 *
#Origin       1 454948996 454948996   7.869 0.0127 *
#Site:Origin  1   3511689   3511689   0.061 0.8085  
#Residuals   16 924999972  57812498                 

#### pairwise comparison of treatmnet ###
a1 <- aov(Intensity~Treatment,data=cat)
TukeyHSD(x=a1,conf.level=0.95)
#           diff        lwr       upr     p adj
#CN-CC  -9250.112 -23008.32  4508.092 0.2574952
#NC-CC -10376.913 -24135.12  3381.290 0.1774531
#NN-CC -17950.913 -31709.12 -4192.710 0.0088229 *
#NC-CN  -1126.801 -14885.00 12631.402 0.9952827
#NN-CN  -8700.801 -22459.00  5057.402 0.3050078
#NN-NC  -7574.000 -21332.20  6184.204 0.4194165


#pairwise t test with correction
pairwise.t.test(cat$Intensity,cat$Treatment, p.adjust.method="fdr")
#   CC     CN     NC    
#CN 0.0549 -      -     
#NC 0.0323 0.5537 -     
#NN 0.0017 0.0549 0.1200

#######################################################################
##################################
## 1.3. Catalse anlysis3 (3.8.16 data, C->N adjusted)###
cat<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_CatalaseC-N.txt", header=T) 

#ANOVA
anova_cat1<-aov(Intensity~Site+Origin,data=cat)
summary(anova_cat1)
#            Df    Sum Sq   Mean Sq F value  Pr(>F)   
#Site         1 5.645e+08 564519532     9.4 0.00700 **
#Origin       1 6.906e+08 690597624    11.5 0.00347 **
#Residuals   17 1.021e+09  60053401

anova_cat2<-aov(Intensity~Site*Origin,data=cat)
summary(anova_cat2)
#           Df    Sum Sq   Mean Sq F value  Pr(>F)   
#Site         1 564519532 564519532   9.270 0.00772 **
#Origin       1 690597624 690597624  11.340 0.00392 **
#Site:Origin  1  46562075  46562075   0.765 0.39483   
#Residuals   16 974345747  6089660                 

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


#pairwise t test with correction
pairwise.t.test(cat$Intensity,cat$Treatment, p.adjust.method="fdr")
#   CC     CN     NC    
#CN 0.027 -     -    
#NC 0.025 0.822 -    
#NN 0.002 0.146 0.173




#######################################################################