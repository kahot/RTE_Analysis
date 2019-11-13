rm(list=ls(all=TRUE))

#2. PGK  No.1 (adjusted to C->N)
p<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_PGK.txt", header=T) 
str(p)

p<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_PGK3.txt", header=T) 

#check assumption
library(Rcmdr)
leveneTest(Intensity ~ Origin, data=p)
leveneTest(Intensity ~ Origin*Site, data=cyp)
leveneTest(Intensity ~ Site, data=cyp)

###
a <- aov(Intensity~Treatment,data=p)
TukeyHSD(x=a,conf.level=0.95)

#            diff        lwr       upr     p adj
#CN-CC  3238.5216  -955.4385 7432.4817 0.1628917
#NC-CC  -363.6264 -4557.5865 3830.3337 0.9944196
#NN-CC  1178.3736 -3015.5865 5372.3337 0.8516443
#NC-CN -3602.1480 -7796.1080  591.8121 0.1058711
#NN-CN -2060.1480 -6254.1080 2133.8121 0.5140596
#NN-NC  1542.0000 -2651.9601 5735.9601 0.7223394



## Separate by Genotype(Origin) ###
pN<-p[1:10,]
pC<-p[11:20,]

leveneTest(Intensity ~ Site, data=pN)
leveneTest(Intensity ~ Site, data=cypC)

#ANOVA
anova_p1<-aov(Intensity~Site*Origin,data=p)
summary(anova_p1)
#           Df    Sum Sq   Mean Sq F value  Pr(>F)   
#Site         1  3597732  3597732   0.670 0.4252  
#Origin       1  7343352  7343352   1.367 0.2595  
#Site:Origin  1 28566733 28566733   5.318 0.0348 *
#Residuals   16 85954149  5372134

## PGK3 ##
#Df   Sum Sq  Mean Sq F value Pr(>F)  
#Site         1 11910116 11910116   4.183 0.0576 .
#Origin       1  5836123  5836123   2.050 0.1715  
#Site:Origin  1       10       10   0.000 0.9986  
#Residuals   16 45558132  2847383 

op<-par(mfrow=c(2,2)); 
plot(anova_p1); 
par(op)

### Any Site effects on N genotype?
anova_pN<-aov(Intensity~Site,data=pN)
summary(anova_pN)
#            Df   Sum Sq  Mean Sq F value Pr(>F)  
#Site         1  5944410 5944410   2.058  0.189
#Residuals    8 23106080 2888260 

#Any Site effects on C genotype?
anova_pC<-aov(Intensity~Site,data=pC)
summary(anova_pC)
#            Df   Sum Sq  Mean Sq F value Pr(>F)  
#Site         1 26220055 26220055   3.338  0.105
#Residuals    8 62848069  7856009  

## PGK3 ##
#Site         1  5965715 5965715   2.126  0.183
#Residuals    8 22452052 2806506


pT1<-p[1:5,]
pT2<-p[11:15,]
pTT<-rbind(pT1,pT2)
anova_pTT<-aov(Intensity~Origin,data=pTT)
summary(anova_pTT)
#            Df   Sum Sq Mean Sq F value Pr(>F)
#Origin       1 32438675 32438675   4.291 0.0721 .
#Residuals    8 60477138  7559642                 
--- 

#PGK3##
#Origin       1   535920  535920   0.134  0.724
#Residuals    8 31968069 3996009 
  
  
pairwise.t.test(p$Intensity,p$Treatment, p.adjust.method="none")
#   CC    CN    NC   
#CN 0.042 -     -    
#NC 0.807 0.026 -    
#NN 0.433 0.179 0.308

pairwise.t.test(p$Intensity,p$Treatment, p.adjust.method="fdr")
#   CC   CN   NC  
#CN 0.13 -    -   
#NC 0.81 0.13 -   
#NN 0.52 0.36 0.46


pairwise.t.test(p$Intensity,p$Treatment, p.adjust.method="bonferroni")


#######################################################################

#2. PGK #2 (adjusted to N->C )
p<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_PGK2.txt", header=T) 

###
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
#NC-CN -3602.1480 -7562.1909  357.895 0.0813320.
#NN-CN -2060.1480 -6020.1909 1899.895 0.4667498
#NN-NC  1542.0000 -2418.0429 5502.043 0.6863939



## Separate by Genotype(Origin) ###
pN<-p[1:10,]
pC<-p[11:20,]

#ANOVA
anova_p1<-aov(Intensity~Site*Origin,data=p)
summary(anova_p1)
#           Df    Sum Sq   Mean Sq F value  Pr(>F)   
#Site         1 42177687 42177687   8.806 0.00907 **
#Origin       1  2434263  2434263   0.508 0.48617   
#Site:Origin  1  9280664  9280664   1.938 0.18297   
#Residuals   16 76633391  4789587 
op<-par(mfrow=c(2,2)); 
plot(anova_p1); 
par(op)

### Any Site effects on N genotype?
anova_pN<-aov(Intensity~Site,data=pN)
summary(anova_pN)
#            Df   Sum Sq  Mean Sq F value Pr(>F)  
#Site         1  5944410 5944410   2.058  0.189
#Residuals    8 23106080 2888260 

#Any Site effects on C genotype?
anova_pC<-aov(Intensity~Site,data=pC)
summary(anova_pC)
#            Df   Sum Sq  Mean Sq F value Pr(>F)  
#Site         1 45513942 45513942   6.802 0.0312 *
#Residuals    8 53527311  6690914  

pT1<-p[1:5,]
pT2<-p[11:15,]
pTT<-rbind(pT1,pT2)
anova_pTT<-aov(Intensity~Origin,data=pTT)
summary(anova_pTT)
#            Df   Sum Sq Mean Sq F value Pr(>F)
#Origin       1 32438675 32438675   4.291 0.0721 .
#Residuals    8 60477138  7559642                 
--- 
  
pairwise.t.test(p$Intensity,p$Treatment, p.adjust.method="none")
#   CC    CN    NC   
#CN 0.042 -     -    
#NC 0.807 0.026 -    
#NN 0.433 0.179 0.308

pairwise.t.test(p$Intensity,p$Treatment, p.adjust.method="fdr")
#   CC   CN   NC  
#CN 0.043 -     -    
#NC 0.638 0.058 -    
#NN 0.234 0.234 0.338


pairwise.t.test(p$Intensity,p$Treatment, p.adjust.method="bonferroni")

########################################
# PGK(2) Reanalyzed PGK3.22 on 5.11.16
p<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_PGK3.txt", header=T) 

leveneTest(Intensity ~ Site*Origin, data=p)

a <- aov(Intensity~Treatment,data=p)
summary(a)
#            Df   Sum Sq Mean Sq F value Pr(>F)
#Treatment    3 17746248 5915416   2.077  0.144
#Residuals   16 45558132 2847383               
TukeyHSD(x=a,conf.level=0.95)
#           diff        lwr       upr     p adj
#CN-CC 1544.7608 -1508.5697 4598.091 0.4898710
#NC-CC 1081.7621 -1971.5684 4135.093 0.7439791
#NN-CC 2623.7621  -429.5684 5677.093 0.1056425
#NC-CN -462.9987 -3516.3293 2590.332 0.9717545
#NN-CN 1079.0013 -1974.3293 4132.332 0.7454245
#NN-NC 1542.0000 -1511.3306 4595.331 0.4913465

## 4.
p<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_PGK4.txt", header=T) 

leveneTest(Intensity ~ Site*Origin, data=p)

a <- aov(Intensity~Treatment,data=p)
summary(a)
#            Df   Sum Sq Mean Sq F value Pr(>F)
#Treatment    3 17746248 5915416   2.077  0.144
#Residuals   16 45558132 2847383               
TukeyHSD(x=a,conf.level=0.95)

######################################################################
# Reanalyzed on 7.15.16
p<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_PGK3-1.txt", header=T) 

leveneTest(Intensity ~ Site*Origin, data=p)

a <- aov(Intensity~Treatment,data=p)
summary(a)
#            Df   Sum Sq Mean Sq F value Pr(>F)
#Treatment    3 17746248 5915416   2.077  0.144
#Residuals   16 45558132 2847383               
TukeyHSD(x=a,conf.level=0.95)

p<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_PGK3-2.txt", header=T) 

leveneTest(Intensity ~ Site*Origin, data=p)

a <- aov(Intensity~Treatment,data=p)
summary(a)
#            Df   Sum Sq Mean Sq F value Pr(>F)
#Treatment    3 17746248 5915416   2.077  0.144
#Residuals   16 45558132 2847383               
TukeyHSD(x=a,conf.level=0.95)


######################################################################
######################################################################
## 3. CYP1A##
cyp<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_CYP2.txt", header=T) 
str(cyp)

anova_cyp1<-aov(Intensity~Site+Origin,data=cyp)
summary(anova_cyp1)
#            Df    Sum Sq   Mean Sq F value  Pr(>F)   
#Site         1    291829   291829   0.024  0.879
#Origin       1   2022945  2022945   0.165  0.690
#Residuals   17 208377733 12257514              
---
op<-par(mfrow=c(2,2))
plot(anova_cyp1)

anova_cyp2<-aov(Intensity~Treatment,data=cyp)
summary(anova_cyp2)
#            Df    Sum Sq  Mean Sq F value Pr(>F)
#Treatment    3   9154476  3051492   0.242  0.866
#Residuals   16 201538030 12596127 

anova_cyp2<-aov(Intensity~Site*Origin,data=cyp)
summary(anova_cyp2)
#           Df    Sum Sq   Mean Sq F value  Pr(>F)   
#Site         1    291829   291829   0.023  0.881
#Origin       1   2022945  2022945   0.161  0.694
#Site:Origin  1   6839703  6839703   0.543  0.472
#Residuals   16 201538030 12596127     

op<-par(mfrow=c(2,2)); 
plot(anova_cyp2); 
par(op)


cypN<-cyp[1:10,]
cypC<-cyp[11:20,]

### Site effect on Genotype N?
anova_cypN<-aov(Intensity~Site,data=cypN)
summary(anova_cypN)
#            Df   Sum Sq  Mean Sq F value Pr(>F)  
#Site         1  2152960 2152960   0.342  0.575
#Residuals    8 50343000 6292875 

#No genotype effects at Site CW.
anova_cypC<-aov(Intensity~Site,data=cypC)
summary(anova_cypC)
#            Df   Sum Sq  Mean Sq F value Pr(>F)  
#Site         1   4978572  4978572   0.263  0.622
#Residuals    8 151195030 18899379  

pairwise.t.test(cyp$Intensity,cyp$Treatment, p.adjust.method="none")
#   CC    CN    NC   
#CN 0.54 -    -   
#NC 0.43 0.86 -   
#NN 0.70 0.82 0.68
pairwise.t.test(cyp$Intensity,cyp$Treatment, p.adjust.method="fdr")



################################################################
###################################################################
## 4. Ferrochelatase (compare across the membranes adjusting to C-> N or N->C )

f<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_Ferrochelatase.txt", header=T) 
str(f)

#check assumption
library(Rcmdr)
leveneTest(Intensity ~ Origin, data=f)
leveneTest(Intensity ~ Origin*Site, data=f)
leveneTest(Intensity ~ Site, data=f)
#ANOVA
anova_f<-aov(Intensity~Site*Origin, data=f)
summary(anova_f)
#            Df    Sum Sq   Mean Sq F value   Pr(>F)    
#Site         1 135933361 135933361   4.655 0.0465 * 
#Origin       1 273229908 273229908   9.357 0.0075 **
#Site:Origin  1 154897103 154897103   5.305 0.0350 * 
#Residuals   16 467196589  29199787  

########
pairwise.t.test(f$Intensity,f$Treatment, p.adjust.method="none")
#   CC    CN    NC   
#CN 0.9193 -      -     
#NC 0.0016 0.0020 -     
#NN 0.5329 0.6004 0.006
pairwise.t.test(f$Intensity,f$Treatment, p.adjust.method="fdr")
#   CC    CN    NC   
#CN 0.919 -     -    
#NC 0.006 0.006 -    
#NN 0.720 0.720 0.012



#### pairwise comparison of treatmnet ###
a1 <- aov(Intensity~Treatment,data=f)
posthoc <- TukeyHSD(x=a1,conf.level=0.95)
posthoc
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






################################################################
################################################################
# 5. Transgelin (only compare within genotype)
t<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_transgelin.txt", header=T) 
str(t)

t<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_transgelin2.txt", header=T) 



#check assumption
library(Rcmdr)
leveneTest(Intensity ~ Origin, data=t)
leveneTest(Intensity ~ Origin*Site, data=t)
leveneTest(Intensity ~ Site, data=t)


## Separate by genotype ###
tN<-t[1:10,]
tC<-t[11:20,]

leveneTest(Intensity ~ Site, data=tN)
leveneTest(Intensity ~ Site, data=tC)


anova_tN<-aov(Intensity~Site,data=tN)
summary(anova_tN)
#            Df   Sum Sq  Mean Sq F value Pr(>F)  
#Site         1 53038090 53038090   7.352 0.0266 *
#Residuals    8 57714120  7214265 

## Transgelin2
#            Df   Sum Sq  Mean Sq F value Pr(>F)  
#Site         1 93208090 93208090   7.771 0.0236 *
#Residuals    8 95950120 11993765 

#genotype effects at Site CW.
anova_tC<-aov(Intensity~Site,data=tC)
summary(anova_tC)
#            Df   Sum Sq  Mean Sq F value Pr(>F)  
#Site         1  42849000 42849000   0.692  0.429
#Residuals    8 495072000 61884000 

## Transgelin2
#            Df    Sum Sq  Mean Sq F value Pr(>F)
#Site         1   3025000  3025000   0.045  0.838
#Residuals    8 542584000 67823000

pairwise.t.test(t$Intensity,t$Treatment, p.adjust.method="none")



#######################################################################
# 6. Calmodulin binding protein 55kDa

cam<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_CaM55.txt", header=T) 
str(cam)

#check assumption
library(Rcmdr)
leveneTest(Intensity ~ Site, data=cam)

## Separate by transplant site ###
cN<-cam[1:10,]
cC<-cam[11:20,]

leveneTest(Intensity ~ Origin, data=cN)
leveneTest(Intensity ~ Origin, data=cC)


anova_cN<-aov(Intensity~Origin,data=cN)
summary(anova_cN)
#            Df   Sum Sq  Mean Sq F value Pr(>F)  
#Origin       1 2.633e+09 2.633e+09   1.295  0.288
#Residuals    8 1.627e+10 2.033e+0                

#genotype effects at Site CW.
anova_cC<-aov(Intensity~Origin,data=cC)
summary(anova_cC)
#            Df   Sum Sq  Mean Sq F value Pr(>F)  
#Origin       1 3.834e+06   3833802   0.029  0.868
#Residuals    8 1.043e+09 13043435

## 2-way anova
a1 <- aov(Intensity~Origin + Site,data=cam)
summary(a1)
#            Df    Sum Sq   Mean Sq F value  Pr(>F)   
#Origin       1 1.419e+09 1.419e+09   1.302 0.26970   
#Site         1 9.556e+09 9.556e+09   8.767 0.00875 **
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


plot(a2)

pairwise.t.test(cam$Intensity,cam$Treatment, p.adjust.method="none")
#   CC    CN    NC   
#CN 0.953 -     -    
#NC 0.195 0.215 -    
#NN 0.010 0.012 0.138

pairwise.t.test(cam$Intensity,cam$Treatment, p.adjust.method="fdr")
#  CC    CN    NC   
#CN 0.953 -     -    
#NC 0.258 0.258 -    
#NN 0.035 0.035 0.258




#######################################################################
#######################################################################






#7. Actin #2 (adjusted to N->C )
a<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_actin.txt", header=T) 

## Separate by Genotype(Origin) ###
aN<-a[1:10,]
aC<-a[11:20,]

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



### Any Site effects on N genotype?
anova_aN<-aov(Intensity~Site,data=aN)
summary(anova_aN)
#            Df   Sum Sq  Mean Sq F value Pr(>F)  
#Site         1  13202010 13202010   0.671  0.436
#Residuals    8 157349000 19668625 

#Any Site effects on C genotype?
anova_aC<-aov(Intensity~Site,data=aC)
summary(anova_aC)
#            Df   Sum Sq  Mean Sq F value Pr(>F)  
#Site         1  18803171 18803171   1.103  0.324
#Residuals    8 136381643 17047705  

aT1<-a[1:5,]
aT2<-a[11:15,]
aTT<-rbind(aT1,aT2)
anova_aTT<-aov(Intensity~Origin,data=aTT)
summary(anova_aTT)
#            Df   Sum Sq Mean Sq F value Pr(>F)
#Origin       1  26816595 26816595   1.126   0.32
#Residuals    8 190510220 23813777                  
--- 
  
pairwise.t.test(a$Intensity,a$Treatment, p.adjust.method="none")
#   CC    CN    NC   
#CN 0.327 -     -    
#NC 0.847 0.244 -    
#NN 0.312 0.056 0.409

pairwise.t.test(a$Intensity,a$Treatment, p.adjust.method="fdr")
#   CC   CN   NC  
#CN 0.043 -     -    
#NC 0.638 0.058 -    
#NN 0.234 0.234 0.338


pairwise.t.test(p$Intensity,p$Treatment, p.adjust.method="bonferroni")



#######################################################################
#######################################################################

## 8. SOD 5.18.16 Data  ####

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

##### 1-way anova

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



###  Perform either pairwise t. test with some correction, or Turkey HSD Test to compare which one differ which)
pairwise.t.test(SOD$Intensity,SOD$Treatment, p.adjust.method="none")
#   CC    CN    NC   
#CN 0.2841 -      -     
#NC 0.9333 0.2501 -     
#NN 0.0094 0.0838 0.0078
pairwise.t.test(SOD$Intensity,SOD$Treatment, p.adjust.method="fdr")
#   CC     CN     NC    
#CN 0.341 -     -    
#NC 0.933 0.341 -    
#NN 0.028 0.168 0.028

#SOD2 Results:
#   CC     CN     NC    
#CN 0.0482 -      -     
#NC 0.2442 0.2442 -     
#NN 0.0023 0.0976 0.0147




### Turkey's HSD Test
a1 <- aov(Intensity~Site*Origin,data=SOD)
posthoc <- TukeyHSD(x=a1,conf.level=0.95)
posthoc
#$`Site:Origin`
#                     diff        lwr       upr     p adj
#SiteN:C-CW:C     20639.873 -32641.340  73921.09 0.6897617
#CW:N-CW:C        -1583.318 -54864.531  51697.90 0.9997699
#SiteN:N-CW:C     54976.682   1695.469 108257.90 0.0419512 *
#CW:N-SiteN:C    -22223.191 -75504.404  31058.02 0.6395669
#SiteN:N-SiteN:C  34336.809 -18944.404  87618.02 0.2900912
#SiteN:N-CW:N     56560.000   3278.787 109841.21 0.0355508 *

#SOD2 RESULTS:

#$`Site:Origin`
#                     diff        lwr       upr     p adj
#SiteN:C-CW:C     43179.79  -6411.830  92771.40 0.0996328 
#CW:N-CW:C        20956.60 -28635.021  70548.21 0.6302132
#SiteN:N-CW:C     77516.60  27924.979 127108.21 0.0019635 **
#CW:N-SiteN:C    -22223.19 -71814.807  27368.43 0.5866509
#SiteN:N-SiteN:C  34336.81 -15254.807  83928.43 0.2357535
#SiteN:N-CW:N     56560.00   6968.384 106151.62 0.0227514 *


##Comparing G8-CW only
sod<-read.table("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_SOD_CW.txt", header=T)  #C->C normalized to O->N samples
anova_sod<-aov(Intensity~Site,data=sod)
summary(anova_sod)
#            Df    Sum Sq   Mean Sq F value Pr(>F)
#Site         1 1.849e+08 184900000   0.519  0.492
#Residuals    8 2.852e+09 356525500


################################################################

#read data
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

a1<- aov(Intensity~Origin + Site,data=h60)
summary(a1)
#            Df    Sum Sq  Mean Sq F value  Pr(>F)   
#Origin       1  44096112 44096112   5.648 0.03030 * 
#Site         1  75506599 75506599   9.670 0.00674 **
#Residuals   16 124928543  7808034