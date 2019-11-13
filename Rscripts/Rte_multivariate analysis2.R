rm(list=ls(all=TRUE))

rte1 <- read.csv("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_corallite1.csv", header=T)
library("car")

###  Multivariate Normality Test using MVN  ###########
library(MVN)
rte2<-rte1[,3:21]
#  1. Mardia’s multivariate skewness and kurtosis coefficients
test1<-mardiaTest(rte2)
test1
#   Mardia's Multivariate Normality Test 
#--------------------------------------- 
#  data : rte2 

#g1p            : 102.7123 
#chi.skew       : 1711.871 
#p.value.skew   : 4.969561e-12 

#g2p            : 413.0736 
#z.kurtosis     : 2.490999 
#p.value.kurt   : 0.01273844 

#chi.small.skew : 1768.481 
#p.value.small  : 5.516629e-15 

#Result          : Data are not multivariate normal. 
#--------------------------------------- 

test2<-hzTest(rte2)
test2


#### Check normality for individual dependent variables  ######
library(Rcmdr)
leveneTest(X1 ~ Samples, data=rte1)
leveneTest(X2 ~ Samples, data=rte1)
leveneTest(X3 ~ Samples, data=rte1)
leveneTest(X4 ~ Samples, data=rte1)
leveneTest(X5 ~ Samples, data=rte1)
leveneTest(X6 ~ Samples, data=rte1)
leveneTest(X7 ~ Samples, data=rte1)
#Levene's Test for Homogeneity of Variance (center = median)
#      Df F value    Pr(>F)    
#group  9  4.2793 0.0001204 ***
#      90                      

leveneTest(X8 ~ Samples, data=rte1)
#      Df F value  Pr(>F)  
#group  9  2.3684 0.01883 *

leveneTest(X9 ~ Samples, data=rte1)
leveneTest(X10 ~ Samples, data=rte1)
leveneTest(X11 ~ Samples, data=rte1)

leveneTest(X12 ~ Samples, data=rte1)
leveneTest(X13 ~ Samples, data=rte1)
leveneTest(X14 ~ Samples, data=rte1)
leveneTest(X15 ~ Samples, data=rte1)
leveneTest(X16 ~ Samples, data=rte1)
leveneTest(X17 ~ Samples, data=rte1)
leveneTest(X18 ~ Samples, data=rte1)       
#      Df F value Pr(>F)  
#group  9   1.764 0.0862 .
#90 
leveneTest(X19 ~ Samples, data=rte1)

## X7, X8, X18 don't have homogeneous variance ##
## Data trnasformation using boxcox:

# 1.X7
leveneTest(log(X7) ~ Samples, data=rte1) 
#        Df F value  Pr(>F)  
#  group  9  2.0085 0.04721 

library(MASS)            
boxcox(X7 ~ Samples+Site, data=rte1)
boxcox(X7 ~ Samples, data=rte1)
#group  9  2.0085 0.04721 *
leveneTest(X7^(0.009) ~ Samples, data=rte1)
#group  9  2.0079 0.04729 *

## X8
leveneTest(log(X8) ~ Samples, data=rte1)
# p=0.02707 *
boxcox(X8 ~ Samples, data=rte1,lambda=seq(-0.1,1,0.01))
leveneTest(X8^(0.66) ~ Samples, data=rte1)
#group  9  2.3301 0.0208 *
leveneTest(sqrt(X8) ~ Samples, data=rte1)
# p=0.02192 *
leveneTest(X8^4 ~ Samples, data=rte1)

## X18
leveneTest(log(X18) ~ Samples, data=rte1)
# group  9  0.8418 0.5798*



#####  convert three variables #####
rtec<-rte1

rtec[,9]<-(-log(rte1[,9]))
rtec[,10]<-(rte1[,10])^(0.65)
rtec[,20]<-(log(rte1[,20]))

library(MVN)
rtec2<-rtec[,3:21]
#  1. Mardia’s multivariate skewness and kurtosis coefficients
testc<-mardiaTest(rtec2)
testc

rte.t<-rtec2[,-20]
rte.t<-rte.t[,-10]
rte.t<-rte.t[,-9]

testt<-mardiaTest(rte.t)
testt



### run MANOVA, Discriminant analysis etc.  ###
library("car")
### model fitting
mod.rtec<- lm(cbind(X1, X2, X3,X4, X5, X6,X7,X8, X9, X10,X11, X12, X13,X14,X15,X16,X17,X18,X19)~ Site, data=rtec)
mod.rtec
summary(mod.rtec)
#lm(formula = X19 ~ Site, data = rtec)

#Residuals:
#  Min         1Q     Median         3Q        Max 
#-0.0242200 -0.0063867  0.0000325  0.0058671  0.0294467 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 0.045715   0.001440  31.737   <2e-16 ***
#  SiteN       0.003505   0.002037   1.721   0.0885 . 
manova.rtec <- Anova(mod.rtec)
summary(manova.rtec)
#Multivariate Tests: Site
#                  Df test stat approx F num Df den Df     Pr(>F)    
#Pillai            1 0.5513767 5.174911     19     80 8.0381e-08 ***
#Wilks             1 0.4486233 5.174911     19     80 8.0381e-08 ***
#Hotelling-Lawley  1 1.2290415 5.174911     19     80 8.0381e-08 ***
#Roy               1 1.2290415 5.174911     19     80 8.0381e-08 ***

##Canonical discriminant analysis ###

library("candisc")
library("heplots")
rtec.can<-candisc(mod.rtec)
plot(rtec.can, scale=6)

#rtec.can2<-candisc(mod.rtec, data=rtec)
#plot(rtec.can2, scale=8, fill=TURE)


heplot(rtec.can, scale=6, fill=TRUE)


### Canonical Discrimannt Analysis ###
library("candisc")
library("heplots")
mod.rte1<- lm(cbind(X1, X2, X3,X4, X5, X6,X7,X8, X9, X10,X11, X12, X13,X14,X15,X16,X17,X18,X19)~ Site, data=rtec)

rte.can<-candisc(mod.rte1)
plot(rte.can) 
plot(rte.can, scale=)

rte1.can2<-candisc(mod.rte1, data=rad1)
plot(rad1.can2, scale=8)

rte.can3 <- candisc(mod.rte1, data=rte1, ndim=2)
plot(rte.can3)


heplot(rte.can, scale=6, fill=TRUE)



