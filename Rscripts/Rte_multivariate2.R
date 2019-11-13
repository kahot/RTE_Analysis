rm(list=ls(all=TRUE)
library("car")
library("candisc")
library("heplots")
library(vegan)
library(lme4)

rte2 <- read.csv("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_corallite2.csv", header=T)
#scatterplotMatrix(rte2[3:28])
str(rte2)
as.integer(rte2$denticle_no)->rte2$denticle_no
as.integer(rte2$pali_height)->rte2$pali_height
#as.factor(rte2$lateral_pairs)->rte2$plateral_pairs
#as.factor(rte2$triplet_spacing)->rte2$triplet_spacing
#as.factor(rte2$corallite_shape )->rte2$corallite_shape
#as.factor(rte2$columella_shape)->rte2$columella_shape


scatterplotMatrix(~ pali_no + lateral_pairs + triplet_spacing + pali_height + corallite_shape + denticle_no + columella_shape | Site,
                  data=rte2, smooth=FALSE, reg.line=FALSE, ellipse=TRUE,
                  by.groups=TRUE, diagonal="none")


mod.rte2<- lm(cbind(X1, X2, X3,X4, X5, X6,X7,X8, X9, X10,X11, X12, X13,X14,X15,X16,X17,X18,X19, pali_no, lateral_pairs, triplet_spacing,pali_height, corallite_shape, denticle_no, columella_shape)~ Site, data=rte2)
class(mod.rte2)
mod.rte2
summary(mod.rte2)
manova.rte2 <- Anova(mod.rte2)
manova.rte2
#Type II MANOVA Tests: Pillai test statistic
#Df test stat approx F num Df den Df   Pr(>F)    
#Site  1   0.71655   7.0977     26     73 1.74e-11 ***


can.rte2<-candisc(mod.rte2, term="Site")
plot(can.rte2, scale=6, fill=TURE)

#One-way MANOVA
manova.rte2<- manova(cbind(X1, X2, X3,X4, X5, X6,X7,X8, X9, X10,X11, X12, X13,X14,X15,X16,X17,X18,X19, pali_no, lateral_pairs, triplet_spacing,pali_height, corallite_shape, denticle_no, columella_shape)~ Site, data=rte2)
summary(manova.rte2)
#          Df  Pillai approx F num Df den Df    Pr(>F)    
#Site       1 0.70731   6.7852     26     73 4.831e-11 ***
#Residuals 98 




#########

canL.rte2 <-candiscList(mod.rte2)
names(canL.rte2)
plot(canL.rte2, type="n", ask=FALSE)
heplot(canL.rte2$Site, scale=2)
plot.candisc(canL.rte2)


## http://strata.uga.edu/6370/lecturenotes/discriminantFunctionAnalysis.html
## Discriminant Function Analysis:

pairs(rte2[,3:8])
pairs(rte2[,9:13])
pairs(rte2[,14:19])
pairs(rte2[,20:24])
pairs(rte2[,25:28])

## Discriminant Analysis w/o transformation ##

lda.rte2<-lda(Site ~ X1+  X2+  X3+ X4+  X5+  X6+ X7+ X8+  X9+  X10+ X11+  X12+  X13+ X14+ X15+ X16+ X17+ X18+ X19+  pali_no+  lateral_pairs+  triplet_spacing+ pali_height+  corallite_shape+  denticle_no+  columella_shape, data=rte2)
lda.rte2
ldap.rte2<-predict(lda.rte2)
ldap.rte2

plot(lda.rte2,type="both")

#evaluating the DFA
tab <- rbind(rte2$Site[2:101,], ldap.rte2$class)
tab  #pretty good fit
ldap.rte2$class


### NMDS ###
require(vegan)
rte2_data<-rte2[3:28]  #extract the variables only
rte2.dist=vegdist(rte2_data)
mds.rte1<-metaMDS(rte2.dist,distance = "bray", enginek=c("isoMDS"), autotransform=FALSE,trymax=50, wascores=TRUE,expand = TRUE )
plot(mds.rte1)
#Run 0 stress 0.2025448 
#Run 1 stress 0.2321227 
#Run 2 stress 0.2022664 
#... New best solution
#... procrustes: rmse 0.005038447  max resid 0.03089456 
#Run 3 stress 0.2025444 
#... procrustes: rmse 0.004979483  max resid 0.03082379 
#Run 4 stress 0.2022669 
#... procrustes: rmse 0.00101912  max resid 0.007137436 
#*** Solution reached
plot(mds.rte1)

N=which(rte2$Site=="N") 
C=which(rte2$Site=="C") 

Nrte1=mds.rte1$points[-C,]  
Crte1=mds.rte1$points[-N,]
plot(mds.rte1, type="n")  
points(Nrte1, col="red")
points(Crte1, col="blue")


mds.rte2<-metaMDS(rte2_data)
mds.rte2  
#Distance: bray 
#Dimensions: 2 
#Stress:     0.2022609 
#Stress type 1, weak ties
plot(mds.rte2, type="n")
points(mds.rte2, display="sites", cex=0.8, pch=21)
text(mds.rte2, display="species", cex=0,7, col="blue")

N=which(rte2$Site=="N") 
C=which(rte2$Site=="C") 

Ns=mds.rte2$points[-C,]  
Cs=mds.rte2$points[-N,]

#plot the ordination using text
ordiplot(mds.rte2, type="t", display="sites")	
points(Ns, col="red")
points(Cs, col="blue")


#### PCA ###
rte2_s<-as.data.frame(scale(rte2[3:28]))
#sapply(rte2_s,mean)
#sapply(rte2_s,sd)

rte2.pca<-prcomp(rte2_s)
summary(rte2.pca)
#rte2.pca$sdev  #SD of each component

#Deciding How Many Principal Components to Retain
screeplot(rte2.pca, type="lines")
(rte2.pca$sdev)^2 #Kaiser’s criterion  upto PC8 has variacne >1
rte2.pca$rotation[,1]  #loadings
rte2.pca$rotation[,2]  #loadigns for PC2

plot(rte2.pca$x[,1],rte2.pca$x[,2]) 
text(rte2.pca$x[,1],rte2.pca$x[,2], rte2$Site, cex =0.6, pos=4, col="blue")


text(rte2.pca$x[,1],rte2.pca$x[,2], rte2$Site, cex =0.6, pos=4, col="blue")


plot(rte2.pca$x[,1],rte2.pca$x[,3]) 
text(rte2.pca$x[,1],rte2.pca$x[,3], rte2$Site, cex =0.6, pos=4, col="blue")

plot(rte2.pca$x[,2],rte2.pca$x[,3]) 
text(rte2.pca$x[,2],rte2.pca$x[,3], rte2$Site, cex =0.6, pos=4, col="blue")

#########################################
#########################################
#########################################
Running each numerical variables to normality test.
library(Rcmdr)
leveneTest(X1 ~ Samples, data=rte2)

result <- vector("list",25)
for (i in 3:28)
{ result[[i]]<-leveneTest(rte2[,i] ~ rte2$Samples)
}
# Levene's Test for Homogeneity of Variance (center = median)
#         Df F value Pr(>F)
#X1 group  9  0.8889 0.5385
#X2 group  9  1.0748 0.3891
#3  group  9  0.4187  0.922
#4  group  9  1.6165 0.1224
#5  group  9  1.0284 0.4239
#6  group  9  1.0592 0.4006
#7  group  9  4.2793 0.0001204 ***
#8  group  9  2.3684 0.01883 *
#9  group  9  1.3462 0.2248
#10 group  9  0.7743 0.6403
#11 group  9  1.5669 0.1373
#12 group  9  1.0965 0.3735
#13 group  9  0.7393 0.6719
#14 group  9  1.1851 0.3142
#15 group  9   0.828 0.5921
#16 group  9  1.4245 0.1895
#17 group  9  0.9783 0.4634
#18 group  9   1.764 0.0862 .
#19 group  9  0.5195 0.8569

#20  group  9  2.1429 0.03363 *
#21 group  9  2.3853 0.01803 *
#22 group  9  5.0826 1.501e-05 ***
#23 group  9  7.2372 7.601e-08 ***
#24 group  9   1.801 0.07881 .
#25 group  9  2.1282 0.03491 *
#26 group  9  2.2422 0.02608 *

scatterplotMatrix(rte2[9:10])
#figure out the transformation method
library(MASS)            
boxcox(X7~ Samples, data =rte2)  # \lambda that maximizes the log-lik is very close to 0. for lambda=0, the transformation is log(x)
#log transformation of X7
logX7<-log(rte2$X7)
plot(logX7)
leveneTest(logX7 ~ rte2$Samples)
#Levene's Test for Homogeneity of Variance (center = median)
#Df F value  Pr(>F)  
#group  9  2.0085 0.04721 *


#transform ^0.2 then, log transformation
sqX7<-(rte2$X7)^0.2
tX7<-log((rte2$X7)^0.2)
plot(tX7)
leveneTest(tX7 ~ rte2$Samples)
#Levene's Test for Homogeneity of Variance (center = median)
#Df F value  Pr(>F)  
#group  9  2.0085 0.04721 *



boxcox((X7)^0.2 ~ Samples, data =rte2)
boxcox(logX7*(-1) ~ rte2$Samples)
leveneTest(sqrt(X7) ~ Samples, data =rte2)
leveneTest((X7)^0.2 ~ Samples, data =rte2)
#log is the best. 
rte.t<-rte2
rte.t[,9]<-logX7
rte.t<-rte.t[,1:21]

#scaleing and run PCA
rte.ts<-as.data.frame(scale(rte.t[3:21]))
rte.ts.pca<-prcomp(rte.ts)
plot(rte.ts.pca$x[,1],rte.ts.pca$x[,2]) 
text(rte.ts.pca$x[,1],rte.ts.pca$x[,2], rte.t$Site, cex =0.6, pos=4, col="blue")
biplot(rte.ts.pca)
summary(rte.ts.pca)

screeplot(rte.ts.pca, type="lines")
#Eigenvalues (variance) for each PC
#Kaiser-Guttmann criterion- axis with eigenvalue less than 1 should be rejected
(rte.ts.pca$sdev)^2
#[1] 5.21968331 2.75765801 1.81207229 1.52770457 1.23336814 0.98786463 0.79437693 0.76653465
#[9] 0.69058992 0.58466472 0.53172040 0.41525723 0.40008123 0.31296818 0.27662245 0.26033778
#[17] 0.17931535 0.14984172 0.09933847

#obtain loadings for the PC1 for 19 variables:
rte.ts.pca$rotation[,1]

#Plot N and C in separate colors
N=which(rte2$Site=="N") 
C=which(rte2$Site=="C") 

#Nonly=rte.ts.pca$x[-C,]  
#Conly=rte.ts.pca$x[-N,]
plot(rte.ts.pca$x[,1],rte.ts.pca$x[,2], ylab='PC2', xlab='PC1') 
text(rte.ts.pca$x[-C,1],rte.ts.pca$x[-C,2], rte.t$Sample[-C], cex =0.6, pos=4, col="red")
text(rte.ts.pca$x[-N,1],rte.ts.pca$x[-N,2], rte.t$Sample[-N], cex =0.6, pos=4, col="blue")

#Plot PC1 x PC3
plot(rte.ts.pca$x[,1],rte.ts.pca$x[,3], ylab='PC3', xlab='PC1') 
text(rte.ts.pca$x[-C,1],rte.ts.pca$x[-C,3], rte.t$Sample[-C], cex =0.6, pos=4, col="red")
text(rte.ts.pca$x[-N,1],rte.ts.pca$x[-N,3], rte.t$Sample[-N], cex =0.6, pos=4, col="blue")





#change the sample names to numbers only
write.csv(rte.t, file="/Users/kahotisthammer/Dropbox/R/RTE1/rte.t.csv", row.names = TRUE, col.names = TRUE)
rte.tn<-read.csv("/Users/kahotisthammer/Dropbox/R/RTE1/rte.t_numbered.csv")

plot(rte.ts.pca$x[,1],rte.ts.pca$x[,2]) 
text(rte.ts.pca$x[-C,1],rte.ts.pca$x[-C,2], rte.tn$Sample[-C], cex =0.6, pos=4, col="red")
text(rte.ts.pca$x[-N,1],rte.ts.pca$x[-N,2], rte.tn$Sample[-N], cex =0.6, pos=4, col="blue")

plot(rte.ts.pca$x[,1],rte.ts.pca$x[,2]) 
text(rte.ts.pca$x[-C,1],rte.ts.pca$x[-C,2], rte.t$Sample[-C], cex =0.6, pos=4, col="red")
text(rte.ts.pca$x[-N,1],rte.ts.pca$x[-N,2], rte.t$Sample[-N], cex =0.6, pos=4, col="blue")



plot(rte.ts.pca$x[,1],rte.ts.pca$x[,3]) 
text(rte.ts.pca$x[,1],rte.ts.pca$x[,3], rte.t$Site, cex =0.6, pos=4, col="blue")

plot(rte.ts.pca$x[,3],rte.ts.pca$x[,2]) 
text(rte.ts.pca$x[,3],rte.ts.pca$x[,2], rte.t$Site, cex =0.6, pos=4, col="blue")

#Different PCA plot
library(ggplot2)
scores = as.data.frame(rte.ts.pca$x)
ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text(colour = "red", alpha = 0.8, size = 4) +
  ggtitle("PCA plot of Corallite Morphometrics")

#one-wat MANOVA
mod.rte.t1<- lm(cbind(X1, X2, X3,X4, X5, X6,X7,X8, X9, X10,X11, X12, X13,X14,X15,X16,X17,X18,X19)~ Site, data=rte.t)
mod.rte.t1
summary(mod.rte.t1)

(manova.rte.t1 <- Anova(mod.rte.t1))
#Type II MANOVA Tests: Pillai test statistic
#     Df test stat approx F num Df den Df    Pr(>F)    
#Site  1   0.54867   5.1187     19     80 9.824e-08 ***

###### Plotting
library(stylo)
require(graphics)

samples<-rte2[1,]
write(samples, file="Users/kahotisthammer/Dropbox/R/RTE1/samples.txt",ncolumns=ncol(samples))
