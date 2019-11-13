rm(list=ls(all=TRUE))

rte1 <- read.csv("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_corallite1.csv", header=T)
library("car")
scatterplotMatrix(rte1[3:21])

scatterplotMatrix(~ X1 + X2 + X3 + X4 + X5+ X6+X7+X8+X9+X10| Site,
                  data=rte1, smooth=FALSE, reg.line=FALSE, ellipse=TRUE,
                  by.groups=TRUE, diagonal="none")
scatterplotMatrix(~ X11 + X12 + X13 + X14 + X15+ X16+X17+X18+X19| Site,
                  data=rte1, smooth=FALSE, reg.line=FALSE, ellipse=TRUE,
                  by.groups=TRUE, diagonal="none")

par(mfrow=c(2, 2))
for (response in c("X8")) + Boxplot(rte1[,response] ~ Site, data=rte1, ylab=response)

plot(rte1$X1,rte1$X2)
text(rte1$X1, rte1$X2, rte1$Site, cex=0.7, pos=4, col="red")
par(mfrow=c(1, 1))

### model fitting
mod.rte1<- lm(cbind(X1, X2, X3,X4, X5, X6,X7,X8, X9, X10,X11, X12, X13,X14,X15,X16,X17,X18,X19)~ Site, data=rte1)

class(mod.rte1)
mod.rte1
summary(mod.rte1)

(manova.rte1 <- Anova(mod.rte1))
#Type II MANOVA Tests: Pillai test statistic
#Df test stat approx F num Df den Df    Pr(>F)    
#Site  1   0.55144   5.1763     19     80 7.999e-08 ***

summary(manova.rte1)


### CANDISC ###

library("candisc")
library("heplots")
rte1.can<-candisc(mod.rte1, term="Samples")
plot(rte1.can, scale=6, fill=TURE)

rte1.can2<-candisc(mod.rte1, data=rte1)
plot(rte1.can2, scale=8, fill=TURE)

rte1.can3 <- candisc(mod.rte1, data=rte1, ndim=1)
plot(rte1.can3)


heplot(rte1.can, scale=6, fill=TRUE)





calcWithinGroupsVariance <- function(variable,groupvariable)
{
  # find out how many values the group variable can take
  groupvariable2 <- as.factor(groupvariable[[1]])
  levels <- levels(groupvariable2)
  numlevels <- length(levels)
  # get the mean and standard deviation for each group:
  numtotal <- 0
  denomtotal <- 0
  for (i in 1:numlevels)
  {
    leveli <- levels[i]
    levelidata <- variable[groupvariable==leveli,]
    levelilength <- length(levelidata)
    # get the standard deviation for group i:
    sdi <- sd(levelidata)
    numi <- (levelilength - 1)*(sdi * sdi)
    denomi <- levelilength
    numtotal <- numtotal + numi
    denomtotal <- denomtotal + denomi
  }
  # calculate the within-groups variance
  Vw <- numtotal / (denomtotal - numlevels)
  return(Vw)
}

calcWithinGroupsVariance(rte1[3],rte1[2])
calcWithinGroupsVariance(rte1[4],rte1[2])


calcBetweenGroupsVariance <- function(variable,groupvariable)
{
  # find out how many values the group variable can take
  groupvariable2 <- as.factor(groupvariable[[1]])
  levels <- levels(groupvariable2)
  numlevels <- length(levels)
  # calculate the overall grand mean:
  grandmean <- mean(variable)
  # get the mean and standard deviation for each group:
  numtotal <- 0
  denomtotal <- 0
  for (i in 1:numlevels)
  {
    leveli <- levels[i]
    levelidata <- variable[groupvariable==leveli,]
    levelilength <- length(levelidata)
    # get the mean and standard deviation for group i:
    meani <- mean(levelidata)
    sdi <- sd(levelidata)
    numi <- levelilength * ((meani - grandmean)^2)
    denomi <- levelilength
    numtotal <- numtotal + numi
    denomtotal <- denomtotal + denomi
  }
  # calculate the between-groups variance
  Vb <- numtotal / (numlevels - 1)
  Vb <- Vb[[1]]
  return(Vb)
}

calcBetweenGroupsVariance(rte1[3],rte1[2])  #not working

# Standardizing the variables
rte_s<-as.data.frame(scale(rte1[3:21]))
#check the mean and variacne
sapply(rte_s,mean)
sapply(rte_s,sd)

#run PCA
rte1.pca<-prcomp(rte_s)
summary(rte1.pca)
rte1.pca$sdev  #SD of each component

#Deciding How Many Principal Components to Retain
screeplot(rte1.pca, type="lines")

# or Kaiser’s criterion: only retain principal components for which the variance is above 1 (when principal component analysis was applied to standardised data)
(rte1.pca$sdev)^2

#Loadings for the Principal Components
#for PC1
rte1.pca$rotation[,1]

plot(rte1.pca$x[,1],rte1.pca$x[,2]) 
text(rte1.pca$x[,1],rte1.pca$x[,2], rte1$Site, cex =0.6, pos=4, col="blue")

plot(rte1.pca$x[,1],rte1.pca$x[,3]) 
text(rte1.pca$x[,1],rte1.pca$x[,3], rte1$Samples, cex =0.6, pos=4, col="blue")

plot(rte1.pca$x[,2],rte1.pca$x[,3]) 
text(rte1.pca$x[,2],rte1.pca$x[,3], rte1$Samples, cex =0.6, pos=4, col="blue")



library("MASS") 
#rte1.lda <- lda(rte1$Site ~ rte1$X1 + rte1$X2+ rte1$X3+ rte1$X4+ rte1$X5+ rte1$X6+ rte1$X7+ rte1$X8+ rte1$X9+ rte1$X10+ rte1$X11+ rte1$X12+ rte1$X13 + rte1$X14 + rte1$X15 + rte1$X16 + rte1$X17+ rte1$X18 + rte1$X19)  

rte1.lda <- lda(Site ~ X1 + X2+ X3+ X4+ X5+ X6+ X7+ X8+ X9+ X10+ X11+ X12+ X13 +X14 +X15 +X16 +X17+X18 +X19, data =rte1)  
rte1.lda
rte1.lda.values <- predict(rte1.lda)
ldahist(data = rte1.lda.values$x[,], g=class) #doesn't work
plot(rte1.lda.values$x[],rte1.lda.values$x[])
text(rte1.lda.values$x[],rte1.lda.values$x[], class,cex=0.7,pos=4,col="red")


#######################################
#############  calculate SD  ########
rm(list=ls(all=TRUE))
rte1 <- read.csv("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_corallite1.csv", header=T)

p<-c(3:21)

n1sd<-numeric()
for (i in 3:21) {
  n1sd[i]<-sd(rte1[1:10,i])
}
n2sd<-numeric()
for (i in 3:21) {
  n2sd[i]<-sd(rte1[11:20,i])
}
n3sd<-numeric()
for (i in 3:21) {
  n3sd[i]<-sd(rte1[21:30,i])
}
n4sd<-numeric()
for (i in 3:21) {
  n4sd[i]<-sd(rte1[31:40,i])
}
n5sd<-numeric()
for (i in 3:21) {
  n5sd[i]<-sd(rte1[41:50,i])
}
SD<-cbind(n1sd,n2sd,n3sd,n4sd,n5sd)


#Average
n1av<-numeric()
for (i in 3:21) {
  n1av[i]<-ave(rte1[1:10,i])
}
n2av<-numeric()
for (i in 3:21) {
  n2av[i]<-ave(rte1[11:20,i])
}
n3av<-numeric()
for (i in 3:21) {
  n3av[i]<-ave(rte1[21:30,i])
}

n4av<-numeric()
for (i in 3:21) {
  n4av[i]<-ave(rte1[31:40,i])
}

n5av<-numeric()
for (i in 3:21) {
  n5av[i]<-ave(rte1[41:50,i])
}

AVE<-cbind(n1av,n2av,n3av,n4av,n5av)
N<-cbind(AVE,SD)
write.csv(N, file="/Users/kahotisthammer/Dropbox/R/N.csv", row.names = TRUE, col.names = TRUE)

######          
c1sd<-numeric()
for (i in 3:21) {
  c1sd[i]<-sd(rte1[51:60,i])
}
c2sd<-numeric()
for (i in 3:21) {
  c2sd[i]<-sd(rte1[61:70,i])
}
c3sd<-numeric()
for (i in 3:21) {
  c3sd[i]<-sd(rte1[71:80,i])
}
c4sd<-numeric()
for (i in 3:21) {
  c4sd[i]<-sd(rte1[81:90,i])
}
c5sd<-numeric()
for (i in 3:21) {
  c5sd[i]<-sd(rte1[91:100,i])
}
SD2<-cbind(c1sd,c2sd,c3sd,c4sd,c5sd)


#Average
c1av<-numeric()
for (i in 3:21) {
  c1av[i]<-ave(rte1[51:60,i])
}
c2av<-numeric()
for (i in 3:21) {
  c2av[i]<-ave(rte1[61:70,i])
}
c3av<-numeric()
for (i in 3:21) {
  c3av[i]<-ave(rte1[71:80,i])
}
c4av<-numeric()
for (i in 3:21) {
  c4av[i]<-ave(rte1[81:90,i])
}
c5av<-numeric()
for (i in 3:21) {
  c5av[i]<-ave(rte1[91:100,i])
}

AVE2<-cbind(c1av,c2av,c3av,c4av,c5av)
C<-cbind(AVE2,SD2)
write.csv(C, file="/Users/kahotisthammer/Dropbox/R/C.csv", row.names = TRUE)


## Using 5 corallites /sample:
#SD 
n1sd<-numeric()
for (i in 3:21) {
  n1sd[i]<-sd(rte1[1:5,i])
}
n2sd<-numeric()
for (i in 3:21) {
  n2sd[i]<-sd(rte1[11:15,i])
}
n3sd<-numeric()
for (i in 3:21) {
  n3sd[i]<-sd(rte1[21:25,i])
}
n4sd<-numeric()
for (i in 3:21) {
  n4sd[i]<-sd(rte1[31:35,i])
}
n5sd<-numeric()
for (i in 3:21) {
  n5sd[i]<-sd(rte1[41:45,i])
}
SD<-cbind(n1sd,n2sd,n3sd,n4sd,n5sd)

#Average
n1av<-numeric()
for (i in 3:21) {
  n1av[i]<-ave(rte1[1:5,i])
}
n2av<-numeric()
for (i in 3:21) {
  n2av[i]<-ave(rte1[11:15,i])
}
n3av<-numeric()
for (i in 3:21) {
  n3av[i]<-ave(rte1[21:25,i])
}

n4av<-numeric()
for (i in 3:21) {
  n4av[i]<-ave(rte1[31:35,i])
}

n5av<-numeric()
for (i in 3:21) {
  n5av[i]<-ave(rte1[41:45,i])
}

AVE<-cbind(n1av,n2av,n3av,n4av,n5av)
N_5<-cbind(AVE,SD)
write.csv(N_5, file="/Users/kahotisthammer/Dropbox/R/N_5.csv", row.names = TRUE, col.names = TRUE)

######          
c1sd<-numeric()
for (i in 3:21) {
  c1sd[i]<-sd(rte1[51:55,i])
}
c2sd<-numeric()
for (i in 3:21) {
  c2sd[i]<-sd(rte1[61:65,i])
}
c3sd<-numeric()
for (i in 3:21) {
  c3sd[i]<-sd(rte1[71:75,i])
}
c4sd<-numeric()
for (i in 3:21) {
  c4sd[i]<-sd(rte1[81:85,i])
}
c5sd<-numeric()
for (i in 3:21) {
  c5sd[i]<-sd(rte1[91:95,i])
}
SD2<-cbind(c1sd,c2sd,c3sd,c4sd,c5sd)


#Average
c1av<-numeric()
for (i in 3:21) {
  c1av[i]<-ave(rte1[51:55,i])
}
c2av<-numeric()
for (i in 3:21) {
  c2av[i]<-ave(rte1[61:65,i])
}
c3av<-numeric()
for (i in 3:21) {
  c3av[i]<-ave(rte1[71:75,i])
}
c4av<-numeric()
for (i in 3:21) {
  c4av[i]<-ave(rte1[81:85,i])
}
c5av<-numeric()
for (i in 3:21) {
  c5av[i]<-ave(rte1[91:95,i])
}

AVE2<-cbind(c1av,c2av,c3av,c4av,c5av)
C_5<-cbind(AVE2,SD2)
write.csv(C_5, file="/Users/kahotisthammer/Dropbox/R/C_5.csv", row.names = TRUE)


