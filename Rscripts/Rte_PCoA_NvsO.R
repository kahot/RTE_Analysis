rm(list=ls(all=TRUE))
library(vegan)
library(ape)
library(labdsv)
#library("car")
#library("candisc")
#library("heplots")
#library(lme4)

rte2 <- read.csv("/Users/kahotisthammer/Dropbox/R/RTE_Analysis/RTE_corallite2.csv", header=T)
str(rte2)
as.integer(rte2$denticle_no)->rte2$denticle_no
as.integer(rte2$pali_height)->rte2$pali_height
str(rte2)
rte2[,23]+1 -> rte2[,23]
rte2[,24]+1 -> rte2[,24]

rte3<-rte2[,3:28]

#dist(rte3, method="euclidean", diag=FALSE, upper=FALSE, p=2) ->rtedist
#rtepcoa<-pcoa(rtedist, correction="none", rn=NULL)
#biplot(rtepcoa)

#Euclidean distance
#rted.eu<-vegdist(rte3, method="euclidean", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
dis.euc <- dist(rte3,'euclidean')
disana(dis.euc)
euc.pco <- pco(dis.euc,k=10)
barplot(euc.pco$eig)
plot(euc.pco)  # not so pretty

#Manhatan distance  #Better
dis.mht <- dist(rte3,"manhattan")
mht.pco <- pco(dis.mht,k=10)
plot(mht.pco)
points(mht.pco,rte3$berrep>0)

#binary
dis.bin <- dist(rte3,"binary")
bin.pco <- pco(dis.bin,k=10)
plot(bin.pco)  #don't use

#compare the plots
plot(euc.pco$eig/sum(euc.pco$eig),type="b",xlab="Axis Number",
     ylab="Fraction of Sum")
lines(mht.pco$eig/sum(mht.pco$eig),type="b",col=2)
lines(bin.pco$eig/sum(bin.pcoe$eig),type="b",col=3)
text(8.0,0.45,'Euclidean')
text(8.0,0.4,"Manhattan",col=2)
text(8.0,0.35,"binary",col=3)

### Show separate sites on the plot ###
C=which(rte2$Site=="C")  
N=which(rte2$Site=="N")

#Euclidean
pC=euc.pco$points[-N,]
pN=euc.pco$points[-C,]
ordiplot(euc.pco, type= "n", display="sites")  
points(pC, col="blue")
points(pN, col='green',pch=18)  

#manhattan
pCm=mht.pco$points[-N,]
pNm=mht.pco$points[-C,]
plot(mht.pco, type= "n", display="sites", title="PCoA_Manhattan")  
points(pCm, col="blue", pch=16)
points(pNm, col='green',pch=15)  
legend('topright', c("Offshore", "Nearshore"), bty="n", pch=c(16,15), col=c("blue", "green"), cex=.9)

######################################################################
#canberra
dis.can <- dist(rte3,"canberra")
can.pco <- pco(dis.can,k=10)
plot(can.pco)

pCc=can.pco$points[-N,]
pNc=can.pco$points[-C,]
plot(can.pco, type= "n", display="sites", title='PCoA_Canberra')  
points(pCc, col="blue",pch=16)
points(pNc, col='green',pch=15)  
legend('topright', c("Offshore", "Nearshore"), bty="n", pch=c(16,15), col=c("blue", "green"), cex=.9)




#bray
dis.bray<-vegdist(rte3, method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm = FALSE)
bray.pco <- pco(dis.bray,k=10)
plot(bray.pco)

pCb=bray.pco$points[-N,]
pNb=bray.pco$points[-C,]
plot(bray.pco, type= "n", display="sites", title='PCoA_Bray')  
points(pCb, col="blue",pch=16)
points(pNb, col='green',pch=15)  
legend('topright', c("Offshore", "Nearshore"), bty="n", pch=c(16,15), col=c("blue", "green"), cex=.9)



____________________________________________________________________________________

#NMDS

rte3.mds <- metaMDS(rte3, dist = "bray") 
rte3.mds
ordiplot(rte3.mds)
ordiplot(rte3.mds, type= "text", display="sites")  

C=which(rte2$Site=="C")  
N=which(rte2$Site=="N")
pC=rte3.mds$points[-N,]
pN=rte3.mds$points[-C,]
ordiplot(rte3.mds, type= "n", display="sites")  
points(pC, pch=16, col='blue')
points(pN, pch=18, col='green')  
legend('topright', c("CW", "SiteN"), bty="n", pch=c(16,18),col=c('blue','green'), cex=.9)
text(rte3.mds, display = "spec", cex=0.7, col="red")
###################################################
# Anosim #

ano<-anosim(rte3, rte2$Site, permutations = 999, distance = "bray", strata = NULL)
ano
#Call:
#anosim(dat = rte3, grouping = rte2$Site, permutations = 999,      distance = "bray") 
#Dissimilarity: bray 

#ANOSIM statistic R: 0.2052 
#Significance: 0.001 

#Permutation: free
#Number of permutations: 999


ano.eu<-anosim(rte3, rte2$Site, permutations = 999, distance = "euclidean", strata = NULL)
ano.eu

#ANOSIM statistic R: 0.1273 
#Significance: 0.001 

#Permutation: free
#Number of permutations: 999


############################################################################  
####Mantel Test###
#morphological distance#

dis.euc

#genetic distance#
mantel(f31.dist, env.dist)
 


############################################################################  