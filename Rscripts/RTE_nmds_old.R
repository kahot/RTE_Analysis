rm(list=ls(all=TRUE))
require(vegan)
require(fields)
require(MASS)
setwd("~/Dropbox/LC-MS/NMDS")

### NN vs NC abacus dataset ###
nnnc<-read.table("NNvsNCdata.csv",sep=",", header = T)
#18  obs. or 5221 variables
nnnc1=apply(nnnc,1,as.numeric)  #This transform the matrix. Can be simply t(nnnc)

#transform the data according the Brook's file
nnnc2<-(nnnc1+1)
nnnc3<-log(nnnc2)
nnnc_nmds1<-metaMDS(nnnc2)
nnnc_nmds1   # Stress= 5.739475e-05
ordiplot(nnnc_nmds1, type= "text", display="sites")  
# 14 seems to be off, so eliminate the #14 (14th row, K17b)
nnnc4<-nnnc3[-14,]
nnnc_nmds2<-metaMDS(nnnc4)
nnnc_nmds2 #stress = 0.1207502 
ordiplot(nnnc_nmds2, type= "text", display="sites")  #better



###### All RTE Samples  ####
abacus <- read.csv("2016_Aug_01_abacus_sample1-39.csv")
abacus<-abacus[,4:125]
rte<-abacus_4_NMDS[,grepl('K10a|K10b|K10c|K11a|K11b|K11c|K12a|K12b|K12c|K13a|K13b|K13c|K14a|K14b|K14c|K3a|K3b|K3c|K4a|K4b|K4c|K5a|K5b|K5c|K6a|K6b|K6c', names(abacus_4_NMDS),)]


rte<-read.csv('RTE.csv', header=T)
rte1=apply(rte,1,as.numeric)
rte_nmds1<-metaMDS(rte1)
rte_nmds1  # stress = 0.1533669
ordiplot(rte_nmds1, type= "text", display="sites")  

rte2<-log(rte1+1)
rte2.2<-data.trans(rte2, method='log', plot=F)
rte_nmds2<-metaMDS(rte2.2, autotransform=F)
rte_nmds2 # stress = 0.1552547
ordiplot(rte_nmds2, type= "text", display="sites")  

rte3<-log(rte1+1, base=10)
rte_nmds3<-metaMDS(rte3, autotransform=F)
ordiplot(rte_nmds3, type= "text", display="sites")  

## Normalization  #don't use
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
range01(rte1)->rte_n
rte_nmds4<-metaMDS(rte_n, autotransform=F)
ordiplot(rte_nmds4, type= "text", display="sites")  
##############################

NN<-c(1,2,3,4,5,6,7,8,9)
CN<-c(10,11,12,13,14,15,16,17,18)
NC<-c(19,20,21,22,23,24,25,26)
CC<-c(27,28,29,30,31,32,33,34,35)

nn=rte_nmds1$points[NN,]	
cn=rte_nmds1$points[CN,]	
nc=rte_nmds1$points[NC,]
cc=rte_nmds1$points[CC,]

#blue=rgb('#0000FF')


ordiplot(rte_nmds1, type= "n", display="sites", xlim=c(-0.35, 0.35), ylim = c)  

ordiplot(rte_nmds1, type= "n", display="sites", xlim=c(-0.35, 0.35), ylim=c(-0.25, 0.4),xlab="NMDS Axis1", ylab="NMDS Axix2")  
points(nn, pch=21, col='green', bg='green', cex=2)
points(cn, pch=21, col='#0000FF', bg='#00FF00',cex=2,lwd=2.5)
points(nc, pch=21, col='green',bg='blue',cex=2,lwd=2.5)
points(cc, pch=21, col='blue', bg='blue',cex=2)



nn=rte_nmds2$points[NN,]	
cn=rte_nmds2$points[CN,]	
nc=rte_nmds2$points[NC,]
cc=rte_nmds2$points[CC,]

ordiplot(rte_nmds2, type= "n", display="sites")  
points(nn, pch=21, col='green', bg='green', cex=2)
points(cn, pch=21, col='#0000FF', bg='#00FF00',cex=2, lwd=2.5)
points(nc, pch=21, col='green',bg='blue',cex=2, lwd=2.5)
points(cc, pch=21, col='blue', bg='blue',cex=2)






#############################################
#NMDS 
#pcb.dist=vegdist(pcb)	

nnnc_nmds1<-metaMDS(nnnc1)
nnnc_nmds1   # Stress=0.02183945 
#Square root transformation
#Wisconsin double standardization

ordiplot(nnnc_nmds1, type= "text", display="sites")  

#plot with treatment info
#PCB treatment
pcbexposed<-c(1,2,3)
control<-c(4,5,6)

t1=pcb_nmds1$points[-control,]	
t2=pcb_nmds1$points[-pcbexposed,]	

ordiplot(pcb_nmds1, type= "n", display="sites")  
points(t1, pch=16, col='red')
points(t2, pch=16, col='blue')  

#############################################
### transformaiton --didn't change the output much
pcb_h=decostand(pcb1,method="hellinger")

pcb_nmds2<-metaMDS(pcb_h)
pcb_nmds2   # Stress=0 
ordiplot(pcb_nmds2, type= "text", display="sites") 









#############################################
### normalization
#############################################

#pcb_n<- scale(pcb1)  #-1 to 1

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
range01(pcb1)->pcb_s

pcb_nmds3<-metaMDS(pcb_s)
pcb_nmds3

ordiplot(pcb_nmds3, type="t", display="sites")	
t1=pcb_nmds3$points[-control,]	
t2=pcb_nmds3$points[-pcbexposed,]	
points(t1, pch=16, col='red')
points(t2, pch=16, col='blue')  



#calculate analysis of similarity statistics 
x<-c(1,1,1,2,2,2)
ano=anosim(pcb, x,permutations=1000)
ano

#ANOSIM statistic R: 0.5926 
#Significance: 0.1 
#Permutation: free
#Number of permutations: 719



#determine singletons -proteins only showed up once in spectral counting among 6 samples
Singletons=pcb[,colSums(pcb)<2] #388 proteins


##################################################################################################
##################################################################################################
#using abacus data file

pcbr<-read.table("pcb_rawd.csv",sep=",")
#6  obs. or 3086 variables
pcbr1<-pcbr[2:5222,]
pcbr2=apply(pcbr1,1,as.numeric)


#############################################
#NMDS 
#pcb.dist=vegdist(pcb)	

pcbr_nmds1<-metaMDS(pcbr2)
pcbr_nmds1   #  0.06213517
#Square root transformation
#Wisconsin double standardization

ordiplot(pcbr_nmds1, type= "text", display="sites")  

#plot with treatment info
#PCB treatment
pcbexposed<-c(1,2,3)
control<-c(4,5,6)

t1=pcb_nmds1$points[-control,]	
t2=pcb_nmds1$points[-pcbexposed,]	

ordiplot(pcb_nmds1, type= "n", display="sites")  
points(t1, pch=16, col='red')
points(t2, pch=16, col='blue')  



pcbr_h=decostand(pcbr2,method="hellinger")

pcbr_nmds2<-metaMDS(pcbr_h)
pcbr_nmds2   # Stress=0 
ordiplot(pcbr_nmds2, type= "text", display="sites") 


############################################
### normalization
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
range01(pcbr2)->pcbr_s

pcbr_nmds2<-metaMDS(pcbr_s)
pcbr_nmds2

ordiplot(pcbr_nmds2, type="t", display="sites")	
t1=pcbr_nmds2$points[-control,]	
t2=pcbr_nmds2$points[-pcbexposed,]	
points(t1, pch=16, col='red')
points(t2, pch=16, col='blue')  


