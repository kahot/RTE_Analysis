#Explore the RTE data with metaMDS

rm(list=ls(all=TRUE))
require(vegan)
require(fields)
require(MASS)
require(fossil)
require(mgcv)
setwd("~/Dropbox/LC-MS/NMDS")


###### All RTE Samples with techinical replicates  ####

#read in the file
rte1<-read.csv('RTE.csv', header=T)
rte1=apply(rte,1,as.numeric)
rte_nmds1<-metaMDS(rte1)
#Square root transformation
#Wisconsin double standardization
rte_nmds1  

#global Multidimensional Scaling using monoMDS

#Data:     wisconsin(sqrt(rte1)) 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.07657125 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘wisconsin(sqrt(rte1))’ 


ordiplot(rte_nmds1, type= "text", display="sites")  

NN<-c(1,2,3,4,5,6,7,8,9)
CN<-c(10,11,12,13,14,15,16,17,18)
NC<-c(19,20,21,22,23,24,25,26)
CC<-c(27,28,29,30,31,32,33,34,35)

nn=rte_nmds1$points[NN,]	
cn=rte_nmds1$points[CN,]	
nc=rte_nmds1$points[NC,]
cc=rte_nmds1$points[CC,]

ordiplot(rte_nmds1, type= "n", display="sites")  

ordiplot(rte_nmds1, type= "n", display="sites", xlim=c(-0.4, 0.35), ylim=c(-0.25, 0.4),xlab="NMDS Axis1", ylab="NMDS Axix2")  
points(nn, pch=21, col='#00FF00', bg="#00FF00CC", cex =3)
points(cn, pch=21, col='#0000FF', bg='#00FF00',cex=3,lwd=3)
points(nc, pch=21, col='green',bg='blue',cex=3,lwd=3)
points(cc, pch=21, col='blue', bg='#0000FF99',cex=3)


# specify Manhattan distance
rte_nmds1m<-metaMDS(rte1,dist = "manhattan")
rte_nmds1m
#Stress:  0.07657125 
ordiplot(rte_nmds1m, type= "text", display="sites")  

#compare Manhattan distance with Bray 
stressplot(rte_nmds1)
stressplot(rte_nmds1m)
plot(procrustes(rte_nmds1,rte_nmds1m))
#almost the same

#log+1 transformation
rte2<-log(rte1+1) 
#rte2.2<-data.trans(rte2, method='log', plot=F)
rte_nmds2<-metaMDS(rte2, autotransform=F)
rte_nmds2 
#Data:     rte2 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.06611574 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘rte2’

ordiplot(rte_nmds2, type= "text", display="sites")  




###################################################################
#####  Using 12 biological replicates (tehc rep pooled)
rte_bio<-read.csv('RTE_BioRep12.csv', header=T)
rte_bio1=t(rte_bio)
rte_bio2<-rte_bio1[-c(1:2),]
rte_bio3<-apply(rte_bio2,2,as.numeric)
rtebio1<-metaMDS(rte_bio3)
ordiplot(rtebio1,type= "text", display="sites" )

rtebio1
#global Multidimensional Scaling using monoMDS

#Data:     wisconsin(sqrt(rte_bio3)) 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.09641986 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘wisconsin(sqrt(rte_bio3))’ 

score_rtebio1<-scores(rtebio1, display = "species") 
write.table(score_rtebio1,file="scores_rtebio1_NMDS.txt")
head(score_rtebio1)

treatment<-c("NN","NN","NN","ON","ON","ON","NO","NO","NO","OO","OO","OO")
#rtedata<-cbind(treatment,rte_bio2)

# 
plot(rtebio1)
ordipointlabel(rtebio1, display='species',add = TRUE, col = "forestgreen")
#Error in box[j, ] : subscript out of bounds


NN<-c(1,2,3)
CN<-c(4,5,6)
NC<-c(7,8,9)
CC<-c(10,11,12)

nn=rtebio1$points[NN,]	
cn=rtebio1$points[CN,]	
nc=rtebio1$points[NC,]
cc=rtebio1$points[CC,]

#blue=rgb('#0000FF')


ordiplot(rtebio1, type= "n", display="sites")  

ordiplot(rtebio1, type= "n", display="sites", xlim=c(-0.35, 0.35), ylim=c(-0.25, 0.3),xlab="NMDS Axis1", ylab="NMDS Axix2")  
points(nn, pch=21, col='green', bg='green', cex=1.5)
points(cn, pch=21, col='#0000FF', bg='#00FF00',cex=1.5,lwd=2)
points(nc, pch=21, col='green',bg='blue',cex=1.5,lwd=2.5)
points(cc, pch=21, col='blue', bg='blue',cex=1.5)


### run ANOSIM  ####

group<-c('N','N','N','N','N','N','O','O','O','O','O','O')
ano=anosim(rte_bio3,group,permutations = 5000)
ano
#Call:
#  anosim(dat = rte_bio3, grouping = group) 
#Dissimilarity: bray 

#ANOSIM statistic R: 0.7074 
#Significance: 0.0031994 

#Permutation: free
#Number of permutations: 5000


ano2=anosim(rte_bio3,treatment, permutations = 5000)
ano2

#Call:
#  anosim(dat = rte_bio3, grouping = treatment) 
#Dissimilarity: bray 

#ANOSIM statistic R: 0.5864 
#Significance: 0.00019996 

#Permutation: free
#Number of permutations: 5000


# Compare OO vs NO and NN vs ON
ncorals<-rte_bio3[1:6,]
ocorals<-rte_bio3[7:12,]
n_treat=c('NN','NN','NN','ON',"ON","ON")
o_treat=c("NO","NO","NO","OO","OO","OO")

anoN=anosim(ncorals,n_treat)
anoN
#ANOSIM statistic R: -0.07407 
#Significance: 0.7 
#Permutation: free
#Number of permutations: 719

anoO=anosim(ocorals,o_treat)
anoO
#ANOSIM statistic R: 0.5556 
#Significance: 0.1 

nn.corals1<-rte_bio3[1:3,]
on.corals1<-rte_bio3[4:6,]
no.corals1<-rte_bio3[7:9,]
oo.corals1<-rte_bio3[10:12,]
n.corals<-rbind(nn.corals1, no.corals1)
o.corals<-rbind(on.corals1,oo.corals1)
siten<-c('NN','NN','NN',"NO","NO","NO")
siteo<-c('ON','ON','ON',"OO","OO","OO")

#Permutation: free
#Number of permutations: 719
#ANOSIM statistic R: 0.5556 
#Significance: 0.1 

anosim(n.corals,siteo)
#nosim(dat = n.corals, grouping = siteo) 
#Dissimilarity: bray 
#ANOSIM statistic R: 0.5556 
#Significance: 0.1 
#Permutation: free
#Number of permutations: 719

anosim(o.corals,siteo)
#anosim(dat = o.corals, grouping = siteo) 
#Dissimilarity: bray 

#ANOSIM statistic R: 1
#Significance: 0.1 

#Permutation: free
#Number of permutations: 719








## try PCA  # doesn't look right

pca1<-prcomp(log(rte_bio3+1))
biplot(pca1)
rda(log(rte_bio3+1))->pca_rte1
plot(pca_rte1, type='t')
pca_rte1
#              Inertia Rank
#Total            1166     
#Unconstrained    1166   11
#Inertia is variance 

#Eigenvalues for unconstrained axes:
#  PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
#313.37 178.28 150.86 116.44 106.45  80.59  67.40  53.35  43.42  33.81  22.26 


biplot(pca_rte1)

##############################
## Normalization  #don't use!! 
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

#smaller points

ordiplot(rte_nmds1, type= "n", display="sites", xlim=c(-0.35, 0.35), ylim=c(-0.25, 0.4),xlab="NMDS Axis1", ylab="NMDS Axix2")  
points(nn, pch=21, col='green', bg='green', cex=1.5)
points(cn, pch=21, col='#0000FF', bg='#00FF00',cex=1.5,lwd=2)
points(nc, pch=21, col='green',bg='blue',cex=1.5,lwd=2.5)
points(cc, pch=21, col='blue', bg='blue',cex=1.5)


######################################################## 
### Find Variables that expalin the ordination axis 
##use rtebio1 ordination results


####@@@@@ we can skip this @@@@@###
#re-read the data file to have the column name to merge together
bio12<-read.csv('RTE_BioRep12.csv', header=F)
tbio12<-t(bio12)
tbio12_2<-tbio12[-c(1:2),]


#read Biol Process metadata file
BP1<-read.csv("BP_pvalue.csv", header=F, na.strings=c("","NA"))
tBP1<-t(BP1)
colnames(tBP1) = tBP1[1, ]
tBP1<-tBP1[-1,]
tBP2<-tBP1[,-1]

tBP3<-apply(tBP2,2,as.numeric)
#combine the two datafile
RTE_bp=merge(tbio12_2,tBP1, by.x=1, by.y=1)
#write.csv(RTE_bp, file="RTE_bp_combined.csv")
#extract the metada portion
bp<-RTE_bp[,4065:4185]
bp<-apply(bp,2,as.numeric)
#write.csv(bp,file="bp_metadata.csv")
data.frame(tBP3)->tBP4
####@@@@@ we can skip this till here @@@@@###




##%%%%  1. Using Biological Process P-values 
### 1.1 replaced non-significant P-value  and missing cells with 1
#data formatting
biopro1<-read.csv("BP_pvalue_replaced.csv", header=F)
biopro1<-t(biopro1)
colnames(biopro1) = biopro1[1, ]
biopro1<-biopro1[-1,]
biopro1.1<-biopro1[,-1]

biopro1.1<-apply(biopro1.1,2,as.numeric)
biopro1.1<-data.frame(biopro1.1)

## selected BP (not-replaced)
biopro1<-read.csv("BP_pvalue_replaced_selected1.csv", header=F)
biopro1<-t(biopro1)
colnames(biopro1) = biopro1[1, ]
biopro1<-biopro1[-1,]
biopro1.1<-biopro1[,-1]

biopro1.1<-apply(biopro1.1,2,as.numeric)
biopro1.1<-data.frame(biopro1.1)

#replaced the temrs with numbers
biopro1<-read.csv("BP_pvalue_replaced_selected3.2.csv", header=F)
biopro1<-t(biopro1)
colnames(biopro1) = biopro1[1, ]
biopro1<-biopro1[-1,]
biopro1.1<-biopro1[,-1]
biopro1.1<-apply(biopro1.1,2,as.numeric)
biopro1.1<-data.frame(biopro1.1)


## selected BP, ns P-values replaced with 1 (better to use non-replaced one)
#biopro1<-read.csv("BP_pvalue_replaced_selected2.csv", header=F)

#run ENVFIT  
en1<-envfit(rtebio1,biopro1.1, choices=1:2, w=weights(rtebio1), na.rm = TRUE)
en1

NN<-c(1,2,3)
CN<-c(4,5,6)
NC<-c(7,8,9)
CC<-c(10,11,12)

nn=rtebio1$points[NN,]	
cn=rtebio1$points[CN,]	
nc=rtebio1$points[NC,]
cc=rtebio1$points[CC,]

ordiplot(rtebio1, type= "n", display="sites")  

ordiplot(rtebio1, type= "n", display="sites", xlim=c(-0.35, 0.35), ylim=c(-0.25, 0.3),xlab="NMDS Axis1", ylab="NMDS Axix2")  
points(nn, pch=21, col='green', bg='green', cex=2)
points(cn, pch=21, col='#0000FF', bg='#00FF00',cex=2,lwd=2)
points(nc, pch=21, col='green',bg='blue',cex=2,lwd=2.5)
points(cc, pch=21, col='blue', bg='blue',cex=2)

plot(en1, cex=0.3)


####  1.3 run with P-balue BP data with blnak filled with 1 ###
bp.f<-read.csv("BP_pvalue_filled.csv", header=F)
bp.f<-t(bp.f)
colnames(bp.f) = bp.f[1, ]
bp.f<-bp.f[-1,]
bp.f.1<-bp.f[,-1]

bp.f.1<-apply(bp.f.1,2,as.numeric)
bp.f.1<-data.frame(bp.f.1)


en_filled<-envfit(rtebio1, bp.f.1, choices=1:2, w=weights(rtebio1), na.rm = TRUE)
en_filled

### ENVFIT with PCA results

en_filled_pca<-envfit(pca_rte1, bp.f.1, choices=1:2, w=weights(rtebio1), na.rm = TRUE)
en_filled_pca


## Selected 33 Biol. processes (+1 detox) above to re-run the envfit & RDA
bp.s<-read.csv("BP_pvalue_selected.csv", header=F)
bp.s<-t(bp.s)
colnames(bp.s) = bp.s[1, ]
bp.s<-bp.s[-1,]
bp.s.1<-bp.s[,-1]

bp.s.1<-apply(bp.s.1,2,as.numeric)
bp.s.1<-data.frame(bp.s.1)


en_selected<-envfit(rtebio1,bp.s.1, choices=1:2, w=weights(rtebio1), na.rm = TRUE)
en_selected

ordiplot(rtebio1, type= "n", display="sites")  
points(nn, pch=21, col='green', bg='green', cex=2.5)
points(cn, pch=21, col='#0000FF', bg='#00FF00',cex=2.5,lwd=2.5)
points(nc, pch=21, col='green',bg='blue',cex=2.5,lwd=2.5)
points(cc, pch=21, col='blue', bg='blue',cex=2.5)

plot(en_selected, cex=.6)



##run RDA using the signficant BP 
bio12<-read.csv('RTE_BioRep12.csv', header=F)
tbio12<-t(bio12)
tbio12_2<-tbio12[-c(1:2),]

RTE_aba<-tbio12_2[,-1]
RTE_aba<-apply(RTE_aba,2,as.numeric)
RTE_aba<-data.frame(RTE_aba)

RDA_bp1<-rda(RTE_aba ~ ., data=bp.s.1,na.action=na.omit )

anova(RDA_bp1,step=1000,perm.max=1000)
#Number of permutations: 999
#         Df Variance      F Pr(>F)   
#Model     3  1032702 3.9223  0.022 *
#Residual  8   702109                
                 

anova(RDA_bp1, by="term")   #, permu=1000
# See file RDA_results_RTE_BP_Pvalues_selected.txt





##%%%%  2. Using Biological Process number of poteins 
#numer of proteins
BP2<-read.csv("BP_noprotein_eliminated.csv", header=F,na.strings=c("","NA"))
BP2<-t(BP2)
colnames(BP2) = BP2[1, ]
BP2<-BP2[-1,]

BP2.1<-BP2[,-1]
BP2.1<-apply(BP2.1,2,as.numeric)
BP2.1<-data.frame(BP2.1)


#run envfit with no of proteins data
en2<-envfit(rtebio1,BP2.1, choices=1:2, w=weights(rtebio1), na.rm = TRUE)
en2

ordiplot(rtebio1, type= "n", display="sites")  
points(nn, pch=21, col='green', bg='green', cex=2.5)
points(cn, pch=21, col='#0000FF', bg='#00FF00',cex=2.5,lwd=2.5)
points(nc, pch=21, col='green',bg='blue',cex=2.5,lwd=2.5)
points(cc, pch=21, col='blue', bg='blue',cex=2.5)

plot(en2, cex=.6)



### ENVFIT with PCA 
en2_pca<-envfit(pca_rte1, BP2.1, choices=1:2, w=weights(rtebio1), na.rm = TRUE)
en2_pca


########################################################################################
#create the spectral count data for ENVFIT
#library(data.table)
# read the protein list file for Biological Processes
spec_list<-read.csv('spec_count_list1.csv', header=F)
#transpose 
list2<-t(spec_list)
list2.2<-data.frame(list2)
#write.table(list2, file='list2.csv')


#read the spectral count data
count<-read.csv("RTE_spec_count.csv", header = F)
#nnnc<-read.csv('NNNC.count.list.csv', header =F)

a=1
bp_summary<-list()

for (i in 1:121) {
  a=i
  bp<-list2.2[,a]
  bp_count<-count[count$V1 %in% bp, ]
  bp_sum<-apply(bp_count[,2:13],2,sum)
  bp_summary[[i]]<-bp_sum
}

summary<-do.call("rbind",bp_summary)
#write.table(summary,'BP_count_sum2.txt')
  
summary.t<-apply(summary,1,as.numeric)

en3<-envfit(rtebio1,summary.t, choices=1:2, w=weights(rtebio1), na.rm = TRUE)
en3

ordiplot(rtebio1, type= "n", display="sites")  
points(nn, pch=21, col='green', bg='green', cex=2.5)
points(cn, pch=21, col='#0000FF', bg='#00FF00',cex=2.5,lwd=2.5)
points(nc, pch=21, col='green',bg='blue',cex=2.5,lwd=2.5)
points(cc, pch=21, col='blue', bg='blue',cex=2.5)

plot(en3, cex=.6, col='blue')

#### with selected 25 BP

plot(en3, cex=.7, col='blue',p.max=0.1)  # this is the same plot as below with 3n4


summary_select<-summary.t[,c(15,16,28,29,30,31,33,34,35,36,38,39,45,46,60,67,68,90,91,92,95,96,107,108,109)]
#bp25<-read.csv("BP-spec.count-selected.csv", header=F,na.strings=c("","NA"))
#bp25<-t(bp25)
#colnames(bp25) = bp25[1, ]
#bp25<-bp25[-1,]
#bp25.1<-bp25[,-1]
#bp25.1<-apply(bp25.1,2,as.numeric)
#bp25.1<-data.frame(bp25.1)

# include p<0.1
summary_select<-summary.t[,c(9,15,16,27,28,29,30,31,32,33,34,35,36,37,38,39,45,46,55,58,59,60,67,68,69,70,90,91,92,95,96,103,107,108,109,119,120)]


en4<-envfit(rtebio1,summary_select, choices=1:2, w=weights(rtebio1), na.rm = TRUE)
en4

ordiplot(rtebio1, type= "n", display="sites",xlim=c(-0.35, 0.35), ylim=c(-0.25, 0.3),xlab="NMDS Axis1", ylab="NMDS Axix2")  
points(nn, pch=21, col='green', bg='green', cex=2.5)
points(cn, pch=21, col='#0000FF', bg='#00FF00',cex=2.5,lwd=2.5)
points(nc, pch=21, col='green',bg='blue',cex=2.5,lwd=2.5)
points(cc, pch=21, col='blue', bg='blue',cex=2.5)
#ordiellipse(rtebio1,treatment,col='darkgrey')

plot(en4, cex=.7)
ordiellipse(rtebio1,treatment, kind = "se", conf = 0.95, col='grey')
ordispider(rtebio1,treatment, col = "cyan", label= TRUE)

#with(summary_select, ordiellipse(rtebio1,treatment, kind = "se", conf = 0.95))
#with(summary_select, ordispider(rtebio1,treatment, col = "grey", label= TRUE))
#with(dune.env, ordihull(dune.ca, Management, col="blue", lty=2))


treatment





#RDA_bp1<-rda(RTE_aba ~ ., data=bp25.1,na.action=na.omit )
#anova(RDA_bp1,step=1000,perm.max=1000)



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


#########################################
### Run envfit with all original proteins #####
en_protein<-envfit(rtebio1,rte_bio3, choices=1:2, w=weights(rtebio1), na.rm = TRUE)
en_protein

envfit.data<-data.frame((en_protein$vectors)$arrows, (en_protein$vectors)$r, (en_protein$vectors)$pvals)

write.table(envfit.data,"envfit_protein_rte.txt")
