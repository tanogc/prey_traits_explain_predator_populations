###################################################################################################
#                                                                                                               
# Supplementary Data 1                                                                                          
#                                                                                                               
# R script to classify invertebrate prey into functional groups based on effect traits                        
#                                                                                                               
# Populations of high-value predators reflect the traits of their prey                                          
#                                                                                                               
# Gutiérrez-Cánovas, C., Worthington, T.A., Jâms, I.B., Noble, D.G., Perkins, D.M., 
# Vaughan, I.P., Woodward, G., Ormerod, S.J. & I. Durance                                                                      
#                                                                                                               
# Code written by Cayetano Gutierrez-Canovas, email: cayeguti@um.es                                                                                                                                                            
###################################################################################################

# working directory
setwd("your_folder")
plot_folder<-paste(getwd(),"/plots/",sep="")

## Loading libraries
library(FD)
library(ecodist)
library(vegan)
library(ade4)
library(ape)
library(sqldf)
library(plyr)

## Loading additional functions
source("0_quality_funct_space_fromdist.R")

# Loading trait data
e_tr<-read.table("e_traits_all.txt",h=T,sep="\t") # Loading taxa x effect trait matrix 
r_tr<-read.table("r_traits_all.txt",h=T,sep="\t") # Loading taxa x effect trait matrix 

# arrangements
species<-e_tr$genus # storing taxon names
rownames(e_tr)<-species
e.tr<-e_tr[,7:ncol(e_tr)]
r.tr<-r_tr[,6:(ncol(r_tr)-1)]
extra<-e_tr[,5:6]

# Effect traits: Preparing traits - converting into percentage values
#traits.blo<-c(4,2,5,4,4,5,4,5,2,5,3,7,9)
#e.tr<-prep.fuzzy(e.tr,traits.blo) # transform fuzzy codes to percentages

# Reduced version
traits.blo<-c(4,2,4,4,5,5,7,9)
e.tr<-prep.fuzzy(e.tr[,-c(7:11,25:28,34:43)],traits.blo) # transform fuzzy codes to percentages

colSums(e.tr,na.rm=T)>0 -> sel.t # removing categories with sum=0
any(colSums(e.tr,na.rm=T)==0)

# Response traits: Preparing traits - converting into percentage values
traits.blo<-c(2,3,4,8,4,5,4)
which(rowSums(r.tr,na.rm=T)!=0)->no_r_tr
r.tr<-prep.fuzzy(r.tr[no_r_tr,],traits.blo) # transform fuzzy codes to percentages
rownames(r.tr)<-species[no_r_tr]
  
names(e.tr)
names(r.tr)

# Gower dissimilarity matrices

e.tr.ktab <- ktab.list.df(list(extra[1],extra[2],e.tr))
e.dist <- dist.ktab(e.tr.ktab, type= c("Q","N","F"))

r.tr.ktab <- ktab.list.df(list(r.tr))
r.dist <- dist.ktab(r.tr.ktab, type= c("F"))

# Effect trait functional space
# Estimating the optimum number of dimensions. Also the functional space to be used in fric 
qual_fs<-quality_funct_space_fromdist(e.dist,  nbdim=10)
qual_fs$meanSD<0.01 # 8D seems to be an appropiate number of dimensions
qual_fs$meanSD

# Storing effect trait functional space
write.table(qual_fs$fpc,"fpc_all.txt",sep="\t")

# Classifying functional groups
# checking for negative eigenvalues
dudi.pco(e.dist,scannf = F,nf=8)->e.pco
length(which(e.pco$eig<0)) # No negative eigenvalues licenced us to use Ward clustering method

# Saving the response trait space
write.table(e.pco$li,"e_pco.txt",sep="\t")

# Total and cummulative variance explained by each axes
e.pco$eig[1:8]/sum(e.pco$eig)
sum(e.pco$eig[1:8]/sum(e.pco$eig))

# Storing two-first axes
e.pco$li[,1]->axis_1
e.pco$li[,2]->axis_2
names(axis_2)<-names(axis_1)<-rownames(e.pco$li)

# Which taxa is at the axis extremes
round(sort(axis_1),2)
round(sort(axis_2),2)

# Spearman rank correlation between effect traits and effect trait-space two first axes
round(cor(data.frame(drift_index=extra[,1],e.tr),e.pco$li[,1:8],method="spearman",
          use="pairwise.complete.obs"),2)->fpc.axes
write.table(fpc.axes,"fpc_axes.txt",sep="\t") # Exporting correlations

# plot dengrogram of species based on effect traits
e.clust <- hclust(e.dist, method = "ward.D2")

cor(cophenetic(e.clust),e.dist)

# Number of groups and species assignation to groups
cut.g <- 4
e.gr <- cutree(e.clust, k = cut.g)

# Taxa included in each group
rownames(e.tr[which(e.gr==1),])
rownames(e.tr[which(e.gr==2),])
rownames(e.tr[which(e.gr==3),])
rownames(e.tr[which(e.gr==4),])

# Assigning colours for groups
fgr.col<-c("blue","dark green","red","violet")

# plot basic tree (10 x 8.5)
pdf(file=paste(plot_folder,"func_cluster.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=8.5,height=17)
plot(as.phylo(e.clust), cex = 0.8, label.offset= 0.05, tip.color = fgr.col[e.gr])
dev.off()

# Exporting functional groups
write.table(as.data.frame(e.gr),"fun_groups.txt",sep="\t")

# How groups spread in the functional space
pdf(file=paste(plot_folder,"func_space_1vs2.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=7.5,height=7.5)
s.class(e.pco$li[,1:2],fac = as.factor(e.gr), col=fgr.col,clabel=1.75, cpoint=2.25,pch=c(17,16,15,18)[e.gr], cellipse = 0)
dev.off()

pdf(file=paste(plot_folder,"func_space_1vs3.pdf",sep=""),onefile=T,width=7.5,height=7.5)
s.class(e.pco$li[,c(1,3)],fac = as.factor(e.gr), col=fgr.col,clabel=1.75, cpoint=2.25)
dev.off()

# Mean effect trait values for each group

e.tr.q<-data.frame(e.gr,drift_index=extra$drift_index,e.tr)

ddply(e.tr.q,.(e.gr),function(x) colMeans(x,na.rm = T))->av_tr

t(av_tr[,-1])->av_tr
colnames(av_tr)<-paste("g",1:4,sep="")

high_cal<-rep(NA,cut.g)

table(e.tr[which(e.gr==1),2])[1]/length(e.tr[which(e.gr==1),2])->high_cal[1]
table(e.tr[which(e.gr==2),2])[1]/length(e.tr[which(e.gr==2),2])->high_cal[2]
table(e.tr[which(e.gr==3),2])[1]/length(e.tr[which(e.gr==3),2])->high_cal[3]
table(e.tr[which(e.gr==4),2])[1]/length(e.tr[which(e.gr==4),2])->high_cal[4]

av_tr[-1,]<-100*av_tr[-1,]        
av_tr<-round(av_tr,0)
av_tr<-rbind(av_tr,round(high_cal*100,0))
rownames(av_tr)[nrow(av_tr+1)]<-"high_calcium"

write.table(av_tr,"FG_traits.txt",sep="\t")

# Response traits

# Response trait functional space
# Estimating the optimum number of dimensions. Also the functional space to be used in fric 

qual_fs<-quality_funct_space_fromdist(r.dist, nbdim=10)
qual_fs$meanSD<0.01 # 7D seems to be an appropiate number of dimensions
qual_fs$meanSD

# checking for negative eigenvalues
dudi.pco(r.dist,scannf = F,nf=7)->r.pco
length(which(r.pco$eig<0)) # No negative eigenvalues licenced us to use Ward clustering method

# Saving the response trait space
write.table(r.pco$li,"r_pco.txt",sep="\t")

# Variance explained by the two first axes
round(r.pco$eig[1:7]/sum(r.pco$eig),2)
cumsum(round(r.pco$eig[1:7]/sum(r.pco$eig),3))

# Two first axes
r.pco$li[,1]->axis_1
r.pco$li[,2]->axis_2
names(axis_2)<-names(axis_1)<-rownames(r.pco$li)

# Which taxa is at the axis extremes
sort(axis_1)
sort(axis_2)

# Spearman rank correlation between effect traits and effect trait-space two first axes
round(cor(r.tr,r.pco$li[,1:7],method="spearman",use="pairwise.complete.obs"),2)->fpc.r.axes
fpc.r.axes

# Averaging response traits for each effect group
r.tr.q<-data.frame(e.gr[no_r_tr],r.tr)
ddply(r.tr.q,.(e.gr),function(x) colMeans(x,na.rm = T))->av_tr

# Transposing matrix
round(t(av_tr[,-c(1:2)])*100,0)->av_tr
colnames(av_tr)<-paste("g",1:4,sep="")

# Saving results
write.table(av_tr,"av_res_traits.txt",sep="\t")

