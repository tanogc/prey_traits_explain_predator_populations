##########################################################################################################
#                                                                                                         
# Supplementary Data 1                                                                                    
#                                                                                                         
# R script to quantify prey biodiversity metrics and prepare response variables (predator populations)                                                         
#                                                                                                         
# Populations of high-value predators reflect the traits of their prey                                          
#                                                                                                               
# Gutiérrez-Cánovas, C., Worthington, T.A., Jâms, I.B., Noble, D.G., Perkins, D.M., 
# Vaughan, I.P., Woodward, G., Ormerod, S.J. & I. Durance                                                                      
#                                                                                                               
# Code written by Cayetano Gutierrez-Canovas, email: cayeguti@um.es                                                                                                                                                            
##########################################################################################################

# working directory
setwd("your_folder")

## Loading libraries
library(FD)
library(vegan)
library(ade4)
library(sqldf)
library(reshape2)
library(ape)
library(plyr)

## Loading additional functions
source("0_quality_funct_space_fromdist.R")
source("0_FD_functions.R")

## Loading invertebrate abundance data
inv<-read.table("inv_ab.txt",h=T,sep="\t") # Loading invert data
e.gr<-read.table("fun_groups.txt",h=T,sep="\t") # Loading effect groups

## Loading environmental data
env_f<-read.table("env_fish.txt",h=T,sep="\t") # fish
env_b<-read.table("env_bird.txt",h=T,sep="\t") # bird

q<-sapply(env_f,class)=="numeric" | sapply(env_f,class)=="integer"# selecting quantitative variables
par(mfrow=c(3,3))
for (i in which(q==T)) hist(env_f[,i], main=names(env_f)[i])
par(mfrow=c(1,1))

env_f<-data.frame(site_code=env_f$site_code, alt=env_f$ALT25M, ph=env_f$pH)
env_b<-data.frame(site_code=env_b$site_code, alt=env_b$ALT25M, ph=env_b$pH)

## Loading fish and bird data
fish<-read.table("fish_dat.txt",h=T,sep="\t") # Loading fish data
bird<-read.table("bird_dat.txt",h=T,sep="\t") # Loading bird data
fish_mean<-read.table("fish_mean_dat.txt",h=T,sep="\t") # Loading fish summed data

fish<-fish[,-c(6:9)] # removing wet weights

## Matching invertebrate and predator data (fish and birds). Matching sites were used to select the fish samples.

# Bird data
bird_dat <- sqldf("select * from inv,bird where inv.site_code = bird.site_code")

# Selecting those samples with bird surveys
bird_dat<-bird_dat[which(bird_dat$birds==1),]
bird_env<-env_b[which(bird_dat$birds==1),]
any(bird_dat$site_code==bird_env$site_code)==F

# Fish data
fish_dat <- sqldf("select * from fish,inv where fish.site_code = inv.site_code") # selecting all the samples
fish_env <- sqldf("select * from env_f,inv where env_f.site_code = inv.site_code")

# Selecting those samples with a unique fish survey
trout_dat<-fish_dat[which(fish_dat$trout==1),]
salmon_dat<-fish_dat[which(fish_dat$salmon==1),]
trout_env<-fish_env[which(fish_env$trout==1),]
salmon_env<-fish_env[which(fish_env$salmon==1),]

## Preparing invertebrate data

# Selecting the invertebrate abundance samples
trout_inv_ab<-trout_dat[,17:105]
salmon_inv_ab<-salmon_dat[,17:105]
bird_inv_ab<-bird_dat[,8:96]

# Preparing fish and bird data
wagtail<-bird_dat$grey_wagtail
dipper<-bird_dat$dipper
trout<-trout_dat$trout_0+trout_dat$trout_1
salmon<-salmon_dat$salmon_0+salmon_dat$salmon_1

# Loading trait data
e_tr<-read.table("inv_e_traits_all.txt",h=T,sep="\t") # Loading taxa x effect trait matrix 
r_tr<-read.table("inv_r_tr_all.txt",h=T,sep="\t") # Loading taxa x effect trait matrix 

# arrangements
species<-e_tr$genus # storing taxon names
rownames(e_tr)<-species
e.tr<-e_tr[,7:ncol(e_tr)]
r.tr<-r_tr[,6:(ncol(r_tr)-1)]
extra<-e_tr[,5:6]

# Effect traits: Preparing traits - converting into percentage values
traits.blo<-c(4,2,4,4,5,5,7,9)
e.tr<-prep.fuzzy(e.tr[,-c(7:11,25:28,34:43)],traits.blo) # transform fuzzy codes to percentages
colSums(e.tr,na.rm=T)>0 -> sel.t # removing categories with sum=0
any(colSums(e.tr,na.rm=T)==0)

# Response traits: Preparing traits - converting into percentage values
traits.blo<-c(2,3,4,8,4,5,4)
which(rowSums(r.tr,na.rm=T)!=0)->no_r_tr
r.tr<-prep.fuzzy(r.tr[no_r_tr,],traits.blo) # transform fuzzy codes to percentages
rownames(r.tr)<-species[no_r_tr]

# Gower dissimilarity matrices

e.tr.ktab <- ktab.list.df(list(extra[1],extra[2],e.tr))
e.dist <- dist.ktab(e.tr.ktab, type= c("Q","N","F"))

r.tr.ktab <- ktab.list.df(list(r.tr))
r.dist <- dist.ktab(r.tr.ktab, type= c("F"))

## Loading invertebrate trait data
e.pco<-read.table("e_pco.txt",h=T,sep="\t") # Loading effect traits
r.pco<-read.table("r_pco.txt",h=T,sep="\t") # Loading response traits

## Subseting effect groups for each dataset
e.gr.t<-e.gr[intersect(rownames(e.pco), colnames(trout_inv_ab)),1]
e.gr.s<-e.gr[intersect(rownames(e.pco), colnames(salmon_inv_ab)),1]
e.gr.b<-e.gr[intersect(rownames(e.pco), colnames(bird_inv_ab)),1]

#############################
# Functional identity (FI)  #
#############################

par(mfrow = c(1,1))

## Functional group abundances (FG)

# Fish
calc.FR(trout_inv_ab,e.gr.t)->FG.ab.trout
calc.FR(salmon_inv_ab,e.gr.s)->FG.ab.salmon

# Birds

calc.FR(bird_inv_ab,e.gr.b)->FG.ab.bird

########################
# Functional diversity #
########################

# Effect trait metrics (FD-e)

# Functional richness of effect traits
fric_3d(trout_inv_ab,e.pco,8,prec="Qt")->FDe_t # trout
fric_3d(salmon_inv_ab,e.pco,8,prec="Qt")->FDe_s # salmon
fric_3d(bird_inv_ab,e.pco,8,prec="Qt")->FDe_b # birds

# Is there any NAs
any(is.na(FDe_t))
any(is.na(FDe_s))
any(is.na(FDe_b))
FDe_t[is.na(FDe_t)]<-0 # assigning 0 instead of NA
FDe_s[is.na(FDe_s)]<-0 # assigning 0 instead of NA
FDe_b[is.na(FDe_b)]<-0 # assigning 0 instead of NA

# Response trait diversity (FD-r)

# Functional richness of effect traits
fric_3d(trout_inv_ab,r.pco,7,prec="Qt")->FDr_t # trout
fric_3d(salmon_inv_ab,r.pco,7,prec="Qt")->FDr_s # salmon
fric_3d(bird_inv_ab,r.pco,7,prec="Qt")->FDr_b # birds

# Is there any NAs
any(is.na(FDr_t))
any(is.na(FDr_s))
any(is.na(FDr_b))
FDr_t[is.na(FDr_t)]<-0 # assigning 0 instead of NA
FDr_s[is.na(FDr_s)]<-0 # assigning 0 instead of NA
FDr_b[is.na(FDr_b)]<-0 # assigning 0 instead of NA

#################################################################
# Storing variables   
#################################################################
# Defining datasets

# Bird dataset
bird.set<-data.frame(sites=as.factor(bird_dat$site_code),
                      year=as.factor(bird_dat$year), 
                      wagtail,
                      dipper,
                      FG2=log(FG.ab.bird$ab.fgrs$FR2+1),
                      FG3=log(FG.ab.bird$ab.fgrs$FR3),
                      FDe=FDe_b,
                      FDr=log(FDr_b+0.01),
                      ric=specnumber(bird_inv_ab),
                      abun=log(rowSums(bird_inv_ab)),
                      bird_env[,-1])

# Trout dataset
trout.set<-data.frame(sites=as.factor(trout_dat$site_code),
                    season=as.factor(trout_dat$season),
                    year=as.factor(trout_dat$year), 
                    trout=trout,
                    FG2=log(FG.ab.trout$ab.fgrs$FR2+1),
                    FG3=log(FG.ab.trout$ab.fgrs$FR3),
                    FDe=FDe_t,
                    FDr=log(FDr_t+0.01),
                    ric=specnumber(trout_inv_ab),
                    abun=log(rowSums(trout_inv_ab)),
                    trout_env[,2:4])

# Salmon dataset
salmon.set<-data.frame(sites=as.factor(salmon_dat$site_code),
                      season=as.factor(salmon_dat$season),
                      year=as.factor(salmon_dat$year), 
                      salmon=salmon,
                      FG2=log(FG.ab.salmon$ab.fgrs$FR2+1),
                      FG3=log(FG.ab.salmon$ab.fgrs$FR3),
                      FDe=FDe_s,
                      FDr=log(FDr_s+0.01),
                      ric=specnumber(salmon_inv_ab),
                      abun=log(rowSums(salmon_inv_ab)),
                      salmon_env[,2:4])


# Preparing multiple-predator coincident dataset

# Salmonid mean biomass for each site
fish_mean$salmon0_mean+fish_mean$salmon1_mean->salmon.mean
fish_mean$trout0_mean+fish_mean$trout1_mean->trout.mean

# subsetting bird data to select only the sites with both birds and fish
dipper.all<-dipper[which(bird_dat$both==2)]
wagtail.all<-wagtail[which(bird_dat$both==2)]

trout.all<-fish_dat$trout_0[which(fish_dat$both==2)]+fish_dat$trout_1[which(fish_dat$both==2)]
salmon.all<-fish_dat$salmon_0[which(fish_dat$both==2)]+fish_dat$salmon_1[which(fish_dat$both==2)]

# Exporting coincident dataset (n=32)
write.table(data.frame(wagtail.all, trout.mean, salmon.mean, dipper.all),"coincident_set.txt",sep="\t",row.names=F)

# Multiple-predator index (mpi) dataset
mpi<-scale(wagtail.all)[,1]+scale(dipper.all)[,1]+scale(sqrt(trout.mean))[,1]+scale(sqrt(salmon.mean))[,1]

mpi<-log(mpi-min(mpi)+1)

# correlation between predator numbers
as.dist(round(cor(data.frame(wagtail.all,dipper.all,trout.mean,salmon.mean), method="spearman"),2))

# Averaging salmon dataset (coincident dataset)
ddply(salmon.set,.(sites),function(x) {colMeans(data.frame(x[,c(5:ncol(salmon.set))]))})->salmon.mean

mpi.set<-data.frame(mpi,salmon.mean)

# Exporting datasets

write.table(bird.set,"bird_set.txt",sep="\t",row.names=F)

write.table(trout.set,"trout_set.txt",sep="\t",row.names=F)

write.table(salmon.set,"salmon_set.txt",sep="\t",row.names=F)

write.table(mpi.set,"mpi_set.txt",sep="\t",row.names=F)

##########################################################################################




