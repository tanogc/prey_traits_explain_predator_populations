###########################################################################################################
#                                                                                                          
# Supplementary Data 1                                                                                     
#                                                                                                          
# R script to perform the multi-threshold analysis                                                                                                                             #
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

#setwd("your_folder")
assum_folder<-paste(getwd(),"/plots/assum/",sep="")
res_folder<-paste(getwd(),"/plots/",sep="")

#library(devtools)
#install_github("jebyrnes/multifunc", force=T)

# Loading libraries
library(multifunc)
library(ggplot2)
library(gridExtra)
library(viridis)

###### Loading full datasets
dat<-read.table("coincident_set.txt",h=T,sep="\t")
mpi.set<-read.table("mpi_set.txt",h=T,sep="\t")

#####################
# Multiple threshold analysis
######################

fam_mod="gaussian"

######################
# ric (taxon richness)
######################

mpi_dat<-data.frame(wagtail=wagtail.all,
                   trout=dat$trout.mean,
                   salmon=(fish_mean$salmon0+fish_mean$salmon1),
                   dipper=dat$dipper.all,
                   ric=mpi.set$ric)

var_lab="ric"

es_vars<-qw(wagtail,trout,salmon,dipper)

pop_thresh<-getFuncsMaxed(mpi_dat, es_vars, threshmin=0.01, threshmax=0.99, prepend=c("ric"), maxN=7)

LinearSlopes<-getCoefTab(funcMaxed ~ ric, data=pop_thresh, coefVar="ric", family=fam_mod)

## Number of populations above the threshold 

pdf(file=paste(res_folder,"thres1_ric.pdf",sep=""),onefile=T,width=6,height=4.5)
par(mfrow=c(1,2))
pop_thresh$percent <- 100*pop_thresh$thresholds
ggplot(data=pop_thresh, aes(x=ric, y=funcMaxed, group=percent)) +
  ylab(expression("Number of ES")) +
  xlab(var_lab) +
  stat_smooth(method="glm", method.args = list(family = fam_mod), lwd=0.8, fill=NA, aes(color=percent)) +
  theme_bw(base_size=18) +
  scale_color_gradient(name="Threshold", low="blue", high="red") +
  ylim(0,4) +theme(axis.text=element_text(size=24),
                   axis.title=element_text(size=24,face="bold"))
dev.off()

# Values of the diversity slope at different levels of the threshold

pdf(file=paste(res_folder,"thres2_ric.pdf",sep=""),onefile=T,width=5.5,height=4.5)
popSlopes <- ggplot(LinearSlopes, aes(x=thresholds)) +
  geom_ribbon(fill="grey50", aes(x=thresholds*100, ymin=Estimate-1.96*LinearSlopes[["Std. Error"]], 
                                 ymax=Estimate+1.96*LinearSlopes[["Std. Error"]])) +
  geom_point(aes(x=thresholds*100, y=Estimate)) +
  ylab("Biodiversity effect on ES\n") +
  xlab("\nThreshold (%)") +
  geom_abline(intercept=0, slope=0, lwd=1, linetype=2) +
  theme_bw(base_size=14)

popSlopes+theme(axis.text=element_text(size=20),
               axis.title=element_text(size=20,face="bold"))

dev.off()

# Indices
IDX <- getIndices(LinearSlopes, pop_thresh, funcMaxed ~ ric, divvar = "ric")
IDX

write.table(IDX, "ric_idx.txt",sep="\t") # ric

## Minimum and maximum number of ES sustained

for (i in seq(0.01,0.99, by=0.01)) print(c(i, round((max(predict(glm(funcMaxed ~ ric, data=subset(pop_thresh, pop_thresh$thresholds==i), family=fam_mod)))),1)))

#########################
# abun (taxon abundance)
#########################

mpi_dat<-data.frame(wagtail=wagtail.all,
                   trout=dat$trout.mean,
                   salmon=(fish_mean$salmon0+fish_mean$salmon1),
                   dipper=dat$dipper.all,
                   abun=(mpi.set$abun))

var_lab="abun"

vars<-qw(wagtail,trout,salmon,dipper)

pop_thresh<-getFuncsMaxed(mpi_dat, vars, threshmin=0.01, threshmax=0.99, prepend=c("abun"), maxN=7)

LinearSlopes<-getCoefTab(funcMaxed ~ abun, data=pop_thresh, coefVar="abun", family=fam_mod)

## Number of populations above the threshold 

pdf(file=paste(res_folder,"thres1_abun.pdf",sep=""),onefile=T,width=6,height=4.5)
par(mfrow=c(1,2))
pop_thresh$percent <- 100*pop_thresh$thresholds
ggplot(data=pop_thresh, aes(x=abun, y=funcMaxed, group=percent)) +
  ylab(expression("Number of ES")) +
  xlab(var_lab) +
  stat_smooth(method="glm", method.args = list(family = fam_mod), lwd=0.8, fill=NA, aes(color=percent)) +
  theme_bw(base_size=18) +
  scale_color_gradient(name="Threshold", low="blue", high="red") +
  ylim(0,4) +theme(axis.text=element_text(size=24),
                   axis.title=element_text(size=24,face="bold"))
dev.off()

# Values of the diversity slope at different levels of the threshold

pdf(file=paste(res_folder,"thres2_abun.pdf",sep=""),onefile=T,width=5.5,height=4.5)
popSlopes <- ggplot(LinearSlopes, aes(x=thresholds)) +
  geom_ribbon(fill="grey50", aes(x=thresholds*100, ymin=Estimate-1.96*LinearSlopes[["Std. Error"]], 
                                 ymax=Estimate+1.96*LinearSlopes[["Std. Error"]])) +
  geom_point(aes(x=thresholds*100, y=Estimate)) +
  ylab("Biodiversity effect on ES\n") +
  xlab("\nThreshold (%)") +
  geom_abline(intercept=0, slope=0, lwd=1, linetype=2) +
  theme_bw(base_size=14)

popSlopes+theme(axis.text=element_text(size=20),
               axis.title=element_text(size=20,face="bold"))

dev.off()

# Indices
IDX <- getIndices(LinearSlopes, pop_thresh, funcMaxed ~ abun, divvar = "abun")
IDX

write.table(IDX, "abun_idx.txt",sep="\t") # abun

## Minimum and maximum number of populations sustained

for (i in seq(0.01,0.99, by=0.01)) print(c(i, round(max(predict(glm(funcMaxed ~ abun, data=subset(pop_thresh, pop_thresh$thresholds==i), family=fam_mod))),1)))

#####################################
# FG3 (Abundance of FG3 individuals)
#####################################

mpi_dat<-data.frame(wagtail=wagtail.all,
                   trout=dat$trout.mean,
                   salmon=(fish_mean$salmon0+fish_mean$salmon1),
                   dipper=dat$dipper.all,
                   FG3=exp(mpi.set$FG3))

var_lab="FG3"

vars<-qw(wagtail,trout,salmon,dipper)

pop_thresh<-getFuncsMaxed(mpi_dat, vars, threshmin=0.01, threshmax=0.99, prepend=c("FG3"), maxN=3)

LinearSlopes<-getCoefTab(funcMaxed ~ FG3, data=pop_thresh, coefVar="FG3", family=fam_mod)

## Number of populations above the threshold

pdf(file=paste(res_folder,"thres1_FG3.pdf",sep=""),onefile=T,width=6,height=4.5)
pop_thresh$percent <- 100*pop_thresh$thresholds
ggplot(data=pop_thresh, aes(x=FG3, y=funcMaxed, group=percent)) +
  ylab(expression("Number of populations")) +
  xlab(var_lab) +
  stat_smooth(method="glm", method.args = list(family = fam_mod), lwd=0.8, fill=NA, aes(color=percent)) +
  theme_bw(base_size=18) +
  scale_color_gradient(name="Threshold", low="blue", high="red") +
  ylim(0,4) +theme(axis.text=element_text(size=24),
                   axis.title=element_text(size=24,face="bold"))
dev.off()

# Values of the diversity slope at different levels of the threshold

pdf(file=paste(res_folder,"thres2_FG3.pdf",sep=""),onefile=T,width=5.5,height=4.5)
popSlopes <- ggplot(LinearSlopes, aes(x=thresholds)) +
  geom_ribbon(fill="grey50", aes(x=thresholds*100, ymin=Estimate-1.96*LinearSlopes[["Std. Error"]], 
                                 ymax=Estimate+1.96*LinearSlopes[["Std. Error"]])) +
  geom_point(aes(x=thresholds*100, y=Estimate)) +
  ylab("Biodiversity effect on ES\n") +
  xlab("\nThreshold (%)") +
  geom_abline(intercept=0, slope=0, lwd=1, linetype=2) +
  theme_bw(base_size=14)

popSlopes+theme(axis.text=element_text(size=20),
               axis.title=element_text(size=20,face="bold"))

dev.off()

# Indices
IDX <- getIndices(LinearSlopes, pop_thresh, funcMaxed ~ FG3, divvar = "FG3")
IDX

write.table(IDX, "FG3_idx.txt",sep="\t") # FG3

## Minimum and maximum number of populations sustained

for (i in seq(0.01,0.99, by=0.01)) print(c(i, round((max(predict(glm(funcMaxed ~ FG3, data=subset(pop_thresh, pop_thresh$thresholds==i), family=fam_mod)))),1)))

##############################
# TD-e (Effect trait richness)
##############################

mpi_dat<-data.frame(wagtail=wagtail.all,
                   trout=dat$trout.mean,
                   salmon=(fish_mean$salmon0+fish_mean$salmon1),
                   dipper=dat$dipper.all,
                   FDe=(mpi.set$FDe))

var_lab="FDe"

vars<-qw(wagtail,trout,salmon,dipper)

pop_thresh<-getFuncsMaxed(mpi_dat, vars, threshmin=0.01, threshmax=0.99, prepend=c("FDe"), maxN=7)

LinearSlopes<-getCoefTab(funcMaxed ~ FDe, data=pop_thresh, coefVar="FDe", family=fam_mod)

## Number of populations above the threshold

pdf(file=paste(res_folder,"thres1_FDe.pdf",sep=""),onefile=T,width=6,height=4.5)
pop_thresh$percent <- 100*pop_thresh$thresholds
ggplot(data=pop_thresh, aes(x=FDe, y=funcMaxed, group=percent)) +
  ylab(expression("Number of populations")) +
  xlab(var_lab) +
  stat_smooth(method="glm", method.args = list(family = fam_mod), lwd=0.8, fill=NA, aes(color=percent)) +
  theme_bw(base_size=18) +
  scale_color_gradient(name="Threshold", low="blue", high="red") +
  ylim(0,4) + theme(axis.text=element_text(size=24),
                   axis.title=element_text(size=24,face="bold"))

dev.off()

# Values of the diversity slope at different levels of the threshold

pdf(file=paste(res_folder,"thres2_FDe.pdf",sep=""),onefile=T,width=5.5,height=4.5)
popSlopes <- ggplot(LinearSlopes, aes(x=thresholds)) +
  geom_ribbon(fill="grey50", aes(x=thresholds*100, ymin=Estimate-1.96*LinearSlopes[["Std. Error"]], 
                                 ymax=Estimate+1.96*LinearSlopes[["Std. Error"]])) +
  geom_point(aes(x=thresholds*100, y=Estimate)) +
  ylab("Biodiversity effect on ES\n") +
  xlab("\nThreshold (%)") +
  geom_abline(intercept=0, slope=0, lwd=1, linetype=2) +
  theme_bw(base_size=14)

popSlopes+theme(axis.text=element_text(size=20),
               axis.title=element_text(size=20,face="bold"))
dev.off()

# Indices
IDX <- getIndices(LinearSlopes, pop_thresh, funcMaxed ~ FDe, divvar = "FDe")
IDX

write.table(IDX, "FDe_idx.txt",sep="\t") # FDe

########################################################
## Minimum and maximum number of population sustained
########################################################

for (i in seq(0.01,0.99, by=0.01)) print(c(i, round((max(predict(glm(funcMaxed ~ FDe, data=subset(pop_thresh, pop_thresh$thresholds==i), family=fam_mod)))),1)))

mt<-matrix(c(72, 99, 69, 74, 
                47, 58, 48, 60, 
                50, 31, 30, 33,
                23, 0, 0, 0),4,4)

colnames(mt)<-1:4

pdf(file=paste(res_folder,"es_mt.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=6,height=6)
par(mfrow=c(1,1), mar = c(5.5, 5.5, 5.5, 1))
par(cex.axis=1.9,cex.lab=1.75)


barplot(es_mt,beside=T,ylim=c(0,100),xlab="Population sizes above threshold",ylab="Relative population size", 
        col=rep(viridis(4)[4:1],4),axis.lty="solid")

legend("topright", c("FG3 biomass", "Effect trait richness" ,"Taxon richness", "Taxon abundance"), 
       fill=viridis(4)[4:1],bty="n", cex=1.4)
dev.off()

