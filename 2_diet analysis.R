##########################################################################################################
# Supplementary Data 1                                                                                    
#                                                                                                         
# R script to analyse diet preferences of riverine predators                                                            
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

assum_folder<-paste(getwd(),"/plots/assum/",sep="")
res_folder<-paste(getwd(),"/plots/",sep="")

# Loading libraries
library(reshape2)
library(car)
library(plyr)
library(sqldf)
library(nlme)
library(lme4)
library(lmerTest)
library(multcomp)
library(sandwich)
library(MuMIn)

# Dietary data from literature
diet<-read.table("aggregated_diet.txt",h=T,sep="\t", encoding = "UTF-8") # Loading fish abun data

# Overall values for raw aquatic diets

rowSums(diet[which(diet$predator=="trout"),c("FG1","FG2","FG3","FG4")])->trout_raw
rowSums(diet[which(diet$predator=="salmon"),c("FG1","FG2","FG3","FG4")])->salmon_raw
rowSums(diet[which(diet$predator=="wagtail"),c("FG1","FG2","FG3","FG4")])->wagtail_raw
rowSums(diet[which(diet$predator=="dipper"),c("FG1","FG2","FG3","FG4")])->dipper_raw

### Plots 
par(mar = c(4.3, 4.7, 2.8, 1))
# portrait, 10, 16
par(mfrow=c(2,2),cex=1, cex.axis=1.5,cex.lab=1.75,cex.main=1.5)

# Distribution of the aquatic prey contribution to each predator diet

hist(wagtail_raw,xlab="% aquatic preys",main="wagtail")
hist(dipper_raw,xlab="% aquatic preys",main="dipper")
hist(trout_raw,xlab="% aquatic preys",main="trout")
hist(salmon_raw,xlab="% aquatic preys",main="salmon")

# Median contribution of aquatic preys to each predator diet

median(wagtail_raw) #24.62
median(dipper_raw) # 94.90
median(trout_raw) # 76.64
median(salmon_raw) # 80.7

# Mean contribution of each FG to each predator diet

ddply(diet,.(predator), function(x) mean(x$FG1))
ddply(diet,.(predator), function(x) mean(x$FG2))
ddply(diet,.(predator), function(x) mean(x$FG3))
ddply(diet,.(predator), function(x) mean(x$FG4))

# Standard Errors of mean contribution of each FG to each predator diet

ddply(diet,.(predator), function(x) sd(x$FG1)/sqrt(nrow(x)))
ddply(diet,.(predator), function(x) sd(x$FG2)/sqrt(nrow(x)))
ddply(diet,.(predator), function(x) sd(x$FG3)/sqrt(nrow(x)))
ddply(diet,.(predator), function(x) sd(x$FG4)/sqrt(nrow(x)))

# Preparing data for analysis

fg1<-melt(diet[c(1:6,7)], value.name="FG1")
fg2<-melt(diet[c(1:6,8)], value.name="FG2")
fg3<-melt(diet[c(1:6,9)], value.name="FG3")
fg4<-melt(diet[c(1:6,10)], value.name="FG4")

# dataset
dat<-data.frame(diet[c(1:6)], 
                FG=rep(paste("FG",1:4,sep=""), each=nrow(diet)),
                per=car::logit(unlist(diet[,7:10])),
                per_raw=(unlist(diet[,7:10])))

dat$reference<-factor(dat$reference)

wdat<-dat[which(dat$predator=="wagtail"),]
tdat<-dat[which(dat$predator=="trout"),]
sdat<-dat[which(dat$predator=="salmon"),]
ddat<-dat[which(dat$predator=="dipper"),]

# Fitting models
mod_w<-lmer(per~FG+(1|reference), data=wdat)
mod_t<-lmer(per~FG+(1|reference), data=tdat)
mod_s<-lmer(per~FG+(1|reference), data=sdat)
mod_d<-lmer(per~FG+(1|reference), data=ddat)

par(mfrow=c(1,2))
hist(resid(mod_w,"pearson"))
plot(fitted(mod_w),resid(mod_w),main="Wagtail")

hist(resid(mod_t,"pearson"))
plot(fitted(mod_t),resid(mod_t),main="Wagtail")

hist(resid(mod_s,"pearson"))
plot(fitted(mod_s),resid(mod_s),main="Wagtail")

hist(resid(mod_d,"pearson"))
plot(fitted(mod_d),resid(mod_d),main="Wagtail")

# Main test
aov_w<-anova(mod_w, ddf="Kenward-Roger")
aov_t<-anova(mod_t, ddf="Kenward-Roger")
aov_s<-anova(mod_s, ddf="Kenward-Roger")
aov_d<-anova(mod_d, ddf="Kenward-Roger")

mod_res<-data.frame(matrix(NA, 4, 6))
names(mod_res)<-c("Predator","F","pval","R2m","R2c","n")
mod_res$Predator<-c("wagtail","trout","salmon","dipper")

round(aov_w$`Pr(>F)`,3)[1]->mod_res[1,3] # coefficient
round(aov_w$`F value`,3)[1]->mod_res[1,2] #pvalue
round(r.squaredGLMM(mod_w),2)->mod_res[1,4:5] # r2
nrow(wdat)->mod_res[1,6]

round(aov_t$`Pr(>F)`,3)[1]->mod_res[2,3] # coefficient
round(aov_t$`F value`,3)[1]->mod_res[2,2] #pvalue
round(r.squaredGLMM(mod_t),2)->mod_res[2,4:5] # r2
nrow(tdat)->mod_res[2,6]

round(aov_s$`Pr(>F)`,3)[1]->mod_res[3,3] # coefficient
round(aov_s$`F value`,3)[1]->mod_res[3,2] #pvalue
round(r.squaredGLMM(mod_s),2)->mod_res[3,4:5] # r2
nrow(sdat)->mod_res[3,6]

round(aov_d$`Pr(>F)`,3)[1]->mod_res[4,3] # coefficient
round(aov_d$`F value`,3)[1]->mod_res[4,2] #pvalue
round(r.squaredGLMM(mod_d),2)->mod_res[4,4:5] # r2
nrow(ddat)->mod_res[4,6]

write.table(mod_res, paste(res_folder,"diet_res.txt",sep=""), sep="\t")

# Post-hoc
phoc_w<-summary(glht(mod_w, mcp(FG = "Tukey")))
phoc_t<-summary(glht(mod_t, mcp(FG = "Tukey")))
phoc_s<-summary(glht(mod_s, mcp(FG = "Tukey")))
phoc_d<-summary(glht(mod_d, mcp(FG = "Tukey")))

# post-hoc letters
let_w<-cld(phoc_w, decreasing = T)
let_t<-cld(phoc_t, decreasing = T)
let_s<-cld(phoc_s, decreasing = T)
let_d<-cld(phoc_d, decreasing = T)

### Plots 

# General diet for each predator

pdf(file="diets_new.pdf",,useDingbats=FALSE,onefile=T,width=8,height=9)
par(mar = c(4, 5, 4, 1))

# portrait, 10, 16
par(mfrow=c(2,2),cex=1, cex.axis=1.3,cex.lab=1.5,cex.main=1.5)

#Overall
boxplot(per_raw~FG,data=wdat, ylab="", xlab="", main="wagtail",ylim=c(0,100),
        col=c("blue","dark green","red","violet"),border=c("blue","dark green","red","violet"))
mtext("a", line = 1.8, adj = -0.1, cex = 2, font = 2)
text(x=c(1,2,3,4),y=95, as.character(let_w$mcletters$Letters), cex=1.75)

boxplot(per_raw~FG,data=tdat, ylab="",xlab="",main="Brown trout",ylim=c(0,100),
        col=c("blue","dark green","red","violet"),border=c("blue","dark green","red","violet"))
mtext("b", line = 1.8, adj = -0.1, cex = 2, font = 2)
text(x=c(1,2,3,4),y=95, as.character(let_t$mcletters$Letters), cex=1.75)

boxplot(per_raw~FG,data=sdat, ylab="",xlab="",main="Atlantic salmon",ylim=c(0,100),
        col=c("blue","dark green","red","violet"),border=c("blue","dark green","red","violet"))
mtext("c", line = 1.8, adj = -0.1, cex = 2, font = 2)
text(x=c(1,2,3,4),y=95, as.character(let_s$mcletters$Letters), cex=1.75)

boxplot(per_raw~FG,data=ddat, ylab="",xlab="",main="Eurasian dipper",ylim=c(0,100),
        col=c("blue","dark green","red","violet"),border=c("blue","dark green","red","violet"))
mtext("d", line = 1.8, adj = -0.1, cex = 2, font = 2)
text(x=c(1,2,3,4),y=95, as.character(let_d$mcletters$Letters), cex=1.75)

dev.off()

