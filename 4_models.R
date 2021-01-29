##########################################################################################################
#                                                                                                         
# Supplementary Data 1                                                                                    
#                                                                                                         
# R script to relate biodiversity dimensions with single and multiple predator population sizes                                                                                                                          #
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

# Loading libraries
library(plyr)
library(nlme)
library(MuMIn)
library(arm)
library(viridis)
library(usdm)

# Loading datasets
bird.set<-read.table("bird_set.txt",h=T,sep="\t")
trout.set<-read.table("trout_set.txt",h=T,sep="\t")
salmon.set<-read.table("salmon_set.txt",h=T,sep="\t")
mpi.set<-read.table("mpi_set.txt",h=T,sep="\t")

# Transforming variables
q<-sapply(bird.set,class)=="numeric" | sapply(bird.set,class)=="integer"# selecting quantitative variables

par(mfrow=c(2,2))
for (i in which(q==T)) hist(bird.set[,i], main=names(bird.set)[i])
par(mfrow=c(1,1))

q<-sapply(trout.set,class)=="numeric" | sapply(trout.set,class)=="integer"# selecting quantitative variables
par(mfrow=c(2,2))
for (i in which(q==T)) hist(trout.set[,i], main=names(trout.set)[i])
par(mfrow=c(1,1))

q<-sapply(salmon.set,class)=="numeric" | sapply(salmon.set,class)=="integer"# selecting quantitative variables

par(mfrow=c(3,3))
for (i in which(q==T)) hist(salmon.set[,i], main=names(salmon.set)[i])
par(mfrow=c(1,1))

q<-sapply(mpi.set,class)=="numeric" | sapply(mpi.set,class)=="integer"# selecting quantitative variables
par(mfrow=c(3,3))
for (i in which(q==T)) hist(mpi.set[,i], main=names(mpi.set)[i])
par(mfrow=c(1,1))

# Variable transformation
sqrt(trout.set$trout)->trout.set$trout
sqrt(salmon.set$salmon)->salmon.set$salmon

sqrt(bird.set$FDe)-> bird.set$FDe
sqrt(trout.set$FDe)-> trout.set$FDe
sqrt(salmon.set$FDe)-> salmon.set$FDe
sqrt(mpi.set$FDe)-> mpi.set$FDe

# Creating FGflow variable

#bird.set$FGflow<-sqrt(I(exp(bird.set$FG2)+exp(bird.set$FG3)))
#trout.set$FGflow<-sqrt(I(exp(trout.set$FG2)+exp(trout.set$FG3)))
#salmon.set$FGflow<-sqrt(I(exp(salmon.set$FG2)+exp(salmon.set$FG3)))
#mpi.set$FGflow<-sqrt(I(exp(mpi.set$FG2)+exp(mpi.set$FG3)))

# standarising variables

bird.set2<-bird.set
trout.set2<-trout.set
salmon.set2<-salmon.set
mpi.set2<-mpi.set

bird.set2[,-c(1:4)]<-scale(bird.set[,-c(1:4)])
trout.set2[,-c(1:4)]<-scale(trout.set[,-c(1:4)])
salmon.set2[,-c(1:4)]<-scale(salmon.set[,-c(1:4)])
mpi.set2[,-c(1:2)]<-scale(mpi.set[,-c(1:2)])


# Birds
b.cor<-as.dist(round(cor(bird.set2[,-c(1:4)]),2))
write.table(as.matrix(b.cor),paste(res_folder,"bird_cor.txt", sep=""),sep="\t",row.names=F)

vifstep(bird.set2[,c("FG2","FDr")],th=2)
vifstep(bird.set2[,c("FG3","FDr","alt")],th=2)
vifstep(bird.set2[,c("FG2","FDe","ph","alt")],th=2)
vifstep(bird.set2[,c("FG3","FDe","ph","alt")],th=2)
vifstep(bird.set2[,c("ric","ph","alt")],th=2)
vifstep(bird.set2[,c("abun","ph","alt")],th=2)

# Trout
t.cor<-as.dist(round(cor(trout.set2[,-c(1:4,13)]),2))
write.table(as.matrix(t.cor),paste(res_folder,"trout_cor.txt", sep=""),sep="\t",row.names=F)

vifstep(trout.set2[,c("FG2","FDe","alt")],th=2)
vifstep(trout.set2[,c("FG2","FDr","alt")],th=2)
vifstep(trout.set2[,c("FG3","FDe","alt")],th=2)
vifstep(trout.set2[,c("FG3","FDr","alt")],th=2)

vifstep(trout.set2[,c("FDe","abun")],th=2)
vifstep(trout.set2[,c("FDe","ph","alt")],th=2)
vifstep(trout.set2[,c("FDr","abun")],th=2)
vifstep(trout.set2[,c("FDr","ph","alt")],th=2)

vifstep(trout.set2[,c("abun","ric","alt")],th=2)
vifstep(trout.set2[,c("abun","ph","alt")],th=2)
vifstep(trout.set2[,c("ric","ph","alt")],th=2)

# Salmon
s.cor<-as.dist(round(cor(salmon.set2[,-c(1:4,13)]),2))
write.table(as.matrix(s.cor),paste(res_folder,"salmon_cor.txt", sep=""),sep="\t",row.names=F)

vifstep(salmon.set2[,c("FG2","FDe","alt")],th=2)
vifstep(salmon.set2[,c("FG2","FDr","alt")],th=2)
vifstep(salmon.set2[,c("FG2","ric","alt")],th=2)
vifstep(salmon.set2[,c("FG3","FDe","alt")],th=2)
vifstep(salmon.set2[,c("FG3","FDr","alt")],th=2)
vifstep(salmon.set2[,c("FG3","ric")],th=2)

vifstep(salmon.set2[,c("FDe","abun")],th=2)
vifstep(salmon.set2[,c("FDe","ph","alt")],th=2)
vifstep(salmon.set2[,c("FDr","abun")],th=2)
vifstep(salmon.set2[,c("FDr","ph","alt")],th=2)

vifstep(salmon.set2[,c("abun","ric")],th=2)
vifstep(salmon.set2[,c("abun","alt")],th=2)
vifstep(salmon.set2[,c("ric","ph","alt")],th=2)

#mpi
mpi.cor<-as.dist(round(cor(mpi.set2[,-c(1:2,11)]),2))
write.table(as.matrix(mpi.cor),paste(res_folder,"mpi_cor.txt", sep=""),sep="\t",row.names=F)


vifstep(mpi.set2[,c("FG2","FDe","alt")],th=2)
vifstep(mpi.set2[,c("FG2","FDr","alt")],th=2)
vifstep(mpi.set2[,c("FG2","ric","alt")],th=2)
vifstep(mpi.set2[,c("FG3","FDe","alt")],th=2)
vifstep(mpi.set2[,c("FG3","FDr","alt")],th=2)
vifstep(mpi.set2[,c("FG3","ric")],th=2)

vifstep(mpi.set2[,c("FDe","abun")],th=2)
vifstep(mpi.set2[,c("FDe","ph","alt")],th=2)
vifstep(mpi.set2[,c("FDr","abun")],th=2)
vifstep(mpi.set2[,c("FDr","ph","alt")],th=2)

vifstep(mpi.set2[,c("abun","ric")],th=2)
vifstep(mpi.set2[,c("abun","alt")],th=2)
vifstep(mpi.set2[,c("ric","ph","alt")],th=2)


####################################################################################################
# Ranking biodiversity metrics according to their capacity to predict ES
####################################################################################################

# Wagtail

    mod1<-glm(wagtail~FG2,data=bird.set2, family="poisson")
    mod2<-glm(wagtail~FG3,data=bird.set2, family="poisson")
    mod3<-glm(wagtail~FDe,data=bird.set2, family="poisson")
    mod4<-glm(wagtail~FDr,data=bird.set2, family="poisson")
    mod5<-glm(wagtail~ric,data=bird.set2, family="poisson")
    mod6<-glm(wagtail~abun,data=bird.set2, family="poisson")
    mod7<-glm(wagtail~ph+alt,data=bird.set2, family="poisson")
    
    mod8<-glm(wagtail~FG2+FDe,data=bird.set2, family="poisson")
    mod9<-glm(wagtail~FG2+FDr,data=bird.set2, family="poisson")
    mod10<-glm(wagtail~FG2+alt+ph,data=bird.set2, family="poisson")
    mod11<-glm(wagtail~FG2+FDe+alt+ph,data=bird.set2, family="poisson")
    mod12<-glm(wagtail~FG3+FDe,data=bird.set2, family="poisson")
    mod13<-glm(wagtail~FG3+FDr,data=bird.set2, family="poisson")
    mod14<-glm(wagtail~FG3+alt+ph,data=bird.set2, family="poisson")
    mod15<-glm(wagtail~FG3+FDe+alt+ph,data=bird.set2, family="poisson")
    
    mod16<-glm(wagtail~FDe+alt+ph,data=bird.set2, family="poisson")
    mod17<-glm(wagtail~FDr+alt+ph,data=bird.set2, family="poisson")
    
    mod16<-glm(wagtail~ric+alt+ph,data=bird.set2, family="poisson")
    mod17<-glm(wagtail~abun+alt+ph,data=bird.set2, family="poisson")
    
    mod.list<-list(mod1, mod2, mod3, mod4, mod5, mod6,
                   mod7,mod8, mod9, mod10, mod11,mod12,
                   mod13, mod14, mod15, mod16, mod17)
    
    w.res<-model.sel(mod.list, extra=c(r2m=function(x) 1-(deviance(x)/summary(x)$null.deviance)))
 
## Dipper
    
    mod1<-glm(dipper~FG2,data=bird.set2, family="poisson")
    mod2<-glm(dipper~FG3,data=bird.set2, family="poisson")
    mod3<-glm(dipper~FDe,data=bird.set2, family="poisson")
    mod4<-glm(dipper~FDr,data=bird.set2, family="poisson")
    mod5<-glm(dipper~ric,data=bird.set2, family="poisson")
    mod6<-glm(dipper~abun,data=bird.set2, family="poisson")
    mod7<-glm(dipper~ph+alt,data=bird.set2, family="poisson")
    
    mod8<-glm(dipper~FG2+FDe,data=bird.set2, family="poisson")
    mod9<-glm(dipper~FG2+FDr,data=bird.set2, family="poisson")
    mod10<-glm(dipper~FG2+alt+ph,data=bird.set2, family="poisson")
    mod11<-glm(dipper~FG2+FDe+alt+ph,data=bird.set2, family="poisson")
    mod12<-glm(dipper~FG3+FDe,data=bird.set2, family="poisson")
    mod13<-glm(dipper~FG3+FDr,data=bird.set2, family="poisson")
    mod14<-glm(dipper~FG3+alt+ph,data=bird.set2, family="poisson")
    mod15<-glm(dipper~FG3+FDe+alt+ph,data=bird.set2, family="poisson")
    
    mod16<-glm(dipper~FDe+alt+ph,data=bird.set2, family="poisson")
    mod17<-glm(dipper~FDr+alt+ph,data=bird.set2, family="poisson")
    
    mod16<-glm(dipper~ric+alt+ph,data=bird.set2, family="poisson")
    mod17<-glm(dipper~abun+alt+ph,data=bird.set2, family="poisson")
    mod18<-glm(dipper~abun+FDe+alt+ph,data=bird.set2, family="poisson")
    
    mod.list<-list(mod1, mod2, mod3, mod4, mod5, mod6,
                   mod7,mod8, mod9, mod10, mod11,mod12,
                   mod13, mod14, mod15, mod16, mod17, mod18)
    
    d.res<-model.sel(mod.list, extra=c(r2m=function(x) 1-(deviance(x)/summary(x)$null.deviance)))

    
## Trout
    
    mod1<-lme(trout~FG2,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod2<-lme(trout~FG3,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod3<-lme(trout~FDe,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod4<-lme(trout~FDr,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod5<-lme(trout~ric,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod6<-lme(trout~abun,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod7<-lme(trout~ph+alt,data=trout.set2, random = ~ 1 | sites, method="ML")
    
    mod8<-lme(trout~FG2+FDe,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod9<-lme(trout~FG2+FDr,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod10<-lme(trout~FG2+ric,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod11<-lme(trout~FG2+alt,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod12<-lme(trout~FG2+FDe+alt,data=trout.set2, random = ~ 1 | sites, method="ML")
    
    mod13<-lme(trout~FG3+FDe,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod14<-lme(trout~FG3+FDr,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod15<-lme(trout~FG3+ric,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod16<-lme(trout~FG3+alt,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod17<-lme(trout~FG3+FDe+alt,data=trout.set2, random = ~ 1 | sites, method="ML")
    
    mod18<-lme(trout~FDe+alt+ph,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod19<-lme(trout~FDr+alt+ph,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod20<-lme(trout~FDe+abun,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod21<-lme(trout~FDr+abun,data=trout.set2, random = ~ 1 | sites, method="ML")
    
    mod22<-lme(trout~ric+alt+ph,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod23<-lme(trout~abun+alt+ph,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod24<-lme(trout~abun+ric,data=trout.set2, random = ~ 1 | sites, method="ML")
    mod25<-lme(trout~abun+ric+alt,data=trout.set2, random = ~ 1 | sites, method="ML")
    
    mod.list<-list(mod1, mod2, mod3, mod4, mod5, mod6,
                   mod7,mod8, mod9, mod10, mod11,mod12,
                   mod13, mod14, mod15, mod16, mod17, mod18,
                   mod19, mod20, mod21, mod22, mod23, mod24,
                   mod25)
    
    t.res<-model.sel(mod.list, extra=c(r2m=function(x) r.squaredGLMM(x)[1],r2c=function(x) r.squaredGLMM(x)[2]))
    
## Salmon
    
    mod1<-lme(salmon~FG2,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod2<-lme(salmon~FG3,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod3<-lme(salmon~FDe,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod4<-lme(salmon~FDr,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod5<-lme(salmon~ric,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod6<-lme(salmon~abun,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod7<-lme(salmon~ph+alt,data=salmon.set2, random = ~ 1 | sites, method="ML")
    
    mod8<-lme(salmon~FG2+FDe,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod9<-lme(salmon~FG2+FDr,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod10<-lme(salmon~FG2+ric,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod11<-lme(salmon~FG2+alt,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod12<-lme(salmon~FG2+FDe+alt,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod13<-lme(salmon~FG2+FDr+alt,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod14<-lme(salmon~FG2+ric+alt,data=salmon.set2, random = ~ 1 | sites, method="ML")
    
    mod15<-lme(salmon~FG3+FDe,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod16<-lme(salmon~FG3+FDr,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod17<-lme(salmon~FG3+ric,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod18<-lme(salmon~FG3+alt,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod19<-lme(salmon~FG3+FDe+alt,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod20<-lme(salmon~FG3+FDr+alt,data=salmon.set2, random = ~ 1 | sites, method="ML")
    
    mod21<-lme(salmon~FDe+alt+ph,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod22<-lme(salmon~FDr+alt+ph,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod23<-lme(salmon~FDe+abun,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod24<-lme(salmon~FDr+abun,data=salmon.set2, random = ~ 1 | sites, method="ML")
    
    mod25<-lme(salmon~abun+ric,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod26<-lme(salmon~ric+alt,data=salmon.set2, random = ~ 1 | sites, method="ML")
    mod27<-lme(salmon~abun+alt,data=salmon.set2, random = ~ 1 | sites, method="ML")
    
    mod.list<-list(mod1, mod2, mod3, mod4, mod5, mod6,
                   mod7,mod8, mod9, mod10, mod11,mod12,
                   mod13, mod14, mod15, mod16, mod17, mod18,
                   mod19, mod20, mod21, mod22, mod23, mod24,
                   mod25, mod26, mod27)
    
    s.res<-model.sel(mod.list, extra=c(r2m=function(x) r.squaredGLMM(x)[1],r2c=function(x) r.squaredGLMM(x)[2]))
    
    ## mpi
    
    mod1<-lm(mpi~FG2,data=mpi.set2)
    mod2<-lm(mpi~FG3,data=mpi.set2)
    mod3<-lm(mpi~FDe,data=mpi.set2)
    mod4<-lm(mpi~FDr,data=mpi.set2)
    mod5<-lm(mpi~ric,data=mpi.set2)
    mod6<-lm(mpi~abun,data=mpi.set2)
    mod7<-lm(mpi~ph+alt,data=mpi.set2)
    
    mod8<-lm(mpi~FG2+FDe,data=mpi.set2)
    mod9<-lm(mpi~FG2+FDr,data=mpi.set2)
    mod10<-lm(mpi~FG2+ric,data=mpi.set2)
    mod11<-lm(mpi~FG2+alt,data=mpi.set2)
    mod12<-lm(mpi~FG2+FDe+alt,data=mpi.set2)
    mod13<-lm(mpi~FG2+FDr+alt,data=mpi.set2)
    mod14<-lm(mpi~FG2+ric+alt,data=mpi.set2)
    
    mod15<-lm(mpi~FG3+FDe,data=mpi.set2)
    mod16<-lm(mpi~FG3+FDr,data=mpi.set2)
    mod17<-lm(mpi~FG3+ric,data=mpi.set2)
    mod18<-lm(mpi~FG3+alt,data=mpi.set2)
    mod19<-lm(mpi~FG3+FDe+alt,data=mpi.set2)
    mod20<-lm(mpi~FG3+FDr+alt,data=mpi.set2)
    
    mod21<-lm(mpi~FDe+alt+ph,data=mpi.set2)
    mod22<-lm(mpi~FDr+alt+ph,data=mpi.set2)
    mod23<-lm(mpi~FDe+abun,data=mpi.set2)
    mod24<-lm(mpi~FDr+abun,data=mpi.set2)
    
    mod25<-lm(mpi~abun+ric,data=mpi.set2)
    mod26<-lm(mpi~ric+alt,data=mpi.set2)
    mod27<-lm(mpi~abun+alt,data=mpi.set2)
    
    mod.list<-list(mod1, mod2, mod3, mod4, mod5, mod6,
                   mod7,mod8, mod9, mod10, mod11,mod12,
                   mod13, mod14, mod15, mod16, mod17, mod18,
                   mod19, mod20, mod21, mod22, mod23, mod24,
                   mod25, mod26, mod27)
    
    mpi.res<-model.sel(mod.list, extra=c(r2m=function(x) summary(x)$adj.r.squared))
    
     
 
  # model assumptions
  
  for (b in 1:9) 
  {
    mod.list[[b]]->mod
    pdf(file=paste(assum_folder,resp[[a]],"_mod",b,"_assum",".pdf",sep=""),onefile=T,width=6.24,height=4.05)
    par(mfrow=c(1,2),cex=1, cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
    if (mod.type[a]=="lme") m.resid<-resid(mod,type="pearson") else m.resid<-resid(mod)
    hist(m.resid,main="",xlab="residuals")
    mtext(paste(resp[[a]],"_mod",b),line=0.5,at=4,cex=2)
    plot(fitted(mod),m.resid,xlab="fitted values",ylab="residuals",main="")  
    dev.off()
  }
  
    


# saving results
write.table(w.res,paste(res_folder,"w_mod_res.txt", sep=""),sep="\t",row.names=F)
write.table(d.res,paste(res_folder,"d_mod_res.txt", sep=""),sep="\t",row.names=F)
write.table(t.res,paste(res_folder,"t_mod_res.txt", sep=""),sep="\t",row.names=F)
write.table(s.res,paste(res_folder,"s_mod_res.txt", sep=""),sep="\t",row.names=F)
write.table(mpi.res,paste(res_folder,"mpi_mod_res.txt", sep=""),sep="\t",row.names=F)

save.image(paste(getwd(),"model_results.RData",sep=""))

# Creating barplots for r2

fi.r2<-fd.r2<-fe.r2<-env.r2<-rep(NA, 5)

for (i in 1:5){
  
  res.df<-res[[i]]
  
  fi.r2[i]<-max(res.df[which(rownames(res.df)==1 | rownames(res.df)==2 | rownames(res.df)==3),"r2m"])
  fd.r2[i]<-max(res.df[which(rownames(res.df)==4 | rownames(res.df)==5 | rownames(res.df)==6),"r2m"])
  fe.r2[i]<-max(res.df[which(rownames(res.df)==7 | rownames(res.df)==8),"r2m"])
  env.r2[i]<-max(res.df[which(rownames(res.df)==9),"r2m"])

}

r2<-data.frame(fi.r2, fd.r2, fe.r2, env.r2)
rownames(r2)<-resp

es.lab<-c("Wagtail territories", 
          "Dipper territories",
          expression(paste("Salmon biomass (", paste(g),")")),
          expression(paste("Trout biomass (", paste(g),")")),
          "Multi-predator index")

pdf(file=paste(res_folder,"r2_col2.pdf",sep=""),,useDingbats=FALSE,onefile=T,
    width=11,height=19)
par(mfrow=c(3,2), mar = c(5.5, 5.5, 5.5, 1))
par(pch=20,cex=1,cex.axis=1.9,cex.lab=1.9,cex.main=2)


for (i in 1:5) {
  
  barplot(as.numeric(r2[i,])[4:1],
          xlab=expression(italic(r)^2),
          xlim=c(0, max(r2)),
          names.arg=c("Key", "TD", "Mass", "Env")[4:1], 
          col=c("light grey", viridis(3)),
          main=es.lab[i],
          horiz=T)
  
  mtext(letters[i], line = 1.8, adj = -0.1, cex = 3, font = 2)
  
}

dev.off()


resp_plot<-c("Wagtail territories", 
             "Dipper territories",
             expression(paste("Salmon biomass (", paste(g),")")),
             expression(paste("Trout biomass (", paste(g),")")),
             "Multi-predator index")


df.list2<-list(bird.set, bird.set, salmon.set, trout.set, mpi.set)

## Single population plots ##

pdf(file=paste(res_folder,"model_plots_es3.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=12,height=13.75)
par(mfrow=c(4,4), mar = c(4, 3.8, 3, 1))
par(pch=20,cex=1,cex.axis=1,cex.lab=1.2,cex.main=1.25)
text.cex=1

a=0

for (i in 1:4) {
  
    a=a+1
    
    df<-df.list[[a]]
    df2<-df.list2[[a]]
    
    y<-df[,resp[[a]]]
    
    fr3.seq<-seq(min(df2$FG3),max(df2$FG3),length.out=1000)
    er.seq<-seq(min(df2$FDe),max(df2$FDe),length.out=1000)
    ric.seq<-seq(min(df2$ric),max(df2$ric),length.out=1000)
    ab.seq<-seq(min(df2$abun),max(df2$abun),length.out=1000)
    
    if (mod.type[a]=="poisson") {
      mod2<-glm(y~FG3,data=df, family="poisson")
      mod4<-glm(y~FDe,data=df, family="poisson")
      mod7<-glm(y~ric,data=df, family="poisson")
      mod8<-glm(y~abun,data=df, family="poisson")
      
      fi.sum<-summary(mod2)$coef
      fd.sum<-summary(mod4)$coef
      fe1.sum<-summary(mod7)$coef
      fe2.sum<-summary(mod8)$coef
    } 
    
    if (mod.type[a]=="lm") {
      mod2<-lm(y~FG3,data=df)
      mod4<-lm(y~FDe,data=df)
      mod7<-lm(y~ric,data=df)
      mod8<-lm(y~abun,data=df)
      
      fi.sum<-summary(mod2)$coef
      fd.sum<-summary(mod4)$coef
      fe1.sum<-summary(mod7)$coef
      fe2.sum<-summary(mod8)$coef
    } 
    
    if (mod.type[a]=="lme") {
      mod2<-lme(y~FG3,data=df, random = ~ 1 | sites, method="REML")
      mod4<-lme(y~FDe,data=df, random = ~ 1 | sites, method="REML")
      mod7<-lme(y~ric,data=df, random = ~ 1 | sites, method="REML")
      mod8<-lme(y~abun,data=df, random = ~ 1 | sites, method="REML")
    
      fi.sum<-summary(mod2)$tTable[,-3]
      fd.sum<-summary(mod4)$tTable[,-3]
      fe1.sum<-summary(mod7)$tTable[,-3]
      fe2.sum<-summary(mod8)$tTable[,-3]
      
      y<-y^2
    }
    
    SES1 <- bquote(SES == .(round(fi.sum[2,1], 2)))
    pval1=round(fi.sum[2,4], 3)
    if(pval1<0.001) ppval1 <- bquote(italic(P) <.(0.001)) else ppval1 <- bquote(italic(P) == .(pval1))
    fi.R2 <-bquote(italic(r)^2 == .(round(r.squaredGLMM(mod2)[1,1], 2)))
    
    SES2 <- bquote(SES == .(round(fd.sum[2,1], 2)))
    pval2=round(fd.sum[2,4], 3)
    if(pval2<0.001) ppval2 <- bquote(italic(P) <.(0.001)) else ppval2 <- bquote(italic(P) == .(pval2))
    fd.R2 <-bquote(italic(r)^2 == .(round(r.squaredGLMM(mod4)[1,1], 2)))
    
    SES3 <- bquote(SES == .(round(fe1.sum[2,1], 2)))
    pval3=round(fe1.sum[2,4], 3)
    if(pval3<0.001) ppval3 <- bquote(italic(P) <.(0.001)) else ppval3 <- bquote(italic(P) == .(pval3))
    fe1.R2 <-bquote(italic(r)^2 == .(round(r.squaredGLMM(mod7)[1,1], 2)))
    
    SES4 <- bquote(SES == .(round(fe2.sum[2,1], 2)))
    pval4=round(fe2.sum[2,4], 3)
    if(pval4<0.001) ppval4 <- bquote(italic(P) <.(0.001)) else ppval4 <- bquote(italic(P) == .(pval4))
    fe2.R2 <-bquote(italic(r)^2 == .(round(r.squaredGLMM(mod8)[1,1], 2)))
    
    #v.col=c("red","blue","orange")
    v.col=viridis(3)[3:1]
    text.cex=1
    
    ### FR3
    plot(exp(df2$FG3), y, xlab="",ylab="", col=grey.colors(1,0.7), pch=16)
    mtext(resp_plot[i], side=2, line=2, cex=1.3)
    if(i==4) mtext("FG3 abundance", side=1, line=2.25, cex=1.3)
    if(i==1)  mtext("Key", side=3, cex=1.5)
    if(pval1<0.05 & mod.type[a]=="lm") lines(exp(fr3.seq),predict(mod2,data.frame(FG3=scale(fr3.seq))),lwd=3,col=v.col[1])
    if(pval1<0.05 & mod.type[a]=="lme") lines(exp(fr3.seq),predict(mod2,data.frame(FG3=scale(fr3.seq)),level=0)^2,lwd=3,col=v.col[1])
    if(pval1<0.05 & mod.type[a]=="poisson") lines(exp(fr3.seq),exp(predict(mod2,data.frame(FG3=scale(fr3.seq)))),lwd=3,col=v.col[1])
    
    mtext(SES1, side = 3, adj = 0.975, line = -1.05, cex = text.cex)
    mtext(ppval1, side = 3, adj = 0.975, line = -2.05, cex = text.cex)
    mtext(fi.R2, side = 3, adj = 0.975, line = -3.05, cex = text.cex)
    mtext(letters[i], line = 1.8, adj = -0.1, cex = 2, font = 2)
    
        ### ER
    plot(df2$FDe, y, xlab="",ylab="", col=grey.colors(1,0.7), pch=16)
    if(i==1)  mtext("TD", side=3, cex=1.5)
    if(i==4) mtext("Effect trait richness", side=1, line=2.25, cex=1.3)
    if(pval2<0.05 & mod.type[a]=="lm") lines(er.seq,predict(mod4,data.frame(FDe=scale(er.seq))),lwd=3,col=v.col[2])
    if(pval2<0.05 & mod.type[a]=="lme") lines(er.seq,predict(mod4,data.frame(FDe=scale(er.seq)),level=0)^2,lwd=3,col=v.col[2])
    if(pval2<0.05 & mod.type[a]=="poisson") lines(er.seq,exp(predict(mod4,data.frame(FDe=scale(er.seq)))),lwd=3,col=v.col[2])
    
    mtext(SES2, side = 3, adj = 0.975, line = -1.05, cex = text.cex)
    mtext(ppval2, side = 3, adj = 0.975, line = -2.05, cex = text.cex)
    mtext(fd.R2, side = 3, adj = 0.975, line = -3.05, cex = text.cex)
    
    ### ric
    plot(df2$ric, y, xlab="",ylab="", col=grey.colors(1,0.7), pch=16)
    if(i==1)  mtext("Mass", side=3, cex=1.5)
    if(i==4) mtext("Taxon richness", side=1, line=2.25, cex=1.3)
    if(pval3<0.05 & mod.type[a]=="lm") lines(ric.seq,predict(mod7,data.frame(ric=scale(ric.seq))),lwd=3,col=v.col[3])
    if(pval3<0.05 & mod.type[a]=="lme") lines(ric.seq,predict(mod7,data.frame(ric=scale(ric.seq)),level=0)^2,lwd=3,col=v.col[3])
    if(pval3<0.05 & mod.type[a]=="poisson") lines(ric.seq,exp(predict(mod7,data.frame(ric=scale(ric.seq)))),lwd=3,col=v.col[3])
    
    mtext(SES3, side = 3, adj = 0.975, line = -1.05, cex = text.cex)
    mtext(ppval3, side = 3, adj = 0.975, line = -2.05, cex = text.cex)
    mtext(fe1.R2, side = 3, adj = 0.975, line = -3.05, cex = text.cex)
    
    
    ### abun
    plot(exp(df2$abun), y, xlab="",ylab="", col=grey.colors(1,0.7), pch=16)
    if(i==1)  mtext("Mass", side=3, cex=1.5)
    if(i==4) mtext("Taxon abundance", side=1, line=2.25, cex=1.3)
    
    if(pval4<0.05 & mod.type[a]=="lm") lines(exp(ab.seq),predict(mod8,data.frame(abun=scale(ab.seq))),lwd=3,col=v.col[3])
    if(pval4<0.05 & mod.type[a]=="lme") lines(exp(ab.seq),predict(mod8,data.frame(abun=scale(ab.seq)),level=0)^2,lwd=3,col=v.col[3])
    if(pval4<0.05 & mod.type[a]=="poisson") lines(exp(ab.seq),exp(predict(mod8,data.frame(abun=scale(ab.seq)))),lwd=3,col=v.col[3])
    
    mtext(SES4, side = 3, adj = 0.975, line = -1.05, cex = text.cex, font=2)
    mtext(ppval4, side = 3, adj = 0.975, line = -2.05, cex = text.cex, font=2)
    mtext(fe2.R2, side = 3, adj = 0.975, line = -3.05, cex = text.cex, font=2)
    
    
}

dev.off()

## Multi-predator index plot ##

pdf(file=paste(res_folder,"model_plots_multi.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=12,height=3.5)
par(mfrow=c(1,4), mar = c(5, 5, 4, 1))
par(pch=20,cex=1,cex.axis=1,cex.lab=1.2,cex.main=1.25)
text.cex=1

i=5
a=5
  
  df<-df.list[[a]]
  df2<-df.list2[[a]]
  
  y<-df[,resp[[a]]]
  
  fr3.seq<-seq(min(df2$FG3),max(df2$FG3),length.out=1000)
  er.seq<-seq(min(df2$FDe),max(df2$FDe),length.out=1000)
  ric.seq<-seq(min(df2$ric),max(df2$ric),length.out=1000)
  ab.seq<-seq(min(df2$abun),max(df2$abun),length.out=1000)
  
    mod2<-lm(y~FG3,data=df)
    mod4<-lm(y~FDe,data=df)
    mod7<-lm(y~ric,data=df)
    mod8<-lm(y~abun,data=df)
    
    fi.sum<-summary(mod2)$coef
    fd.sum<-summary(mod4)$coef
    fe1.sum<-summary(mod7)$coef
    fe2.sum<-summary(mod8)$coef
 
  SES1 <- bquote(SES == .(round(fi.sum[2,1], 2)))
  pval1=round(fi.sum[2,4], 3)
  if(pval1<0.001) ppval1 <- bquote(italic(P) <.(0.001)) else ppval1 <- bquote(italic(P) == .(pval1))
  fi.R2 <-bquote(italic(r)^2 == .(round(r.squaredGLMM(mod2)[1,1], 2)))
  
  SES2 <- bquote(SES == .(round(fd.sum[2,1], 2)))
  pval2=round(fd.sum[2,4], 3)
  if(pval2<0.001) ppval2 <- bquote(italic(P) <.(0.001)) else ppval2 <- bquote(italic(P) == .(pval2))
  fd.R2 <-bquote(italic(r)^2 == .(round(r.squaredGLMM(mod4)[1,1], 2)))
  
  SES3 <- bquote(SES == .(round(fe1.sum[2,1], 2)))
  pval3=round(fe1.sum[2,4], 3)
  if(pval3<0.001) ppval3 <- bquote(italic(P) <.(0.001)) else ppval3 <- bquote(italic(P) == .(pval3))
  fe1.R2 <-bquote(italic(r)^2 == .(round(r.squaredGLMM(mod7)[1,1], 2)))
  
  SES4 <- bquote(SES == .(round(fe2.sum[2,1], 2)))
  pval4=round(fe2.sum[2,4], 3)
  if(pval4<0.001) ppval4 <- bquote(italic(P) <.(0.001)) else ppval4 <- bquote(italic(P) == .(pval4))
  fe2.R2 <-bquote(italic(r)^2 == .(round(r.squaredGLMM(mod8)[1,1], 2)))
  
  #v.col=c("red","blue","orange")
  v.col=viridis(3)[3:1]
  text.cex=1
  
  ### FR3
  plot(exp(df2$FG3), y, xlab="FG3 abundance",ylab=resp_plot[i], col=grey.colors(1,0.7), pch=16)
  mtext("Key", side=3, cex=1.5)
  if(pval1<0.05 & mod.type[a]=="lm") lines(exp(fr3.seq),predict(mod2,data.frame(FG3=scale(fr3.seq))),lwd=3,col=v.col[1])
  
  mtext(SES1, side = 3, adj = 0.975, line = -6.5, cex = text.cex)
  mtext(ppval1, side = 3, adj = 0.975, line = -7.5, cex = text.cex)
  mtext(fi.R2, side = 3, adj = 0.975, line = -8.5, cex = text.cex)
  
  ### ER
  plot(df2$FDe, y, xlab="Effect trait richness",ylab="", col=grey.colors(1,0.7), pch=16)
  mtext("TD", side=3, cex=1.5)
  if(pval2<0.05 & mod.type[a]=="lm") lines(er.seq,predict(mod4,data.frame(FDe=scale(er.seq))),lwd=3,col=v.col[2])
  
  mtext(SES2, side = 3, adj = 0.975, line = -6.5, cex = text.cex)
  mtext(ppval2, side = 3, adj = 0.975, line = -7.5, cex = text.cex)
  mtext(fd.R2, side = 3, adj = 0.975, line = -8.5, cex = text.cex)
  
  ### ric
  plot(df2$ric, y, xlab="Taxon richness",ylab="", col=grey.colors(1,0.7), pch=16)
  mtext("Mass", side=3, cex=1.5)
  if(pval3<0.05 & mod.type[a]=="lm") lines(ric.seq,predict(mod7,data.frame(ric=scale(ric.seq))),lwd=3,col=v.col[3])
  i
  mtext(SES3, side = 3, adj = 0.975, line = -6.5, cex = text.cex)
  mtext(ppval3, side = 3, adj = 0.975, line = -7.5, cex = text.cex)
  mtext(fe1.R2, side = 3, adj = 0.975, line = -8.5, cex = text.cex)
  
  
  ### abun
  plot(exp(df2$abun), y, xlab="Taxon abundance",ylab="", col=grey.colors(1,0.7), pch=16)
  mtext("Mass", side=3, cex=1.5)
  if(pval4<0.05 & mod.type[a]=="lm") lines(exp(ab.seq),predict(mod8,data.frame(abun=scale(ab.seq))),lwd=3,col=v.col[3])
 
  mtext(SES4, side = 3, adj = 0.975, line = -6.5, cex = text.cex, font=2)
  mtext(ppval4, side = 3, adj = 0.975, line = -7.5, cex = text.cex, font=2)
  mtext(fe2.R2, side = 3, adj = 0.975, line = -8.5, cex = text.cex, font=2)
  

dev.off()


