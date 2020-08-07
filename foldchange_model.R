###################################################
# Function for GLMM neg binomial for cell counts  #
###################################################

##Loading the packages
library(openxlsx)
library(lme4)
library("nlme")
library(phia)
library("ggplot2")
library("reshape")
library(dplyr)
library(plyr)
library(emmeans)

##Combo data
flock_combo <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/FLOCK_percentages_combo.xlsx",sheet = 1)
flock_combo <- as.data.frame(unclass(flock_combo))
flock_combo$Training <- factor(flock_combo$Training,levels=c("Pre","Post"))
flock_cbc <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/FLOCK_percentages_combo.xlsx",sheet = 3)
flock_fin <- flock_combo[,1:6]
ls <- list(list())
for(i in 1:25){
  flock_fin[,i+6] <- ((flock_combo[,i+6]*flock_combo$GR*flock_cbc[,4]*10)/flock_combo$Leukocytes)
  colnames(flock_fin)[i+6] <- paste0("Pop",i)
  res <- paste0("Pop",i)
  data <- flock_fin[,c(1:6,i+6)]
  flock_com <- data[,c(1:3,7)]
  f <- interaction(data[,6],data[,5])
  flock_com <- cbind(flock_com,f)
  flock_com <- flock_com[,c(1,2,3,5,4)]
  colnames(flock_com)[4] <- "Time_point"
  final_flock_combo <- dcast(flock_com,Subject+Condition~Time_point)
  #infini <- which(final_flock_combo$baseline.Pre == 0 | final_flock_combo$baseline.Post == 0)
  #if(infini > 0){
    #if(final_flock_combo$baseline.Pre[infini] == 0){
      #final_flock_combo$baseline.Pre[infini] = 1
    #} else if(final_flock_combo$baseline.Post[infini] == 0){
      #final_flock_combo$baseline.Post[infini] = 1
    #}
  #}
  #fc <- mutate( final_flock_combo,Baseline.Pre.FC=final_flock_combo[,3]/ final_flock_combo[,3],
                #Peak.Pre.FC=final_flock_combo[,4]/ final_flock_combo[,3],
                #Recovery.Pre.FC=final_flock_combo[,5]/ final_flock_combo[,3],
                #Baseline.Post.FC=final_flock_combo[,6]/ final_flock_combo[,6],
                #Peak.Post.FC=final_flock_combo[,7]/ final_flock_combo[,6],
                #Recovery.Post.FC=final_flock_combo[,8]/ final_flock_combo[,6])[,c(1,2,9:14)]
  fc <- mutate( final_flock_combo,Peak.Pre.FC=final_flock_combo[,4]/ final_flock_combo[,3],Recovery.Pre.FC=final_flock_combo[,5]/ final_flock_combo[,3],
                Peak.Post.FC=final_flock_combo[,7]/ final_flock_combo[,6],
                Recovery.Post.FC=final_flock_combo[,8]/ final_flock_combo[,6])[,c(1,2,9:12)]
  fc_melt <- melt(fc,na.rm=TRUE)
  inf <- which(is.infinite(fc_melt$value) == "TRUE")
  if(length(inf) > 0){
    for(j in inf){
    fc_melt <- fc_melt[-inf,]
    } 
  }
  colnames(fc_melt) <- c("Subject","Condition","Time_point","Fold_change")
  
  #quartz()
  #p <- ggplot(data = fc_melt)+geom_boxplot(aes(x=Time_point,y=Fold_change,color=Condition))+
    #ggtitle(paste0("Healhty vs asthmatic across time points for Population ",i))+
    #xlab("Time point")+ylab("Fold change")+theme_grey()+ylim(0,6)
  #print(p)
    fcm <- fc_melt
    fcm <- fcm %>% separate(Time_point,c("Timepoint","Training","rest"))
    fcm$Training <- factor(fcm$Training,levels=c("Pre","Post"))
  #quartz()
  #p<-interaction.ABC.plot(Fold_change,Timepoint,Condition,Training,data = fcm,fun="mean",
                       #title = paste0("Foldchange across timepoints for Population ",i),
                       #xlab="Time_point",ylab="Mean of Foldchange")+theme_grey()
  #print(p)
  formula_1 <- fcm[,"Fold_change"] ~ Condition*Timepoint*Training+(1+Training|Subject)
  #formula_2 <- fc_melt[,"Fold_change"] ~ Condition*Time_point+(1+Time_point|Subject)
  model_1 <- glmer.nb(formula_1,data=fcm,family="negative.binomial",
                       control=glmerControl(optimizer="Nelder_Mead", tolPwrss=1e-2,optCtrl = list(maxfun = 100000)))
  #model_2 <- glmer.nb(formula_2,data=fc_melt,family="negative.binomial")
  ls1 <- emmeans(model_1,~Condition|Timepoint*Training)
  test1 <- contrast(ls1,interaction = "trt.vs.ctrl1",adjust="none")
  if(i==1){
 pval <- cbind(as.data.frame(summary(test1)[1]),as.data.frame(summary(test1)[2]),as.data.frame(summary(test1)[3]),as.data.frame(summary(test1)[8]))
  }else{
   pval <- cbind(pval,as.data.frame(summary(test1)[8]))
 }
}
quartz()
ggplot(data=cbind(fcm,pred=as.data.frame(fitted(model_1,type="response"))))+
  facet_wrap(~Subject)+geom_smooth(aes(x=Timepoint,y=pred),method="lm",alpha=0.3)+
  geom_point(aes(x=Timepoint,y=Fold_change,colour=Training))+
  scale_color_manual(values=c( "#E69F00", "#56B4E9"))+theme_gray()+
  geom_line(aes(x=Timepoint,y=pred,colour=Training,group=Training))+ylim(0,10)+
  ylab("Normalized counts")+xlab("Time point")+ggtitle(paste0("Line plots for pop",i))

  
}
colnames(pval)[4:28] <- paste0("Pop_",seq(1,25,1))

#######################
# Model using weights #
#######################
for(i in 1:25){
  flock_fin[,i+6] <- ((flock_combo[,i+6]*flock_combo$GR*flock_cbc[,4]*10)/flock_combo$Leukocytes)
  colnames(flock_fin)[i+6] <- paste0("Pop",i)
  res <- paste0("Pop",i)
  data <- flock_fin[,c(1:6,i+6)]
  flock_combo_1 <- dcast(data,Subject+Condition+Training~Time_point)
  fc1 <- flock_combo_1[,1:5]
  fc1$Peak.Recovery <- rep("Peak",nrow(fc1))
  fc2 <- flock_combo_1[,c(1:4,6)]
  fc2$Peak.Recovery <- rep("Recovery",nrow(fc2))
  colnames(fc1)[5] <- colnames(fc2)[5] <- "FC"
  flock_combo_2 <- rbind(fc1,fc2)
  fl_c <- as.data.frame(unclass(flock_combo_2))[-c(47,51,68),]
  fl_c$Foldchange <- fl_c[,5]/fl_c[,4]
  inf <- which(fl_c$Foldchange == "Inf" | fl_c$Foldchange == "NaN" )
  if(length(inf) > 0){
    fl_c <- fl_c[-inf,]
  }
  Timepoint <- interaction(fl_c$Peak.Recovery,fl_c$Training)
  fl_c <- cbind(fl_c,Timepoint)
  formula_1 <- fl_c[,"FC"]/fl_c[,"baseline"] ~ Condition*Timepoint+(1|Subject)
  #formula_2 <- fc_melt[,"Fold_change"] ~ Condition*Time_point+(1+Time_point|Subject)
  model_2 <- glmer.nb(formula_1,data=fl_c,family="negative.binomial",weights=baseline,
                      control=glmerControl(optimizer="Nelder_Mead", tolPwrss=1e-2,optCtrl = list(maxfun = 100000)))
  #model_2 <- glmer.nb(formula_2,data=fc_melt,family="negative.binomial")
  ls1 <- emmeans(model_1,~Condition|Peak.Recovery*Training)
  
  test1 <- contrast(ls1,interaction = "trt.vs.ctrl1",adjust="none")
  #quartz()
  #emmip(model_1,Condition~Peak.Recovery|Training,type="response",CI=T,cov.reduce = range)+scale_color_manual(values=c("gray33", "#56B4E9"))+
    #ggtitle(res)+theme_bw()
  
  if(i==1){
    pval <- cbind(as.data.frame(summary(test1)[1]),as.data.frame(summary(test1)[2]),as.data.frame(summary(test1)[3]),as.data.frame(summary(test1)[8]))
  }else{
    pval <- cbind(pval,as.data.frame(summary(test1)[8]))
  }
}

quartz()
ggplot(data=cbind(fl_c,pred=as.data.frame(fitted(model_1,type="response"))))+
  facet_wrap(~Subject)+geom_smooth(aes(x=Peak.Recovery,y=pred),method="lm",alpha=0.3)+
  geom_point(aes(x=Peak.Recovery,y=fl_c$Foldchange,colour=Training))+
  scale_color_manual(values=c( "#E69F00", "#56B4E9"))+theme_gray()+
  geom_line(aes(x=Peak.Recovery,y=pred,colour=Training,group=Training))+ylim(0,10)+
  ylab("Normalized counts")+xlab("Time point")+ggtitle(paste0("Line plots for pop",i))


