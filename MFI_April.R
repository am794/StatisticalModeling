######################
# MFI 25 populations #
# Date: 04/04/2019   #
######################

library(openxlsx)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(reshape)
library(reshape2)

mfi_combo <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/MFI/FLOCK_MFI_all25.xlsx",rowNames=TRUE,colNames = TRUE)
mfi_combo$TimePoint <- interaction(mfi_combo$Time_point,mfi_combo$Training)
colnames(mfi_combo)[6:30] <- paste0("Population",1:25) 

data <- mfi_combo[,c(1:5,31,i+5)]
res <- paste0("Population",i)

mixed <- function(data,res){
  data <- as.data.frame(unclass(data))
  data$Training <- factor(data$Training,levels=c("Pre","Post"))
  data$TimePoint <- factor(data$TimePoint,levels=c("Baseline.Pre","Peak.Pre","Recovery.Pre","Baseline.Post","Peak.Post","Recovery.Post"))
  data$Time_point.cont <- as.numeric(as.factor(data$Time_point))
  data$Training.cont <- as.numeric(as.factor(data$Training))
  data$Condition.cont <- as.factor(as.factor(data$Condition))
  formula_3w <- data[,res] ~ TimePoint+Condition+TimePoint*Condition+(1|Subject)
  model_3w <- lmer(formula_3w,data=data)
  
  #Post-hoc tests
  
  ls6 <- lsmeans(model_3w,~Condition|TimePoint)
  test6 <- emmeans::contrast(ls6,interaction = "trt.vs.ctrl1",adjust="none")
  
  ls7 <- lsmeans(model_3w,~ Condition*TimePoint)
  test7 <- emmeans::contrast(ls7,interaction = "trt.vs.ctrl1",adjust="none")
  
  ls8 <- lsmeans(model_3w,~TimePoint|Condition)
  test8 <- emmeans::contrast(ls8,interaction = "trt.vs.ctrl", adjust="none")
  
  ls9 <- lsmeans(model_3w,~TimePoint*Condition)
  test9 <- emmeans::contrast(ls9,interaction = "trt.vs.ctrl", adjust="none")
  
  tab6 <- cbind(as.data.frame(summary(test6)[1]),as.data.frame(summary(test6)[2]),as.data.frame(summary(test6)[7]))
  tab7 <- cbind(as.data.frame(summary(test7)[1]),as.data.frame(summary(test7)[2]),as.data.frame(summary(test7)[7]))
  tab8 <- cbind(as.data.frame(summary(test8)[1]),as.data.frame(summary(test8)[2]),as.data.frame(summary(test8)[7]))
  colnames(tab6)<-colnames(tab7)<-colnames(tab8)<-c("Contrast","Level2","p-value")
  p_summary <- rbind(tab6,tab7,tab8)
  return(p_summary)
}
quartz()
bb <- cbind(data,pred=fitted(model_3w))
ggplot(data=bb)+facet_wrap(~Subject)+
  geom_point(aes(x=TimePoint,y=data[,7]))+
  geom_line(aes(x=TimePoint,y=pred,group=1))

quartz()
ggplot(data=bb)+geom_boxplot(aes(x=TimePoint,y=bb[,7],color=Condition))+theme_bw()+
  ggtitle(paste0("Boxplot of MFI between asthma and control in ",res))+ylab("MFI")+
  xlab("Time Point")

mfi_pvals <- c()
for(i in 1:25){
  if(i == 1) { 
    mfi_pvals <- assign(paste0("mm_",i),mixed(mfi_combo[,c(1:5,31,i+5)],paste0("Population",i)))
  } else {
    a <- assign(paste0("mm_",i),mixed(mfi_combo[,c(1:5,31,i+5)],paste0("Population",i)))
    mfi_pvals <- cbind(mfi_pvals,a[,3])
  }
}
mfi_pvals <- mfi_pvals[c(1,4,9,14,19),]
colnames(mfi_pvals)[3:27] <- paste0("Pop",seq(1,25,1))
write.xlsx(mfi_pvals,"/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/MFI/MFI_pvals.xlsx",rowNames=TRUE,colNames=T)


###############Fold change
library(openxlsx)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(dplyr)

mfi_combo <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/MFI/FLOCK_MFI_all25.xlsx",rowNames=TRUE,colNames = TRUE)
mfi_combo$TimePoint <- interaction(mfi_combo$Time_point,mfi_combo$Training)
colnames(mfi_combo)[6:30] <- paste0("Population",1:25) 

data <- mfi_combo[,c(1:5,31,i+5)]
res <- paste0("Population",i)
####for baseline post/baseline pre
mixed2 <- function(data){
  data <- as.data.frame(unclass(data))
  data$Training <- factor(data$Training,levels=c("Pre","Post"))
  data$TimePoint <- factor(data$TimePoint,levels=c("Baseline.Pre","Peak.Pre","Recovery.Pre","Baseline.Post","Peak.Post","Recovery.Post"))
  data$Time_point.cont <- as.numeric(as.factor(data$Time_point))
  data$Training.cont <- as.numeric(as.factor(data$Training))
  data$Condition.cont <- as.factor(as.factor(data$Condition))
  data_cast <- dcast(data[-c(78,79),1:7],Subject+Condition ~ Time_point+Training)
  data_fc <- mutate(data_cast,Baseline.Pre.FC=data_cast[,3]/data_cast[,3],Peak.Pre.FC=data_cast[,5]/data_cast[,3],
                    Recovery.Pre.FC=data_cast[,7]/data_cast[,3],
                    Baseline.Post.FC=data_cast[,4]/data_cast[,3],Peak.Post.FC=data_cast[,6]/data_cast[,3],
                    Recovery.Post.FC=data_cast[,8]/data_cast[,3] )
  data_fc_melt <- melt(data_fc[,c(1:2,9:14)])
  data_fc_final <- separate(data_fc_melt,variable, c("TimePoint","Training"))
  formula_3w <- data_fc_final[,5] ~ Condition+TimePoint+Training+Condition*TimePoint*Training+(1|Subject)
  model_3w <- lmer(formula_3w,data=data_fc_final)
  #formula_2w <- data_fc_melt[,5] ~ Training+variable+Condition+Training*variable*Condition+(1|Subject)
  #Post-hoc tests
  formula_2w <- data_fc_final[,5] ~ Condition+TimePoint+Condition*TimePoint+(1|Subject)
  model_2w <- lmer(formula_2w,data=data_fc_final)
  
  ls6 <- lsmeans(model_3w,~Condition*Training|TimePoint)
  test6 <- emmeans::contrast(ls6,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  ls7 <- lsmeans(model_3w,~Condition|TimePoint*Training )
  test7 <- emmeans::contrast(ls7,interaction = "trt.vs.ctrl1",adjust="none")
  
  tab6 <- cbind(as.data.frame(summary(test6)[1]),as.data.frame(summary(test6)[2]),
                as.data.frame(summary(test6)[3]),as.data.frame(summary(test6)[8]))
  tab7 <- cbind(as.data.frame(summary(test7)[1]),as.data.frame(summary(test7)[2]),
                as.data.frame(summary(test7)[3]),as.data.frame(summary(test7)[8]))
  #tab8 <- cbind(as.data.frame(summary(test8)[1]),as.data.frame(summary(test8)[2]),as.data.frame(summary(test8)[7]))
  colnames(tab6)<- colnames(tab7) <- c("Contrast","Level2","Level3","p-value")
  p_summary <- rbind(tab6,tab7)
  return(p_summary)
}

mfi_pvals <- c()
for(i in 1:25){
  if(i == 1) {
    mfi_pvals <- assign(paste0("mm_",i),mixed2(mfi_combo[,c(1:5,31,i+5)]))
  } else {
    a <- assign(paste0("mm_",i),mixed2(mfi_combo[,c(1:5,31,i+5)]))
    mfi_pvals <- cbind(mfi_pvals,a[,4])
  }
}
#mfi_pvals <- mfi_pvals[c(1,4,9,14,19),]
colnames(mfi_pvals)[4:28] <- paste0("Pop",seq(1,25,1))
write.xlsx(mfi_pvals,"/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/MFI/MFI_pvals_june2019.xlsx",rowNames=TRUE,colNames=T)

quartz()
bb <- cbind(data,pred=fitted(model_3w))
ggplot(data=bb)+facet_wrap(~Subject)+
  geom_point(aes(x=TimePoint,y=data[,7]))+
  geom_line(aes(x=TimePoint,y=pred,group=1))

dat <- na.omit(data_fc_final[which(data_fc_final$TimePoint == "Recovery"),])
quartz()
ggplot(data=dat)+geom_boxplot(aes(x=Training,y=value,color=Condition))+theme_bw()+
  ggtitle(paste0("Boxplot of MFI between asthma and control in ",res))+ylab("MFI")+
  xlab("Training")

dat2 <- na.omit(data_fc_final[which(data_fc_final$TimePoint == "Peak" & data_fc_final$Training == "Pre"),])
quartz()
ggplot(data=dat2)+geom_boxplot(aes(x=Condition,y=value))+theme_bw()+
  ggtitle(paste0("Boxplot of MFI between asthma and control in ",res))+ylab("MFI")+
  xlab("Condition")

##CD193

