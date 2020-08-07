
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

###### COMBO #######
flock_combo <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/FLOCK_percentages_combo.xlsx",sheet = 1)
flock_combo <- as.data.frame(unclass(flock_combo))
flock_combo$Training <- factor(flock_combo$Training,levels=c("Pre","Post"))
flock_cbc <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/FLOCK_percentages_combo.xlsx",sheet = 3)
flock_fin <- flock_combo[,1:6]
for(i in 1:25){
  flock_fin[,i+6] <- ((flock_combo[,i+6]*flock_combo$GR*flock_cbc[,4]*10)/flock_combo$Leukocytes)
  colnames(flock_fin)[i+6] <- paste0("Pop",i)
  res <- paste0("Pop",i)
  data <- flock_fin[,c(1:6,i+6)]
  flock_combo_1 <- dcast(data,Subject+Condition~Timepoint)
  fc <- mutate(flock_combo_1,Peak_Pre=flock_combo_1[,4]-flock_combo_1[,3],
               Recovery_Pre=flock_combo_1[,5]-flock_combo_1[,3],
               Baseline_Post=flock_combo_1[,6]-flock_combo_1[,3],
               Peak_Post=flock_combo_1[,7]-flock_combo_1[,3],
               Recovery_Post=flock_combo_1[,8]-flock_combo_1[,3])[,-c(3:8)]
  fc_b <- mutate(flock_combo_1,Baseline_Pre=flock_combo_1[,3]-flock_combo_1[,3],
               Peak_Pre=flock_combo_1[,4]-flock_combo_1[,3],
               Recovery_Pre=flock_combo_1[,5]-flock_combo_1[,3],
               Baseline_Post=flock_combo_1[,6]-flock_combo_1[,3],
               Peak_Post=flock_combo_1[,7]-flock_combo_1[,3],
               Recovery_Post=flock_combo_1[,8]-flock_combo_1[,3])[,-c(3:8)]
  fc_melt <- melt(fc,na.rm=TRUE)
  fc_b_melt <- melt(fc_b,na.rm=TRUE)
  colnames(fc_melt)[3:4] <- c("Time_point","FC_difference")
  colnames(fc_b_melt)[3:4] <- c("Time_point","FC_difference")
  fc_melt$Training <- fc_melt$Time_point
  fc_melt <- fc_melt %>% separate(Training,c("a","Training"))
  fc_melt <- fc_melt[,-6]
  val <- min(fc_melt$FC_difference)
  if(val < 0){
    fc_melt$FC_difference <- fc_melt$FC_difference+abs(val)
    fc_b_melt$FC_difference <- fc_b_melt$FC_difference+abs(val)
  }
  
  fc_melt$Sample <- as.factor(paste0(fc_melt$Subject,"_",fc_melt$Time_point))
  formula <- fc_melt[,4] ~ Condition*Time_point+(1|Subject)
  model_2 <- glmer.nb(formula,data=fc_melt,family="negative.binomial",
                     control=glmerControl(optimizer="Nelder_Mead", tolPwrss=1e-2,optCtrl = list(maxfun = 100000)))
  ls1 <- emmeans(model_2,~Condition|Time_point)
  ls2 <- emmeans(model_2,~Condition|Time_point,type="response")
  test1 <- contrast(ls1,interaction = "trt.vs.ctrl1",adjust="none")
  if(i==1){
    pval <- cbind(as.data.frame(summary(test1)[1]),as.data.frame(summary(test1)[2]),as.data.frame(summary(test1)[3]),as.data.frame(summary(test1)[7]))
  }else{
    pval <- cbind(pval,as.data.frame(summary(test1)[7]))
  }
  dat <- melt(data.frame("Condition"=c("Asthma","Healthy"),"Baseline_Pre"=c(abs(val),abs(val)),
                    "Peak_Pre"=c(as.data.frame(ls2)[c(1,2),3])))
  quartz()
  #p <- emmip(model_2,Condition~Time_point,type="response",CI=T,cov.reduce = range)+scale_color_manual(values=c("gray33", "#56B4E9"))+
  #ggtitle(res)+theme_bw()+labs(title=paste0("Mean and CI across timepoints for Population ",i))+xlab("Time point")+ylab("Normalized counts")+
    #geom_point(aes(x="Baseline_Pre",y=abs(val)),colour="black")+
    #scale_x_discrete(limits=c("Baseline_Pre","Peak_Pre","Recovery_Pre","Baseline_Post","Peak_Post","Recovery_Post"))
  emmip(model_2,Condition~Time_point,type="response",CI=T,cov.reduce = range,lyt=c(1,2))+scale_linetype_manual(values=c("solid", "dashed"))+
  ggtitle(res)+theme_bw()+labs(title=paste0("Mean and CI across timepoints for Population ",i))+xlab("Time point")+ylab("Normalized counts")+
  geom_point(aes(x="Baseline_Pre",y=abs(val)),colour="black",size=0.9)+ylim(0,1300)+
  scale_x_discrete(limits=c("Baseline_Pre","Peak_Pre","Recovery_Pre","Baseline_Post","Peak_Post","Recovery_Post"))
  
  
  quartz()
  ggplot(data = fc_melt)+geom_boxplot(aes(x=Time_point,y=FC_difference,color=Condition))+
  ggtitle(paste0("Healthy vs asthmatic across time points for Population ",i))+
  xlab("Time point")+ylab("Normalized counts")+theme_grey()
  
  d2 <- cbind(fc_melt,pred=as.data.frame(fitted(model_2,type="response")))
  colnames(d2)[6] <- "Pred"
  #quartz()
  q <- ggplot(data=d2)+
    facet_wrap(~Subject)+geom_smooth(aes(x=Time_point,y=d2$Pred),method="lm",alpha=0.3)+
    geom_point(aes(x=Time_point,y=FC_difference))+
    theme_gray()+
    geom_line(aes(x=Time_point,y=d2$Pred,group=Subject))+
    ylab("Normalized counts")+xlab("Time point")+ggtitle(paste0("Line plots for pop",i))
  print(p)
  #print(q)
}
 write.xlsx(pval[,-3],"/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo_Pvals_02_20.xlsx")
 corr_counts_pvals <- apply(pval[,-c(1:3)],2,p.adjust,"BH")


 ############### CD193 PANEL #############
 ##Reading the percentages and normalize them with respect to 
 ##absolute granulocyte count which is a sum of eosinophil, basophil and neutrophil
 ##Populations to keep 2,7,12,13,16 and 18
 flock_cd193 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/CD193/CD193_FLOCK.xlsx",sheet=1)
 flock_cbc_193 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/CD193/CD193_FLOCK.xlsx",sheet = 3) 
 flock_fin_193 <- flock_cd193[,1:6]
 ls <- list(list())
 pop <- c(2,7,12,13,16,18)
 for (j in 1:18) {
   ii <- is.element(j,pop)
   if(ii == TRUE){
     i <- match(j,pop)
     flock_fin_193[,i+6] <- (flock_cd193[,i+6]*flock_cd193$GRpos*flock_cbc_193[,8]*10/flock_cd193$Leukocytes)
     colnames(flock_fin_193)[i+6] <- paste0("Pop",pop[i])
     res <- paste0("Pop",pop[i])
     data <- flock_fin_193[,c(1:6,i+6)]
     flock_cd193_1 <- dcast(data,Subject+Condition~Timepoint)
     fc <- mutate(flock_cd193_1,Peak_Pre=flock_cd193_1[,4]-flock_cd193_1[,3],
                  Recovery_Pre=flock_cd193_1[,5]-flock_cd193_1[,3],
                  Baseline_Post=flock_cd193_1[,6]-flock_cd193_1[,3],
                  Peak_Post=flock_cd193_1[,7]-flock_cd193_1[,3],
                  Recovery_Post=flock_cd193_1[,8]-flock_cd193_1[,3])[,-c(3:8)]
     fc_b <- mutate(flock_cd193_1,Baseline_Pre=flock_cd193_1[,3]-flock_cd193_1[,3],
                    Peak_Pre=flock_cd193_1[,4]-flock_cd193_1[,3],
                    Recovery_Pre=flock_cd193_1[,5]-flock_cd193_1[,3],
                    Baseline_Post=flock_cd193_1[,6]-flock_cd193_1[,3],
                    Peak_Post=flock_cd193_1[,7]-flock_cd193_1[,3],
                    Recovery_Post=flock_cd193_1[,8]-flock_cd193_1[,3])[,-c(3:8)]
     fc_melt <- melt(fc,na.rm=TRUE)
     fc_b_melt <- melt(fc_b,na.rm=TRUE)
     colnames(fc_melt)[3:4] <- c("Time_point","FC_difference")
     colnames(fc_b_melt)[3:4] <- c("Time_point","FC_difference")
     val <- min(fc_melt$FC_difference)
     if(val < 0){
       fc_melt$FC_difference <- fc_melt$FC_difference+abs(val)
       fc_b_melt$FC_difference <- fc_b_melt$FC_difference+abs(val)
     }
     fc_melt$Sample <- as.factor(paste0(fc_melt$Subject,"_",fc_melt$Time_point))
     formula <- fc_melt[,4] ~ Condition*Time_point+(1|Subject)
     model_2 <- glmer.nb(formula,data=fc_melt,family="negative.binomial",
                         control=glmerControl(optimizer="Nelder_Mead", tolPwrss=1e-2,optCtrl = list(maxfun = 100000)))
     ls1 <- emmeans(model_2,~Condition|Time_point)
     ls2 <- emmeans(model_2,~Condition|Time_point,type="response")
     test1 <- contrast(ls1,interaction = "trt.vs.ctrl1",adjust="none")
     if(i==1){
       pval <- cbind(as.data.frame(summary(test1)[1]),as.data.frame(summary(test1)[2]),as.data.frame(summary(test1)[3]),as.data.frame(summary(test1)[7]))
     }else{
       pval <- cbind(pval,as.data.frame(summary(test1)[7]))
     }
     dat <- melt(data.frame("Condition"=c("Asthma","Healthy"),"Baseline_Pre"=c(abs(val),abs(val)),
                            "Peak_Pre"=c(as.data.frame(ls2)[c(1,2),3])))
     quartz()
     p <- emmip(model_2,Condition~Time_point,type="response",CI=T,cov.reduce = range)+scale_color_manual(values=c("gray33", "#56B4E9"))+
       ggtitle(res)+theme_bw()+labs(title=paste0("Mean and CI across timepoints for Population ",j," in CD193 panel"))+xlab("Time point")+ylab("Normalized counts")+
       geom_point(aes(x="Baseline_Pre",y=abs(val)),colour="black")+
       scale_x_discrete(limits=c("Baseline_Pre","Peak_Pre","Recovery_Pre","Baseline_Post","Peak_Post","Recovery_Post"))
     
     d2 <- cbind(fc_melt,pred=as.data.frame(fitted(model_2,type="response")))
     colnames(d2)[6] <- "Pred"
     #quartz()
     q <- ggplot(data=d2)+
       facet_wrap(~Subject)+geom_smooth(aes(x=Time_point,y=d2$Pred),method="lm",alpha=0.3)+
       geom_point(aes(x=Time_point,y=FC_difference))+
       theme_gray()+
       geom_line(aes(x=Time_point,y=d2$Pred,group=Subject))+
       ylab("Normalized counts")+xlab("Time point")+ggtitle(paste0("Line plots for pop",j))
     print(p)
     #print(q)
   }
 }
 
 write.xlsx(pval[,-3],"/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/CD193/Baseline_diff/",rownames=T)

 
 ### CD203 ###
 ############### CD203 PANEL #############
 ##Reading the percentages and normalize them with respect to 
 ##absolute granulocyte count which is a sum of eosinophil, basophil and neutrophil
 ##Populations to keep 2,7,12,13,16 and 18
 flock_cd203 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/CD203/CD203_FLOCK.xlsx",sheet=1)
 flock_cbc_203 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/CD203/CD203_FLOCK.xlsx",sheet = 3) 
 flock_fin_203 <- flock_cd203[,1:6]
 ls <- list(list())
 pop2 <- c(2,6,8,12,13,17)
 for (j in 1:17) {
   ii <- is.element(j,pop2)
   if(ii == TRUE){
     i <- match(j,pop2)
     flock_fin_203[,i+6] <- (flock_cd203[,i+6]*flock_cbc_203[,8]*flock_cd203$GRpos*10/flock_cd203$Leukocytes)
     colnames(flock_fin_203)[i+6] <- paste0("Pop",pop2[i])
     res <- paste0("Pop",pop2[i])
     data <- flock_fin_203[,c(1:6,i+6)]
     flock_cd203_1 <- dcast(data,Subject+Condition~Timepoint)
     fc <- mutate(flock_cd203_1,Peak_Pre=flock_cd203_1[,4]-flock_cd203_1[,3],
                  Recovery_Pre=flock_cd203_1[,5]-flock_cd203_1[,3],
                  Baseline_Post=flock_cd203_1[,6]-flock_cd203_1[,3],
                  Peak_Post=flock_cd203_1[,7]-flock_cd203_1[,3],
                  Recovery_Post=flock_cd203_1[,8]-flock_cd203_1[,3])[,-c(3:8)]
     fc_b <- mutate(flock_cd203_1,Baseline_Pre=flock_cd203_1[,3]-flock_cd203_1[,3],
                    Peak_Pre=flock_cd203_1[,4]-flock_cd203_1[,3],
                    Recovery_Pre=flock_cd203_1[,5]-flock_cd203_1[,3],
                    Baseline_Post=flock_cd203_1[,6]-flock_cd203_1[,3],
                    Peak_Post=flock_cd203_1[,7]-flock_cd203_1[,3],
                    Recovery_Post=flock_cd203_1[,8]-flock_cd203_1[,3])[,-c(3:8)]
     fc_melt <- melt(fc,na.rm=TRUE)
     fc_b_melt <- melt(fc_b,na.rm=TRUE)
     colnames(fc_melt)[3:4] <- c("Time_point","FC_difference")
     colnames(fc_b_melt)[3:4] <- c("Time_point","FC_difference")
     val <- min(fc_melt$FC_difference)
     if(val < 0){
       fc_melt$FC_difference <- fc_melt$FC_difference+abs(val)
       fc_b_melt$FC_difference <- fc_b_melt$FC_difference+abs(val)
     }
     fc_melt$Sample <- as.factor(paste0(fc_melt$Subject,"_",fc_melt$Time_point))
     formula <- fc_melt[,4] ~ Condition*Time_point+(1|Subject)
     model_2 <- glmer.nb(formula,data=fc_melt,family="negative.binomial",
                         control=glmerControl(optimizer="Nelder_Mead", tolPwrss=1e-2,optCtrl = list(maxfun = 100000)))
     ls1 <- emmeans(model_2,~Condition*Time_point)
     ls2 <- emmeans(model_2,~Condition|Time_point,type="response")
     test1 <- contrast(ls1,interaction = "trt.vs.ctrl1",adjust="none")
     if(i==1){
       pval <- cbind(as.data.frame(summary(test1)[1]),as.data.frame(summary(test1)[2]),as.data.frame(summary(test1)[3]),as.data.frame(summary(test1)[7]))
     }else{
       pval <- cbind(pval,as.data.frame(summary(test1)[7]))
     }
     dat <- melt(data.frame("Condition"=c("Asthma","Healthy"),"Baseline_Pre"=c(abs(val),abs(val)),
                            "Peak_Pre"=c(as.data.frame(ls2)[c(1,2),3])))
     quartz()
     p <- emmip(model_2,Condition~Time_point,type="response",CI=T,cov.reduce = range)+scale_color_manual(values=c("gray33", "#56B4E9"))+
       ggtitle(res)+theme_bw()+labs(title=paste0("Mean and CI across timepoints for Population ",j," in CD203 panel"))+xlab("Time point")+ylab("Normalized counts")+
       geom_point(aes(x="Baseline_Pre",y=abs(val)),colour="black")+
       scale_x_discrete(limits=c("Baseline_Pre","Peak_Pre","Recovery_Pre","Baseline_Post","Peak_Post","Recovery_Post"))
     
     d2 <- cbind(fc_melt,pred=as.data.frame(fitted(model_2,type="response")))
     colnames(d2)[6] <- "Pred"
     #quartz()
     q <- ggplot(data=d2)+
       facet_wrap(~Subject)+geom_smooth(aes(x=Time_point,y=d2$Pred),method="lm",alpha=0.3)+
       geom_point(aes(x=Time_point,y=FC_difference))+
       theme_gray()+
       geom_line(aes(x=Time_point,y=d2$Pred,group=Subject))+
       ylab("Normalized counts")+xlab("Time point")+ggtitle(paste0("Line plots for pop",j))
     print(p)
     #print(q)
   }
 }
 write.xlsx(pval[,-3],"/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/CD203/Baseline_diff/",rownames=T)
  