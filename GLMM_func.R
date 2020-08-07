###################################################
# Function for GLMM neg binomial for cell counts  #
# 2w:2 way interactions for two way comparisons   #
# 3w:3 way interactions for three way comparisons #
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

############ COMBO PANEL #############
##Reading the percentages and normalize them with respect to CBC data
flock <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/FLOCK_percentages_combo.xlsx",sheet = 1)
flock_cbc <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/FLOCK_percentages_combo.xlsx",sheet = 3)
flock_fin <- flock[,1:6]
ls <- list(list())

for(i in 1:25){
  flock_fin[,i+6] <- (flock[,i+6]*flock_cbc[,4]*10)
  colnames(flock_fin)[i+6] <- paste0("Pop",i)
  res <- paste0("Pop",i)
  data <- flock_fin[,c(1:6,i+6)]
  #ls[[i]] <- gen_MM(data,res)
  if(i == 1){
    glmm_pvals <- assign(paste0("p_",i),gen_MM(data,res))
  } else if(i > 1 && i < 26){
    a <- assign(paste0("p_",i),gen_MM(data,res))
    glmm_pvals <- cbind(glmm_pvals,a[,4])
  }
}  
colnames(glmm_pvals)[4:28] <- paste0("Pop",seq(1,25,1))
write.xlsx(glmm_pvals,"/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/p-vals_combo_rs_percs.xlsx",rownames=T)

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
  flock_fin_193[,i+6] <- (flock_cd193[,i+6]*flock_cbc_193[,8]*10)
  colnames(flock_fin_193)[i+6] <- paste0("Pop",pop[i])
  res <- paste0("Pop",pop[i])
  data <- flock_fin_193[,c(1:6,i+6)]
  if(i == 1){
    glmm_pvals_cd193 <- assign(paste0("p_",pop[i]),gen_MM(data,res))
  } else {
    a <- assign(paste0("p_",pop[i]),gen_MM(data,res))
    glmm_pvals_cd193 <- cbind(glmm_pvals_cd193,a[,4])
  }
}
}
colnames(glmm_pvals_cd193)[4:9] <- paste0("Pop",pop)
write.xlsx(glmm_pvals_cd193,"/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/p-vals_cd193_rs2_percs.xlsx",rownames=T)

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
    flock_fin_203[,i+6] <- (flock_cd203[,i+6]*flock_cbc_203[,8]*10)
    colnames(flock_fin_203)[i+6] <- paste0("Pop",pop2[i])
    res <- paste0("Pop",pop2[i])
    data <- flock_fin_203[,c(1:6,i+6)]
    if(i == 1){
      glmm_pvals_cd203 <- assign(paste0("p_",pop2[i]),gen_MM(data,res))
    } else {
      a <- assign(paste0("p_",pop2[i]),gen_MM(data,res))
      glmm_pvals_cd203 <- cbind(glmm_pvals_cd203,a[,4])
    }
  }
}
colnames(glmm_pvals_cd203)[4:9] <- paste0("Pop",pop2)
write.xlsx(glmm_pvals_cd203,"/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/p-vals_cd203_rs2_percs.xlsx",rownames=T)

Training_Pre <- subset(data,data$Training=="Pre")
Training_Post <- subset(data,data$Training=="Post")

quartz()
with(Training_Pre ,{ 
  par(mfrow = c(1,2)) 
  interaction.plot(Time_point,factor(Condition),Pop18,type = "b",legend="F",col=c(4,6),trace.label = "Pre") 
  
  #legend(3,300, c("Asthmatic","Control"),col = c(4,6),text.col = "green4", lty = c(2,1,3), pch = c(4, 6),merge = TRUE, bty = "n")
  })
with(Training_Post ,{ 
  #par(mfrow = c(2,2)) 
  interaction.plot(Time_point,factor(Condition),Pop18,type = "b",legend="F",col=c(4,6),trace.label = "Post") 
  
  #legend(3,300, c("Asthmatic","Control"),col = c(4,6),text.col = "green4", lty = c(2,1,3), pch = c(4, 6),merge = TRUE, bty = "n")
})
quartz()
interaction.ABC.plot(Pop9,Training,Condition,Time_point,data = data,fun="mean",title = paste0("Baseline vs recovery mean for ",res),
                     xlab="Time_point",ylab="Mean of Normalized counts")+
  scale_color_manual(values=c( "#E69F00", "#56B4E9","#773F07"))

quartz()
ggplot(data = data)+geom_boxplot(aes(x=Training,y=Pop17,color=Condition))+facet_grid(~Time_point)+
  ggtitle(paste0("Healhty vs asthmatic across time points for ",res))+
  xlab("Time point")+ylab("Normalized counts")

+
  scale_color_manual(values=c( "#E69F00", "#56B4E9","#773F07"))

quartz()
interaction.ABC.plot(Pop9,Time_point,Condition,Training,data = data,fun="mean",title = paste0("Baseline vs recovery mean for ",res),
                     xlab="Time_point",ylab="Mean of Normalized counts")

quartz()
ggplot(data = data)+geom_boxplot(aes(x=Time_point,y=Pop12,color=Condition))+facet_grid(~Training)+
  ggtitle(paste0("Healhty vs asthmatic across time points for ",res))+
  xlab("Time point")+ylab("Normalized counts")

quartz()
emmip(model_8w,~ Condition+Time_point|Training,type="response",CI=T,cov.reduce = F)+scale_color_manual(values=c("gray33", "#56B4E9"))+
  ggtitle(res)+theme_bw()

##GLMM  re-normalized
gen_MM <- function(data,res){
  lst <- list(list())
  data <- as.data.frame(unclass(data))
  data$Training <- factor(data$Training,levels=c("Pre","Post"))
  data$Interaction <- as.numeric(as.factor(interaction(data$Time_point,data$Training)))
  data$Time_point.cont <- as.numeric(as.factor(data$Time_point))
  data$Training.cont <- as.numeric(as.factor(data$Training))
  data$Condition.cont <- as.factor(as.factor(data$Condition))
  #data <- data[,-4]
  #formula_6w <- data[,res] ~ Time_point+Training+Condition+Time_point*Training+Time_point*Condition+
   # Training*Condition+Training*Condition*Time_point+(1+Training|Subject)+(1|Subject)
  
  #formula_8w <- data[,res] ~ Time_point+Training+Condition+Time_point*Training+Time_point*Condition+
    #Training*Condition+Training*Condition*Time_point+(1+Interaction|Subject)+(1|Subject)
  
  formula_8w <- data[,res] ~ Time_point+Training+Condition+Time_point*Training+Time_point*Condition+
  Training*Condition+Training*Condition*Time_point+(1+Time_point.cont+Training.cont|Subject)
  
  formula_7w <- data[,res] ~ Time_point+Training+Condition+Time_point*Training+Time_point*Condition+
    Training*Condition+Training*Condition*Time_point+(1+Training.cont|Subject)
  
  #model_6w <- glmer.nb(formula_6w,data=data,family="negative.binomial",
                       #control=glmerControl(optimizer="Nelder_Mead", tolPwrss=1e-2,optCtrl = list(maxfun = 100000)))
  model_8w <- glmer.nb(formula_8w,data=data,family="negative.binomial",
                       control=glmerControl(optimizer="Nelder_Mead", tolPwrss=1e-2,optCtrl = list(maxfun = 100000)))
  
  #model_7w <- glmer.nb(formula_7w,data=data,family="negative.binomial",
                       #control=glmerControl(optimizer="Nelder_Mead", tolPwrss=1e-2,optCtrl = list(maxfun = 100000)))
  
  ls6 <- emmeans(model_8w,~Condition|Training*Time_point)
  test6 <- emmeans::contrast(ls6,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  
  ls7 <- lsmeans(model_8w,~Condition*Time_point|Training)
  test7 <- emmeans::contrast(ls7,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  
  ls8 <- lsmeans(model_8w,~Time_point|Condition*Training)
  test8 <- emmeans::contrast(ls8,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  
  ls9 <- lsmeans(model_8w,~Condition*Training|Time_point)
  test9 <- emmeans::contrast(ls9,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  
  ls10 <- lsmeans(model_8w,~Time_point*Condition|Training)
  test10 <- emmeans::contrast(ls10,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  
  ls11 <- lsmeans(model_8w,~Training*Condition|Time_point)
  test11 <- emmeans::contrast(ls11,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  
  l <- as.data.frame(rbind(ls7))
  ggplot(l, aes(x=Time_point,y=emmean,linetype=Condition,group=Condition))+ 
    geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE),colour="black", 
                  width=0.1,size=0.7,stat="identity", position = pd) +
    geom_line(position=pd,size=1) + facet_grid(~Training)+
    geom_point(position=pd,size=3,shape=21)+
    scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    ylim(c(0,15))+theme_bw()+
    theme(legend.justification=c(1,0),
          legend.position=c(1,0),
          plot.title = element_text(size = 12)) + ylim(0,7)+
    ylab("Estimated Marginal Means")+
    xlab("Time Point")+
    ggtitle(paste0("Mean-SE plot for normalized counts of ", res))
  
  
  ## Gender differences(for population 12 and 15 but no significant difference)
  #formula_7w <- data[,res] ~ Time_point+Training+Gender+Time_point*Training+Time_point*Gender+
    #Training*Gender+Training*Gender*Time_point+(1+Training|Subject)+(1|Subject)
  
  #model_7w <- glmer.nb(formula_7w,data=data,family="negative.binomial")
  #ls7 <- lsmeans(model_5w,~Gender|Time_point)
  #test7 <- emmeans::contrast(ls7,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  
  #formula_5w <- data[,res] ~ Time_point+Training+Gender+Time_point*Gender+
    #(1+Training|Subject)+(1|Subject)
  
  #model_5w <- glmer.nb(formula_5w,data=data,family="negative.binomial",
                       #control=glmerControl(optimizer="bobyqa", tolPwrss=1e-3,optCtrl = list(maxfun = 100000)))
  
  tab6 <- cbind(as.data.frame(summary(test6)[1]),as.data.frame(summary(test6)[2]),as.data.frame(summary(test6)[3]),summary(test6)[,8])
  tab7 <- cbind(as.data.frame(summary(test7)[1]),as.data.frame(summary(test7)[2]),as.data.frame(summary(test7)[3]),summary(test7)[,8])
  tab8 <- cbind(as.data.frame(summary(test8)[1]),as.data.frame(summary(test8)[2]),as.data.frame(summary(test8)[3]),summary(test8)[,8])
  tab9 <- cbind(as.data.frame(summary(test9)[1]),as.data.frame(summary(test9)[2]),as.data.frame(summary(test9)[3]),summary(test9)[,8])
  tab10 <- cbind(as.data.frame(summary(test10)[1]),as.data.frame(summary(test10)[2]),as.data.frame(summary(test10)[3]),summary(test10)[,8])
  tab11 <- cbind(as.data.frame(summary(test11)[1]),as.data.frame(summary(test11)[2]),as.data.frame(summary(test11)[3]),summary(test11)[,8])
 
   #tab7 <- cbind(as.data.frame(summary(test7)[1]),as.data.frame(summary(test7)[2]),as.data.frame(summary(test7)[3]),summary(test7)[,8])
  #colnames(tab6)<-colnames(tab7)<-c("Contrast","Level2","Level3","p-value")
  #p_summary <- rbind(tab6,tab7)
  #return(p_summary)
  colnames(tab6) <- colnames(tab7) <- colnames(tab8) <- colnames(tab9) <- colnames(tab10) <- colnames(tab11) <- c("Contrast","Level2","Level3","p-value")
  p_summary <- rbind(tab6,tab7,tab8,tab9,tab10,tab11)
  return(p_summary)
}
quartz()
emmip(model_8w,Condition~Time_point|Training,type="response",CI=F,cov.reduce = F)+scale_color_manual(values=c("gray33", "#56B4E9"))+
  ggtitle(res)+theme_bw()

##Function for generating negative binomial models with two way and three way interaction terms
gen_mixed_models <- function(data,res){
  lst <- list(list())
  data <- as.data.frame(unclass(data))
  data$Training <- factor(data$Training,levels=c("Pre","Post"))

  formula_2w <- data[,res] ~ Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
    Training*Gender+Training*Condition+Gender*Condition+(1|Subject)+(1+Training|Subject)
  
  formula_3w <- data[,res] ~ Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
    Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Condition*Gender+
    Training*Gender*Time_point+Condition*Gender*Time_point+(1+Training|Subject)+
    (1|Subject)
  
  formula_4w <- data[,res] ~ Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Condition+
    Training*Condition+Training*Condition*Time_point+(1+Training|Subject)
  
  formula_5w <- data[,res] ~ Time_point+Training+Condition+Time_point*Training+Time_point*Condition+
    Training*Condition+Training*Condition*Time_point+(1+Training|Subject)
  
  formula_6w <- data[,res] ~ Time_point+Training+Condition+Time_point*Training+Time_point*Condition+
    Training*Condition+Training*Condition*Time_point+(1+Training|Subject)+(1|Subject)
  
  model_2w <- glmer.nb(formula_2w,data=data,family="negative.binomial")
  model_3w <- glmer.nb(formula_3w,data=data,family="negative.binomial")
  model_4w <- glmer.nb(formula_4w,data=data,family="negative.binomial")
  model_5w <- glmer.nb(formula_5w,data=data,family="negative.binomial")
  model_6w <- glmer.nb(formula_6w,data=data,family="negative.binomial")
  
  #Post-hoc tests
  ls1 <- emmeans(model_2w,~Time_point|Training)
  test1 <- emmeans::contrast(ls1,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  ls2 <- emmeans(model_2w,~Training|Time_point)
  test2 <- emmeans::contrast(ls2,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  ls3 <- emmeans(model_2w,~Condition|Time_point)
  test3 <- emmeans::contrast(ls3,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  
  ls4 <- lsmeans(model_4w,~Time_point|Training*Condition)
  test4 <- emmeans::contrast(ls4,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  ls5 <- lsmeans(model_3w,~Gender|Training*Time_point)
  test5 <- emmeans::contrast(ls5,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  ls6 <- lsmeans(model_6w,~Condition|Training*Time_point)
  test6 <- emmeans::contrast(ls6,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  ls7 <- lsmeans(model_5w,~ Training|Time_point*Condition)
  test7 <- emmeans::contrast(ls7,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  
  ls8 <- lsmeans(model_3w,~ Condition*Training|Time_point)
  test8 <- emmeans::contrast(ls8,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  ls9 <- lsmeans(model_3w,~ Time_point*Training|Condition)
  test9 <- emmeans::contrast(ls9,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  
  tab1 <- cbind(as.data.frame(summary(test1)[,1]),as.data.frame(summary(test1)[2]),as.data.frame(noquote(rep('-',4))),summary(test1)[,7])
  tab2 <- cbind(as.data.frame(summary(test2)[,1]),as.data.frame(summary(test2)[2]),as.data.frame(noquote(rep('-',3))),summary(test2)[,7])
  tab3 <- cbind(as.data.frame(summary(test3)[,1]),as.data.frame(summary(test3)[2]),as.data.frame(noquote(rep('-',3))),summary(test3)[,7])
  tab4 <- cbind(as.data.frame(summary(test4)[1]),as.data.frame(summary(test4)[2]),as.data.frame(summary(test4)[3]),summary(test4)[,8])
  tab5 <- cbind(as.data.frame(summary(test5)[1]),as.data.frame(summary(test5)[2]),as.data.frame(summary(test5)[3]),summary(test5)[,8])
  tab6 <- cbind(as.data.frame(summary(test6)[1]),as.data.frame(summary(test6)[2]),as.data.frame(summary(test6)[3]),summary(test6)[,8])
  tab7 <- cbind(as.data.frame(summary(test7)[1]),as.data.frame(summary(test7)[2]),as.data.frame(summary(test7)[3]),summary(test7)[,8])
  tab8 <- cbind(as.data.frame(summary(test8)[1]),as.data.frame(summary(test8)[2]),as.data.frame(summary(test8)[3]),summary(test8)[,8])
  tab9 <- cbind(as.data.frame(summary(test9)[1]),as.data.frame(summary(test9)[2]),as.data.frame(summary(test9)[3]),summary(test9)[,8])
  colnames(tab1)<-colnames(tab2)<-colnames(tab3)<-colnames(tab4)<-colnames(tab5)<-colnames(tab6)<-colnames(tab7)<-colnames(tab8)<-colnames(tab9)<-c("Contrast","Level2","Level3","p-value")
  p_summary <- rbind(tab1,tab2,tab3,tab4,tab5,tab6,tab7,tab8,tab9)
  
  #Merge them into one list
  lst[[1]] <- p_summary
  lst[[2]] <- model_2w
  lst[[3]] <- model_3w
  lst[[4]] <- list(ls1,ls2,ls3,ls4,ls4,ls6,ls7,ls8,ls9)
  
  return(lst)
}

##Subject wise line plots for mixed effects model
quartz()
plot1 <- ggplot(data=cbind(data,pred=fitted(model_6w,type="response")))+
    facet_wrap(~Subject)+geom_smooth(aes(x=Time_point,y=pred),method="lm",alpha=0.3)+
    geom_point(aes(x=Time_point,y=data[,res],colour=Training))+
    scale_color_manual(values=c( "#E69F00", "#56B4E9"))+
    geom_line(aes(x=Time_point,y=pred,colour=Training,group=Training))+
    ylab("Normalized counts")+xlab("Time point")+ggtitle(paste0("Line plots for ",res))
  
##Point and errorbar plots
summ_res <- Rmisc::summarySE(data, measurevar=res, groupvars=c("Training","Time_point"))
pd <- position_dodge(0.3)
quartz()
pre_vs_post_mse <- 
    ggplot(summ_res, aes(x=Time_point,y=summ_res[,res],linetype=Training,group=Training))+ 
    geom_errorbar(aes(ymin=summ_res[,res]-se, ymax=summ_res[,res]+se),colour="black", 
                  width=0.1,size=0.7,stat="identity", position = pd) +
    geom_line(position=pd,size=1) +
    geom_point(position=pd,size=3,shape=21)+
    scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    ylim(c(0,15))+theme_bw()+
    theme(legend.justification=c(1,0),
          legend.position=c(1,0),
          plot.title = element_text(size = 12)) +
    ylab("MFI")+
    xlab("Time Point")+
    ggtitle(paste0("Mean-SE plot for MFI of ", res))
  
summ_res2 <- Rmisc::summarySE(data,measurevar= colnames(data)[7], groupvars=c("Training","Time_point","Condition"))
quartz()
ggplot(summ_res2, aes(x=Time_point,y=summ_res2[,7],linetype=Training,group=Training))+ 
    facet_grid(~Condition)+
    geom_errorbar(aes(ymin=summ_res2[,7]-se, ymax=summ_res2[,7]+se),colour="black", 
                  width=0.1,size=0.7,stat="identity", position = pd) +
    geom_line(position=pd,size=1) +
    geom_point(position=pd,size=3,shape=21)+
    scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    theme_bw()+
    theme(legend.justification=c(0.9,0.9),
          legend.position=c(0.95,0.95),
          plot.title = element_text(size = 12)) +
    ylab("Normalized counts")+
    xlab("Time Point")+
    ggtitle(paste0("Mean-SE plot for MFI of ", colnames(data2)[7]))
  
  
##Barplots
quartz()
pre_vs_post_msebar <- ggplot(summ_res, aes(x=Time_point,y=summ_res[,res]))+ 
    geom_errorbar(aes(ymin=summ_res[,res]-se, ymax=summ_res[,res]+se,group=Training), 
                  width=0.3,stat="identity", position = position_dodge(0.8)) +
    geom_bar(aes(fill=Training),position=position_dodge(0.8),stat="identity",width=0.5) +
    theme_bw()+
    theme(legend.justification=c(1,1),
          legend.position=c(0.98,0.98),
          plot.title = element_text(size = 12)) +
    ylab("Normalized counts")+
    xlab("Time Point")+
    ggtitle(paste0("Mean-SE plot for normalized counts of ", res))+
    scale_fill_manual(values = c("grey80", "grey30"))
  
  
quartz()
emmip(model_3w,Training~Time_point|Condition,type="response",CI=T)+scale_color_manual(values=c("gray33", "#56B4E9"))+
    theme_bw()
a <- regrid(ref_grid(model_3w),transform = TRUE)
emmip(a,Training~Time_point|Condition,CIs=TRUE,type="response")
plot(model_3w,comparisons=T)
quartz()
emmip(model_3w,Condition~Time_point+Training,type="response",CI=F,cov.reduce = range)+scale_color_manual(values=c("gray33", "#56B4E9"))+
    ggtitle(res)+theme_bw()
quartz()
emmip(model_2w,Training~Time_point,type="response",CI=F,cov.reduce = range)+scale_color_manual(values=c("gray33", "#56B4E9"))+
    ggtitle(res)+theme_bw()
  
  
##Point and errorbar plots
data2 <- cbind(data[,-7],pred=fitted(model_3w,type="response"))
PD <- colnames(data2)[7]
summ_res2 <- Rmisc::summarySE(data2,measurevar= colnames(data2)[7], groupvars=c("Training","Time_point","Condition"))
  
pd <- position_dodge(0.4)
quartz()
ggplot(summ_res2, aes(x=Time_point,y=summ_res2[,7],linetype=Training,group=Training))+ 
    facet_grid(~Condition)+
    geom_errorbar(aes(ymin=summ_res2[,7]-se, ymax=summ_res2[,7]+se),colour="black", 
                  width=0.1,size=0.7,stat="identity", position = pd) +
    geom_line(position=pd,size=1) +
    geom_point(position=pd,size=3,shape=21)+
    scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_bw()+
    theme(legend.justification=c(0.9,0.9),
          legend.position=c(0.95,0.95),
          plot.title = element_text(size = 12)) +
    ylab("Normalized counts")+ylim(c(-100,300))+
    xlab("Time Point")+
    ggtitle(paste0("Mean-SE plot for MFI of ", colnames(data2)[7]))
  
lst[[1]] <- p_summary
lst[[2]] <- model_2w
lst[[3]] <- model_3w
lst[[4]] <- list(ls1,ls2,ls3,ls4,ls4,ls6,ls7,ls8,ls9)
}

quartz()
print(pre_vs_post_mse)


##Reading the percentages and normalize them with respect to CBC data
flock <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/FLOCK_percentages_combo.xlsx",sheet = 1)
flock_cbc <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/FLOCK_percentages_combo.xlsx",sheet = 3)
flock_fin <- flock[,1:6]
data <- flock_fin[,c(1:6,i+6)]
res <- paste0("Pop",i)
for(i in 1:25){
  flock_fin[,i+6] <- (flock[,i+6]*flock_cbc[,4]*10)
}
colnames(flock_fin)[7:31] <- paste0("Pop",seq(1,25,1))
flock_fin <- cbind(flock_fin[,1:6],apply(flock_fin[,7:31],2,round))

#Mixed model analysis - calling the function 
#Ignore population 13 from the analysis
glmm_pvals <- c()
for(i in 1:25){
  if(i == 1){
    glmm_pvals <- assign(paste0("p_",i),gen_mixed_models(flock_fin[,c(1:6,i+6)],paste0("Pop",i))[[1]])
  } else if(i > 1 && i < 13){
    a <- assign(paste0("p_",i),gen_mixed_models(round(flock_fin[,c(1:6,i+6)]),paste0("Pop",i)))
    glmm_pvals <- cbind(glmm_pvals,a[[1]][,4])
  } else if(i > 13){
    a <- assign(paste0("p_",i),gen_mixed_models(flock_fin[,c(1:6,i+5)],paste0("Pop",i)))
    glmm_pvals <- cbind(glmm_pvals,a[[1]][,4])
  }
}
colnames(glmm_pvals)[4:15] <- paste0("Pop",seq(1,12,1))
colnames(glmm_pvals)[16:27] <- paste0("Pop",seq(14,25,1))

wb <- createWorkbook("UCI-pvalues.xlsx")
addWorksheet(wb,"Uncorrected-counts-pvaules")
writeData(wb,sheet = 1,glmm_pvals,rowNames = TRUE)
saveWorkbook(wb,"UCI-pvalues.xlsx")

corr_counts_pvals <- apply(glmm_pvals[,4:27],1,p.adjust,"BH")
write.xlsx(corr_counts_pvals,"counts-pvals-corrected.xlsx")

#Plots
i=2
data <- flock_fin[,c(1:6,i+6)]
res <- paste0("Pop",i)
data <- flock_fin_cd203[,c(1:6,i+6)]

####CD193 data
flock_cd193 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/CD193/CD193_FLOCK.xlsx",sheet=1)
flock_cbc_193 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/CD193/CD193_FLOCK.xlsx",sheet = 3)
norm <- 
flock_fin <- flock_cd193[,1:6]
glmm_pvals <- c()
for(i in 1:18){
  flock_fin[,i+6] <- (flock_cd193[,i+6]*flock_cd193$GRpos*flock_cbc_193[,6]*10)/flock_cd193$Leukocytes
  colnames(flock_fin)[i+6] <- paste0("Pop",i)
  #a <- assign(paste0("p_",i),gen_mixed_models(flock_fin[,c(1:6,i+6)],paste0("Pop",i)))
  #if(i==1){
    #glmm_pvals <- a[[1]]
  #}else{
    #glmm_pvals <- cbind(glmm_pvals,a[[1]][,4])
    #}
}
colnames(glmm_pvals)[4:21] <- paste0("Pop",seq(1,18,1))
corr_193_counts <- apply(glmm_pvals[,4:21],1,p.adjust,"BH")

write.xlsx(glmm_pvals,"FLOCK_CD193_pvals.xlsx")
write.xlsx(corr_193_counts,"flock_cd193_corrected.xlsx",rownames=T)

#####CD203 data
flock_cd203 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/CD203/CD203_FLOCK.xlsx",sheet = 1)
flock_cbc_203 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/CD203/CD203_FLOCK.xlsx",sheet = 3)
flock_fin_cd203 <- flock_cd203[,1:6]
glmm_pvals <- c()
for(i in 1:17){
  #flock_fin_cd203[,i+6] <- (flock_cd203[,i+6]*flock_cd203$GRpos*flock_cbc_203[,5]*10)/flock_cd203$Leukocytes
  flock_fin_cd203[,i+6] <- (flock_cd203[,i+6]*flock_cd203$GRpos)/flock_cd203$Leukocytes
  colnames(flock_fin_cd203)[i+6] <- paste0("Pop",i)
  #a <- assign(paste0("p_",i),gen_mixed_models(flock_fin_cd203[,c(1:6,i+6)],paste0("Pop",i)))
  #if(i==1){
    #glmm_pvals <- a[[1]]
  #}else{
    #glmm_pvals <- cbind(glmm_pvals,a[[1]][,4])
  #}
}

colnames(glmm_pvals)[4:20] <- paste0("Pop",seq(1,17,1))
corr_203_counts <- apply(glmm_pvals[,4:20],1,p.adjust,"BH")

write.xlsx(glmm_pvals,"FLOCK_CD203_pvals.xlsx")
write.xlsx(corr_203_counts,"flock_cd203_corrected.xlsx",rownames=T)

quartz()
ggplot(data=cbind(data,pred=fitted(model_3w,type="response")))+
  facet_wrap(~Subject)+geom_smooth(aes(x=Time_point,y=pred),method="lm",alpha=0.3)+
  geom_point(aes(x=Time_point,y=data[,res],colour=Training))+
  scale_color_manual(values=c( "#E69F00", "#56B4E9"))+
  geom_line(aes(x=Time_point,y=pred,colour=Training,group=Training))+
  ylab("Normalized counts")+xlab("Time point")+ggtitle(paste0("Line plots for ",res))




