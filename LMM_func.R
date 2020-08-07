###################################################
# Function for Linear mixed model(for MFI values) #
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
library(glmmADMB)

##Combo panel with 25 cell populations (MFI values)
mfi_25 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/MFI/FLOCK_MFI_all25.xlsx",sheet = 1)
colnames(mfi_25)[7:31] <- paste0("Pop",seq(1,25,1))

mfi_pvals <- c()

for(i in 1:25) {
  if(i == 1) { 
    mfi_pvals <- assign(paste0("mm_",i),mixed(mfi_25[,c(1:6,i+6)],paste0("Pop",i)))
  } else {
    a <- assign(paste0("mm_",i),mixed(mfi_25[,c(1:6,i+6)],paste0("Pop",i)))
    mfi_pvals <- cbind(mfi_pvals,a[,4])
  }
}
colnames(mfi_pvals)[4:28] <- paste0("Pop",seq(1,25,1))
write.xlsx(mfi_pvals,"/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/p-vals_rs_combo_mfi.xlsx",rownames=T)


wb <- createWorkbook("UCI")
addWorksheet(wb,"Uncorrected-mfi-pvaules")
writeData(wb,sheet = 1,mfi_pvals,rowNames = TRUE)
saveWorkbook(wb,"UCI-mfi-pvalues.xlsx")

corr_mfi_pvals <- apply(mfi_pvals[,4:28],1,p.adjust,"BH")
write.xlsx(corr_mfi_pvals,"mfi-pvals-corrected.xlsx")


####CD193
mfi_cd193 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/CD193/CD193_FLOCK.xlsx",sheet=2)
mfi_pvals <- c()
for(i in 1:18) {
  if(i == 1) { 
    mfi_pvals <- assign(paste0("mm_",i),mixed(mfi_cd193[,c(1:6,i+6)],paste0("Pop",i)))
  } else {
    a <- assign(paste0("mm_",i),mixed(mfi_cd193[,c(1:6,i+6)],paste0("Pop",i)))
    mfi_pvals <- cbind(mfi_pvals,a[,4])
  }
}
colnames(mfi_pvals)[4:21] <- paste0("Pop",seq(1,18,1))
write.xlsx(mfi_pvals,"/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/p-vals_rs_cd193_mfi.xlsx",rownames=T)

corr_193_counts <- apply(mfi_pvals[,4:21],1,p.adjust,"BH")

write.xlsx(mfi_pvals,"FLOCK_CD193mfi_pvals.xlsx")
write.xlsx(corr_193_counts,"flock_cd193mfi_corrected.xlsx",rownames=T)

data <- mfi_cd203[,c(1:6,i+6)]

##CD203
mfi_cd203 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/CD203/CD203_FLOCK.xlsx",sheet = 2)
mfi_pvals_203 <- c()
for(i in 1:17) {
  if(i == 1) { 
    mfi_pvals_203 <- assign(paste0("mm_",i),mixed(mfi_cd203[,c(1:6,i+6)],paste0("Pop",i)))
  } else {
    a <- assign(paste0("mm_",i),mixed(mfi_cd203[,c(1:6,i+6)],paste0("Pop",i)))
    mfi_pvals_203 <- cbind(mfi_pvals_203,a[,4])
  }
}
colnames(mfi_pvals_203)[4:20] <- paste0("Pop",seq(1,17,1))
write.xlsx(mfi_pvals_203,"/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/Combo/p-vals_rs_cd203_mfi.xlsx",rownames=T)


corr_203_mfi <- apply(mfi_pvals_203[,4:20],1,p.adjust,"BH")
write.xlsx(mfi_pvals_203,"/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/FLOCK_CD203mfi_pvals.xlsx")
write.xlsx(corr_203_mfi,"/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK/flock_cd203mfi_corrected.xlsx",rownames=T)


mixed <- function(data,res){
  lst <- list(list())
  data <- as.data.frame(unclass(data))
  data$Training <- factor(data$Training,levels=c("Pre","Post"))
  data$Training <- factor(data$Training,levels=c("Pre","Post"))
  data$Interaction <- as.numeric(as.factor(interaction(data$Time_point,data$Training)))
  data$Time_point.cont <- as.numeric(as.factor(data$Time_point))
  data$Training.cont <- as.numeric(as.factor(data$Training))
  data$Condition.cont <- as.factor(as.factor(data$Condition))
  
  formula_3w <- data[,res] ~ Time_point+Training+Condition+Time_point*Training+Time_point*Condition+
    Training*Condition+Training*Condition*Time_point+(1+Time_point.cont+Training.cont|Subject)
  
  model_3w <- lmer(formula_3w,data=data)
  
  #Post-hoc tests
  
  ls6 <- lsmeans(model_3w,~Condition|Training*Time_point)
  test6 <- emmeans::contrast(ls6,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  
  ls7 <- lsmeans(model_3w,~ Condition*Time_point|Training)
  test7 <- emmeans::contrast(ls7,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  
  ls8 <- lsmeans(model_3w,~ Time_point|Training*Condition)
  test8 <- emmeans::contrast(ls8,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  
  tab6 <- cbind(as.data.frame(summary(test6)[1]),as.data.frame(summary(test6)[2]),as.data.frame(summary(test6)[3]),summary(test6)[,8])
  tab7 <- cbind(as.data.frame(summary(test7)[1]),as.data.frame(summary(test7)[2]),as.data.frame(summary(test7)[3]),summary(test7)[,8])
  tab8 <- cbind(as.data.frame(summary(test8)[1]),as.data.frame(summary(test8)[2]),as.data.frame(summary(test8)[3]),summary(test8)[,8])
  colnames(tab6)<-colnames(tab7)<-colnames(tab8)<-c("Contrast","Level2","Level3","p-value")
  p_summary <- rbind(tab6,tab7,tab8)
  return(p_summary)
}


##Function for generating negative binomial models with two way and three way interaction terms
mixed_models <- function(data,res){
  lst <- list(list())
  data <- as.data.frame(unclass(data))
  data$Training <- factor(data$Training,levels=c("Pre","Post"))
  
  formula_2w <- data[,res] ~ Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
    Training*Gender+Training*Condition+Gender*Condition+(1+Training|Subject)
  
  formula_3w <- data[,res] ~ Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
    Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Condition*Gender+
    Training*Gender*Time_point+Condition*Gender*Time_point+(1+Training|Subject)
  
  model_2w <- lmer(formula_2w,data=data)
  model_3w <- lmer(formula_3w,data=data)
  
  #Post-hoc tests
  ls1 <- emmeans(model_2w,~Time_point|Training)
  test1 <- emmeans::contrast(ls1,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  ls2 <- emmeans(model_2w,~Training|Time_point)
  test2 <- emmeans::contrast(ls2,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  ls3 <- emmeans(model_2w,~Condition|Time_point)
  test3 <- emmeans::contrast(ls3,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  
  ls4 <- lsmeans(model_3w,~Time_point|Training*Condition)
  test4 <- emmeans::contrast(ls4,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  ls5 <- lsmeans(model_3w,~Gender|Training*Time_point)
  test5 <- emmeans::contrast(ls5,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  ls6 <- lsmeans(model_3w,~Condition|Training*Time_point)
  test6 <- emmeans::contrast(ls6,interaction = "trt.vs.ctrl1",adjust="Dunnett")
  ls7 <- lsmeans(model_3w,~ Training|Time_point*Condition)
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
  
  lst[[1]] <- p_summary
  lst[[2]] <- model_2w
  lst[[3]] <- model_3w
  lst[[4]] <- list(ls1,ls2,ls3,ls4,ls4,ls6,ls7,ls8,ls9)
  
  return(lst)
}

