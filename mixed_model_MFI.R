############### Mixed model for MFI ##############
library(openxlsx)
library(lme4)
library("nlme")
library(phia)
library("ggplot2")
library("reshape")
library(dplyr)
library(plyr)

############### Read the MFI tables ################
mfi_cd193 <- read.xlsx("/Volumes/Samsung_T5/UCI/Asthma/Oct11th/MFI_all_panels.xlsx",sheet=1)

#convert the characters to factors
mfi_cd193 <- as.data.frame(unclass(mfi_cd193))
model0 <- lmer(GRpos_CD193pos~Gender+Time_point+Condition+Training+(1+Training|Subject),data = mfi_cd193,REML = T)
model1 <- lmer(GRpos_CD193pos~Gender*Time_point*Condition*Training+(1+Training|Subject),data = mfi_cd193,REML = T)
summary(model1)
model2 <- lmer(GRpos_CD193pos~Gender*Time_point*Condition*Training+(1|Subject),data = mfi_cd193,REML = T)
summary(model2)
model3 <- lmer(GRpos_CD193pos~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
               Gender*Training+Time_point*Condition+Time_point*Training+Condition*Training+Gender*Time_point*Training+
                 Gender*Time_point*Condition+Time_point*Condition*Training+Gender*Condition*Training+
                 (1+Training|Subject),data = mfi_cd193,REML = T)
model4 <- lmer(GRpos_CD193pos~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
                 Gender*Training+Time_point*Training+Condition*Training+
                 Time_point*Condition+
                 (1+Training||Subject),data = mfi_cd193,REML = T)
model5 <- lmer(GRpos_CD193pos~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
                 Gender*Training+Time_point*Condition+Time_point*Training+Condition*Training+
                 (1+Training|Subject),data = mfi_cd193,REML = T)
summary(model5) 
anova(model3)
fit1 <- lme(GRpos_CD193pos~Gender*Time_point*Condition*Training,random=~+1|Subject,data=mfi_cd193)
anova(fit1)
summary(fit1)
fit2 <- lme(GRpos_CD193pos~Gender*Time_point*Condition*Training,random=~+1|Subject,data=mfi_cd193)
anova(fit2)
summary(fit2)

#Actual values

quartz()
ggplot(mfi_cd193)+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=GRpos_CD193pos,colour=Training))+
  geom_line(data=mfi_cd193,aes(x=Time_point,y=GRpos_CD193pos,colour=Training,group=Training))+
  ylab("MFI of GR+CD193+")+xlab("Time point")+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  ggtitle("CD193 MFI values")

#Predicted values
quartz()
ggplot(mfi_cd193)+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=GRpos_CD193pos,colour=Training))+
  geom_line(data=cbind(mfi_cd193,pred=fitted(model3)),aes(x=Time_point,y=pred,colour=Training,group=Training))+
  ylab("MFI")+xlab("Time point")+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

summ_cd193 <- summarySE(mfi_cd193, measurevar="GRpos_CD193pos", groupvars=c("Training","Time_point"))
pd <- position_dodge(0.4)
quartz()
ggplot(summ_cd193, aes(x=Time_point,y=GRpos_CD193pos,linetype=Training,group=Training))+ 
  geom_errorbar(aes(ymin=GRpos_CD193pos-se, ymax=GRpos_CD193pos+se),colour="black", 
                width=0.1,size=0.7,stat="identity", position = pd) +
  geom_line(position=pd,size=1) +
  geom_point(position=pd,size=3,shape=21)+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  ylim(c(2000,3000))+theme_bw()+
  theme(legend.justification=c(1,0),
        legend.position=c(1,0)) +
  ylab("MFI")+
  xlab("Time Point")+
  ggtitle("Mean-SE plot for GR+ CD193+ MFI")

emmip(model3,Training~Time_point|Condition*Gender)

#Post-hoc tests
cd193mfi_ls1 <- lsmeans(model4,~Time_point|Training)
cd193mfi_test1 <- contrast(cd193mfi_ls1,interaction = "trt.vs.ctrl1")

cd193mfi_ls2 <- lsmeans(model4,~Training|Time_point)
cd193mfi_test2 <- contrast(cd193mfi_ls2,interaction = "trt.vs.ctrl1")

cd193mfi_ls3 <- lsmeans(model3,~Time_point|Training*Condition)
cd193mfi_test3 <- contrast(cd193mfi_ls3,interaction = "trt.vs.ctrl1")

cd193mfi_ls4 <- lsmeans(model3,~Gender|Training*Time_point)
cd193mfi_test4 <- contrast(cd193mfi_ls4,interaction = "trt.vs.ctrl1")

cd193mfi_ls5 <- lsmeans(model3,~Condition|Training*Time_point)
cd193mfi_test5 <- contrast(cd193mfi_ls5,interaction = "trt.vs.ctrl1")

cd193mfi_ls6 <- lsmeans(model3,~ Training|Time_point*Condition)
cd193mfi_test6 <- contrast(cd193mfi_ls6,interaction = "trt.vs.ctrl1")

quartz()
emmip(model3,Training~Condition|Time_point,CIs=TRUE)

quartz()
emmip(model3,Training~Gender|Time_point,CIs=TRUE)



###################CD203####################
mfi_cd203 <- read.xlsx("/Volumes/Samsung_T5/UCI/Asthma/Oct11th/MFI_all_panels.xlsx",sheet=2)
mfi_cd203 <- as.data.frame(unclass(mfi_cd203))
model1_cd203 <- lmer(Grpos_CD203pos~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
                 Gender*Training+Time_point*Condition+Time_point*Training+Condition*Training+Gender*Time_point*Training+
                   Gender*Time_point*Condition+Time_point*Condition*Training+Gender*Condition*Training+
                 (1+Training|Subject),data = mfi_cd203,REML = T)
model2_cd203 <- lmer(Grpos_CD203pos~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
                       Gender*Training+Time_point*Condition+Time_point*Training+Condition*Training+
                       (1+Training|Subject),data = mfi_cd203,REML = T)
quartz()
ggplot(mfi_cd203)+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=Grpos_CD203pos,colour=Training))+
  geom_line(data=mfi_cd203,aes(x=Time_point,y=Grpos_CD203pos,colour=Training,group=Training))+
  ylab("MFI of GR+CD203+")+xlab("Time point")+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  ggtitle("CD203 MFI values")
quartz()
emmip(model1_cd203,Training~Condition|Time_point,CIs=TRUE)

quartz()
emmip(model1_cd203,Training~Gender|Time_point,CIs=TRUE)

summ_cd203 <- summarySE(mfi_cd203, measurevar="Grpos_CD203pos", groupvars=c("Training","Time_point"))
pd <- position_dodge(0.4)
quartz()
ggplot(summ_cd203, aes(x=Time_point,y=Grpos_CD203pos,linetype=Training,group=Training))+ 
  geom_errorbar(aes(ymin=Grpos_CD203pos-se, ymax=Grpos_CD203pos+se),colour="black", 
                width=0.1,size=0.7,stat="identity", position = pd) +
  geom_line(position=pd,size=1) +
  geom_point(position=pd,size=3,shape=21)+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  ylim(c(2000,3000))+theme_bw()+
  theme(legend.justification=c(1,0),
        legend.position=c(1,0)) +
  ylab("MFI")+
  xlab("Time Point")+
  ggtitle("Mean-SE plot for GR+ CD203+ MFI")

#Post-hoc tests
cd203mfi_ls1 <- lsmeans(model2_cd203,~Time_point|Training)
cd203mfi_test1 <- contrast(cd203mfi_ls1,interaction = "trt.vs.ctrl1")

cd203mfi_ls2 <- lsmeans(model2_cd203,~Training|Time_point)
cd203mfi_test2 <- contrast(cd203mfi_ls2,interaction = "trt.vs.ctrl1")

cd203mfi_ls3 <- lsmeans(model1_cd203,~Time_point|Training*Condition)
cd203mfi_test3 <- contrast(cd203mfi_ls3,interaction = "trt.vs.ctrl1")

cd203mfi_ls4 <- lsmeans(model1_cd203,~Gender|Training*Time_point)
cd203mfi_test4 <- contrast(cd203mfi_ls4,interaction = "trt.vs.ctrl1")

cd203mfi_ls5 <- lsmeans(model1_cd203,~Condition|Training*Time_point)
cd203mfi_test5 <- contrast(cd203mfi_ls5,interaction = "trt.vs.ctrl1")

cd203mfi_ls6 <- lsmeans(model1_cd203,~ Training|Time_point*Condition)
cd203mfi_test6 <- contrast(cd203mfi_ls6,interaction = "trt.vs.ctrl1")

quartz()
emmip(model1_cd203,Training~Gender|Time_point,CIs=TRUE,main="Interaction plot of CD203 Male vs Female")
####################Combo CD3 ######################
mfi_combo <- read.xlsx("/Volumes/Samsung_T5/UCI/Asthma/Oct11th/MFI_all_panels.xlsx",sheet=3)
mfi_combo <- as.data.frame(unclass(mfi_combo))
mfi_cd3 <- mfi_combo[,c(1:7,11)]
model1_cd3 <- lmer(GRpos_CD3pos~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
                       Gender*Training+Time_point*Condition+Time_point*Training+Condition*Training+
                       (1+Training|Subject),data = mfi_cd3,REML = T)

model2_cd3 <- lmer(GRpos_CD3pos~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
                 Gender*Training+Time_point*Condition+Time_point*Training+Condition*Training+Gender*Time_point*Training+
                 Gender*Time_point*Condition+Time_point*Condition*Training+Gender*Condition*Training+
                 (1+Training|Subject),data = mfi_cd3,REML = T)

model3_cd3 <- lmer(GRpos_CD3pos~Gender*Time_point*Condition*Training+(1+Training|Subject),data = mfi_cd3,REML = T)

#Post-hoc tests
cd3mfi_ls1 <- lsmeans(model1_cd3,~Time_point|Training)
cd3mfi_test1 <- contrast(cd3mfi_ls1,interaction = "trt.vs.ctrl1")

cd3mfi_ls2 <- lsmeans(model1_cd3,~Training|Time_point)
cd3mfi_test2 <- contrast(cd3mfi_ls2,interaction = "trt.vs.ctrl1")

cd3mfi_ls3 <- lsmeans(model2_cd3,~Time_point|Training*Condition)
cd3mfi_test3 <- contrast(cd3mfi_ls3,interaction = "trt.vs.ctrl1")

cd3mfi_ls4 <- lsmeans(model2_cd3,~Condition|Training*Time_point)
cd3mfi_test4 <- contrast(cd3mfi_ls4,interaction = "trt.vs.ctrl1")

cd3mfi_ls7 <- lsmeans(model2_cd3,~Condition*Training|Time_point)
cd3mfi_test7 <- contrast(cd3mfi_ls7,interaction = "trt.vs.ctrl1")

cd3mfi_ls8 <- lsmeans(model2_cd3,~Gender*Training|Time_point)
cd3mfi_test8 <- contrast(cd3mfi_ls8,interaction = "trt.vs.ctrl1")

cd3mfi_ls5 <- lsmeans(model2_cd3,~Gender|Training*Time_point)
cd3mfi_test5 <- contrast(cd3mfi_ls5,interaction = "trt.vs.ctrl1")

cd3mfi_ls6 <- lsmeans(model2_cd3,~ Training|Time_point*Condition)
cd3mfi_test6 <- contrast(cd3mfi_ls6,interaction = "trt.vs.ctrl1")

summary(model3_cd3)

quartz()
emmip(model2_cd3,Training~Condition|Time_point,CIs=TRUE)

quartz()
emmip(model2_cd3,Training~Gender|Time_point,CIs=TRUE)

quartz()
emmip(model2_cd3,Training~Condition|Time_point*Training,CIs=TRUE)

summ_cd3 <- summarySE(mfi_cd3, measurevar="GRpos_CD3pos", groupvars=c("Training","Time_point"))
pd <- position_dodge(0.4)
quartz()
ggplot(summ_cd3, aes(x=Time_point,y=GRpos_CD3pos,linetype=Training,group=Training))+ 
  geom_errorbar(aes(ymin=GRpos_CD3pos-se, ymax=GRpos_CD3pos+se),colour="black", 
                width=0.1,size=0.7,stat="identity", position = pd) +
  geom_line(position=pd,size=1) +
  geom_point(position=pd,size=3,shape=21)+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_bw()+
  ggtitle("Mean-SE plot for GR+ CD3+ MFI")+
  ylim(2000,3000)+
  theme(legend.justification=c(1,0),
        legend.position=c(1,0)) +
  ylab("MFI")+
  xlab("Time Point")


########################Combo CD14################################
mfi_cd14 <- mfi_combo[,c(1:7,10)]
model1_cd14 <- lmer(GRpos_CD14pos~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
                     Gender*Training+Time_point*Condition+Time_point*Training+Condition*Training+
                     (1+Training|Subject),data = mfi_cd14,REML = T)

model2_cd14 <- lmer(GRpos_CD14pos~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
                     Gender*Training+Time_point*Condition+Time_point*Training+Condition*Training+Gender*Time_point*Training+
                     Gender*Time_point*Condition+Time_point*Condition*Training+Gender*Condition*Training+
                     (1+Training|Subject),data = mfi_cd14,REML = T)

#Post-hoc tests
cd14mfi_ls1 <- lsmeans(model1_cd14,~Time_point|Training)
cd14mfi_test1 <- contrast(cd14mfi_ls1,interaction = "trt.vs.ctrl1")

cd14mfi_ls2 <- lsmeans(model1_cd14,~Training|Time_point)
cd14mfi_test2 <- contrast(cd14mfi_ls2,interaction = "trt.vs.ctrl1")

cd14mfi_ls3 <- lsmeans(model2_cd14,~Time_point|Training*Condition)
cd14mfi_test3 <- contrast(cd203mfi_ls3,interaction = "trt.vs.ctrl1")

cd14mfi_ls4 <- lsmeans(model2_cd14,~Condition|Training*Time_point)
cd14mfi_test4 <- contrast(cd14mfi_ls4,interaction = "trt.vs.ctrl1")

cd14mfi_ls5 <- lsmeans(model2_cd14,~Gender|Training*Time_point)
cd14mfi_test5 <- contrast(cd14mfi_ls5,interaction = "trt.vs.ctrl1")

cd14mfi_ls7 <- lsmeans(model2_cd14,~Gender*Time_point|Training)
cd14mfi_test7 <- contrast(cd14mfi_ls7,interaction = "trt.vs.ctrl1")

cd14mfi_ls6 <- lsmeans(model2_cd14,~ Training|Time_point*Condition)
cd14mfi_test6 <- contrast(cd14mfi_ls6,interaction = "trt.vs.ctrl1")

quartz()
emmip(model2_cd14,Training~Condition|Time_point,CIs=TRUE)

quartz()
emmip(model2_cd14,Training~Gender|Time_point,CIs=TRUE)

quartz()
emmip(model2_cd14,Training~Gender|Time_point*Training,CIs=TRUE)

summ_cd14 <- summarySE(mfi_cd14, measurevar="GRpos_CD14pos", groupvars=c("Training","Time_point"))
pd <- position_dodge(0.4)
quartz()
ggplot(summ_cd14, aes(x=Time_point,y=GRpos_CD14pos,linetype=Training,group=Training))+ 
  geom_errorbar(aes(ymin=GRpos_CD14pos-se, ymax=GRpos_CD14pos+se),colour="black", 
                width=0.1,size=0.7,stat="identity", position = pd) +
  geom_line(position=pd,size=1) +
  geom_point(position=pd,size=3,shape=21)+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_bw()+
  ggtitle("Mean-SE plot for GR+ CD14+ MFI")+
  ylim(2000,3000)+
  theme(legend.justification=c(1,0),
        legend.position=c(1,0)) +
  ylab("MFI")+
  xlab("Time Point")

#############################Combo CD16pos##########################
mfi_cd16pos <- mfi_combo[,c(1:7,9)]
model1_cd16pos <- lmer(GRpos_CD16pos~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
                         Gender*Training+Time_point*Condition+Time_point*Training+Condition*Training+
                         (1+Training|Subject),data = mfi_cd16pos,REML = T)

model2_cd16pos <- lmer(GRpos_CD16pos~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
                      Gender*Training+Time_point*Condition+Time_point*Training+Condition*Training+Gender*Time_point*Training+
                      Gender*Time_point*Condition+Time_point*Condition*Training+Gender*Condition*Training+
                      (1+Training|Subject),data = mfi_cd16pos,REML = T)

model3_cd16 <- lmer(GRpos_CD16pos~Gender*Time_point*Condition*Training+(1+Training|Subject),data = mfi_cd16pos,REML = T)

#Post-hoc tests
cd16mfi_ls1 <- lsmeans(model1_cd16pos,~Time_point|Training)
cd16mfi_test1 <- contrast(cd16mfi_ls1,interaction = "trt.vs.ctrl1")

cd16mfi_ls2 <- lsmeans(model1_cd16pos,~Training|Time_point)
cd16mfi_test2 <- contrast(cd16mfi_ls2,interaction = "trt.vs.ctrl1")

cd16mfi_ls3 <- lsmeans(model2_cd16pos,~Time_point|Training*Condition)
cd16mfi_test3 <- contrast(cd16mfi_ls3,interaction = "trt.vs.ctrl1")

cd16mfi_ls4 <- lsmeans(model2_cd16pos,~Condition|Training*Time_point)
cd16mfi_test4 <- contrast(cd16mfi_ls4,interaction = "trt.vs.ctrl1")

cd16mfi_ls5 <- lsmeans(model2_cd16pos,~Gender|Training*Time_point)
cd16mfi_test5 <- contrast(cd16mfi_ls5,interaction = "trt.vs.ctrl1")

cd16mfi_ls6 <- lsmeans(model2_cd16pos,~ Training|Time_point*Condition)
cd16mfi_test6 <- contrast(cd16mfi_ls6,interaction = "trt.vs.ctrl1")

quartz()
emmip(model2_cd16pos,Training~Condition|Time_point,CIs=TRUE)

quartz()
emmip(model2_cd16pos,Training~Gender|Time_point,CIs=TRUE)

quartz()
emmip(model2_cd16pos,Training~Time_point|Condition,CIs=TRUE)

summ_cd16pos <- summarySE(mfi_cd16pos, measurevar="GRpos_CD16pos", groupvars=c("Training","Time_point"))
pd <- position_dodge(0.4)
quartz()
ggplot(summ_cd16pos, aes(x=Time_point,y=GRpos_CD16pos,linetype=Training,group=Training))+ 
  geom_errorbar(aes(ymin=GRpos_CD16pos-se, ymax=GRpos_CD16pos+se),colour="black", 
                width=0.1,size=0.7,stat="identity", position = pd) +
  geom_line(position=pd,size=1) +
  geom_point(position=pd,size=3,shape=21)+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_bw()+
  ggtitle("Mean-SE plot for GR+ CD16+ MFI")+
  ylim(2000,3000)+
  theme(legend.justification=c(1,0),
        legend.position=c(1,0)) +
  ylab("MFI")+
  xlab("Time Point")

################################Combo CD16hi ##############################
mfi_cd16hi <- mfi_combo[,c(1:8)]
model1_cd16hi <- lmer(GRpos_CD16Hi~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
                         Gender*Training+Time_point*Condition+Time_point*Training+Condition*Training+
                         (1+Training|Subject),data = mfi_cd16hi,REML = T)

model2_cd16hi <- lmer(GRpos_CD16Hi~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
                         Gender*Training+Time_point*Condition+Time_point*Training+Condition*Training+Gender*Time_point*Training+
                         Gender*Time_point*Condition+Time_point*Condition*Training+Gender*Condition*Training+
                         (1+Training|Subject),data = mfi_cd16hi,REML = T)

model3_cd16hi <- lmer(GRpos_CD16Hi~Gender*Time_point*Condition*Training+(1+Training|Subject),data = mfi_cd16hi,REML = T)
model_cd16hi <- lmer(GRpos_CD16Hi~Gender*Time_point*Condition*Training+(1|Subject),data = mfi_cd16hi,REML = T)
#Post-hoc tests
cd16hi_mfi_ls1 <- lsmeans(model1_cd16hi,~Time_point|Training)
cd16hi_mfi_test1 <- contrast(cd16hi_mfi_ls1,interaction = "trt.vs.ctrl1")

cd16hi_mfi_ls2 <- lsmeans(model1_cd16hi,~Training|Time_point)
cd16hi_mfi_test2 <- contrast(cd16hi_mfi_ls2,interaction = "trt.vs.ctrl1")

cd16hi_mfi_ls3 <- lsmeans(model2_cd16hi,~ Time_point|Training*Condition)
cd16hi_mfi_test3 <- contrast(cd16hi_mfi_ls3,interaction = "trt.vs.ctrl1")

cd16hi_mfi_ls4 <- lsmeans(model2_cd16hi,~Condition|Training*Time_point)
cd16hi_mfi_test4 <- contrast(cd16hi_mfi_ls4,interaction = "trt.vs.ctrl1")

cd16hi_mfi_ls5 <- lsmeans(model2_cd16hi,~Gender|Training*Time_point)
cd16hi_mfi_test5 <- contrast(cd16hi_mfi_ls5,interaction = "trt.vs.ctrl1")

cd16hi_mfi_ls6 <- lsmeans(model2_cd16hi,~ Training|Time_point*Condition)
cd16hi_mfi_test6 <- contrast(cd16hi_mfi_ls6,interaction = "trt.vs.ctrl1")

quartz()
emmip(model2_cd16hi,Training~Condition|Time_point,CIs=TRUE)

quartz()
emmip(model2_cd16hi,Training~Gender|Time_point,CIs=TRUE)

quartz()
emmip(model2_cd16hi,Training~Time_point|Condition,CIs=TRUE)

summ_cd16hi <- summarySE(mfi_cd16hi, measurevar="GRpos_CD16Hi", groupvars=c("Training","Time_point"))
pd <- position_dodge(0.4)
quartz()
ggplot(summ_cd16hi, aes(x=Time_point,y=GRpos_CD16Hi,linetype=Training,group=Training))+ 
  geom_errorbar(aes(ymin=GRpos_CD16Hi-se, ymax=GRpos_CD16Hi+se),colour="black", 
                width=0.1,size=0.7,stat="identity", position = pd) +
  geom_line(position=pd,size=1) +
  geom_point(position=pd,size=3,shape=21)+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_bw()+
  ggtitle("Mean-SE plot for GR+ CD16Hi MFI")+
  ylim(2000,3000)+
  theme(legend.justification=c(1,0),
        legend.position=c(1,0)) +
  ylab("MFI")+
  xlab("Time Point")
##############Monocytes#####################
mfi_monocytes <- read.xlsx("/Volumes/Samsung_T5/UCI/Asthma/Oct11th/MFI_all_panels.xlsx",sheet=4)
classical_monocytes <- mfi_monocytes[,c(1:8)]
intermediate_monocytes <- mfi_monocytes[,c(1:7,9)]
non_c_monocytes <- mfi_monocytes[,c(1:7,10)]

classical_model1_monocytes <- lmer(classical_monocytes[,8]~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
           Gender*Training+Time_point*Condition+Time_point*Training+Condition*Training+
           (1+Training|Subject),data = classical_monocytes,REML = T)

classical_model2_monocytes <- lmer(classical_monocytes[,8]~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
                        Gender*Training+Time_point*Condition+Time_point*Training+Condition*Training+Gender*Time_point*Training+
                        Gender*Time_point*Condition+Time_point*Condition*Training+Gender*Condition*Training+
                        (1+Training|Subject),data = classical_monocytes,REML = T)

classical_model3_monocytes <- lmer(classical_monocytes[,8]~Gender*Time_point*Condition*Training+(1+Training|Subject),data = classical_monocytes,REML = T)

classical_mfi_ls1 <- lsmeans(classical_model1_monocytes,~Time_point|Training)
classical_mfi_test1 <- contrast(classical_mfi_ls1,interaction = "trt.vs.ctrl1")

classical_mfi_ls2 <- lsmeans(classical_model1_monocytes,~Training|Time_point)
classical_mfi_test2 <- contrast(classical_mfi_ls2,interaction = "trt.vs.ctrl1")

classical_mfi_ls3 <- lsmeans(classical_model2_monocytes,~ Time_point|Training*Condition)
classical_mfi_test3 <- contrast(classical_mfi_ls3,interaction = "trt.vs.ctrl1")

classical_mfi_ls4 <- lsmeans(classical_model2_monocytes,~Condition|Training*Time_point)
classical_mfi_test4 <- contrast(classical_mfi_ls4,interaction = "trt.vs.ctrl1")

classical_mfi_ls5 <- lsmeans(classical_model2_monocytes,~Gender|Training*Time_point)
classical_mfi_test5 <- contrast(classical_mfi_ls5,interaction = "trt.vs.ctrl1")

classical_mfi_ls6 <- lsmeans(classical_model2_monocytes,~ Training|Time_point*Condition)
classical_mfi_test6 <- contrast(classical_mfi_ls6,interaction = "trt.vs.ctrl1")

quartz()
emmip(classical_model2_monocytes,Training~Condition|Time_point,CIs=TRUE)

quartz()
emmip(classical_model2_monocytes,Training~Gender|Time_point,CIs=TRUE)

summ_class <- summarySE(classical_monocytes, measurevar="Classical_Backgated", groupvars=c("Training","Time_point"))
pd <- position_dodge(0.4)
quartz()
ggplot(summ_class, aes(x=Time_point,y=Classical_Backgated,linetype=Training,group=Training))+ 
  geom_errorbar(aes(ymin=Classical_Backgated-se, ymax=Classical_Backgated+se),colour="black", 
                width=0.1,size=0.7,stat="identity", position = pd) +
  geom_line(position=pd,size=1) +
  geom_point(position=pd,size=3,shape=21)+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_bw()+
  ggtitle("Mean-SE plot for GR expression on Classical Monocytes")+
  ylim(2000,3000)+
  theme(legend.justification=c(1,0),
        legend.position=c(1,0)) +
  ylab("MFI")+
  xlab("Time Point")

int_model1_monocytes <- lmer(intermediate_monocytes[,8]~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
                                     Gender*Training+Time_point*Condition+Time_point*Training+Condition*Training+
                                     (1+Training|Subject),data = intermediate_monocytes,REML = T)

int_model2_monocytes <- lmer(intermediate_monocytes[,8]~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
                                     Gender*Training+Time_point*Condition+Time_point*Training+Condition*Training+Gender*Time_point*Training+
                                     Gender*Time_point*Condition+Time_point*Condition*Training+Gender*Condition*Training+
                                     (1+Training|Subject),data = intermediate_monocytes,REML = T)

int_model3_monocytes <- lmer(intermediate_monocytes[,8]~Gender*Time_point*Condition*Training+(1+Training|Subject),data = intermediate_monocytes,REML = T)

int_mfi_ls1 <- lsmeans(int_model1_monocytes,~Time_point|Training)
int_mfi_test1 <- contrast(int_mfi_ls1,interaction = "trt.vs.ctrl1")

int_mfi_ls2 <- lsmeans(int_model1_monocytes,~Training|Time_point)
int_mfi_test2 <- contrast(int_mfi_ls2,interaction = "trt.vs.ctrl1")

int_mfi_ls3 <- lsmeans(int_model2_monocytes,~ Time_point|Training*Condition)
int_mfi_test3 <- contrast(int_mfi_ls3,interaction = "trt.vs.ctrl1")

int_mfi_ls4 <- lsmeans(int_model2_monocytes,~Condition|Training*Time_point)
int_mfi_test4 <- contrast(int_mfi_ls4,interaction = "trt.vs.ctrl1")

int_mfi_ls5 <- lsmeans(int_model2_monocytes,~Gender|Training*Time_point)
int_mfi_test5 <- contrast(int_mfi_ls5,interaction = "trt.vs.ctrl1")

int_mfi_ls6 <- lsmeans(int_model2_monocytes,~ Training|Time_point*Condition)
int_mfi_test6 <- contrast(int_mfi_ls6,interaction = "trt.vs.ctrl1")

quartz()
emmip(int_model2_monocytes,Training~Condition|Time_point,CIs=TRUE)

quartz()
emmip(int_model2_monocytes,Training~Gender|Time_point,CIs=TRUE)

summ_int <- summarySE(intermediate_monocytes, measurevar="Intermediate_Monocytes", groupvars=c("Training","Time_point"))
pd <- position_dodge(0.4)
quartz()
ggplot(summ_int, aes(x=Time_point,y=Intermediate_Monocytes,linetype=Training,group=Training))+ 
  geom_errorbar(aes(ymin=Intermediate_Monocytes-se, ymax=Intermediate_Monocytes+se),colour="black", 
                width=0.1,size=0.7,stat="identity", position = pd) +
  geom_line(position=pd,size=1) +
  geom_point(position=pd,size=3,shape=21)+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_bw()+
  ggtitle("Mean-SE plot for GR expression on Intermediate Monocytes")+
  ylim(2000,3000)+
  theme(legend.justification=c(1,0),
        legend.position=c(1,0)) +
  ylab("MFI")+
  xlab("Time Point")

non_c_model1_monocytes <- lmer(non_c_monocytes[,8]~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
                                     Gender*Training+Time_point*Condition+Time_point*Training+Condition*Training+
                                     (1+Training|Subject),data = non_c_monocytes,REML = T)

non_c_model2_monocytes <- lmer(non_c_monocytes[,8]~Gender+Time_point+Condition+Training+Gender*Time_point+Gender*Condition+
                                     Gender*Training+Time_point*Condition+Time_point*Training+Condition*Training+Gender*Time_point*Training+
                                     Gender*Time_point*Condition+Time_point*Condition*Training+Gender*Condition*Training+
                                     (1+Training|Subject),data = non_c_monocytes,REML = T)

non_c_model3_monocytes <- lmer(non_c_monocytes[,8]~Gender*Time_point*Condition*Training+(1+Training|Subject),data =non_c_monocytes,REML = T)

nonc_mfi_ls1 <- lsmeans(non_c_model1_monocytes,~Time_point|Training)
nonc_mfi_test1 <- contrast(nonc_mfi_ls1,interaction = "trt.vs.ctrl1")

nonc_mfi_ls2 <- lsmeans(non_c_model1_monocytes,~Training|Time_point)
nonc_mfi_test2 <- contrast(nonc_mfi_ls2,interaction = "trt.vs.ctrl1")

nonc_mfi_ls3 <- lsmeans(non_c_model2_monocytes,~ Time_point|Training*Condition)
nonc_mfi_test3 <- contrast(nonc_mfi_ls3,interaction = "trt.vs.ctrl1")

nonc_mfi_ls4 <- lsmeans(non_c_model2_monocytes,~Condition|Training*Time_point)
nonc_mfi_test4 <- contrast(nonc_mfi_ls4,interaction = "trt.vs.ctrl1")

nonc_mfi_ls5 <- lsmeans(non_c_model2_monocytes,~Gender|Training*Time_point)
nonc_mfi_test5 <- contrast(nonc_mfi_ls5,interaction = "trt.vs.ctrl1")

nonc_mfi_ls6 <- lsmeans(non_c_model2_monocytes,~ Training|Time_point*Condition)
nonc_mfi_test6 <- contrast(nonc_mfi_ls6,interaction = "trt.vs.ctrl1")

quartz()
emmip(non_c_model2_monocytes,Training~Condition|Time_point,CIs=TRUE)

quartz()
emmip(non_c_model2_monocytes,Training~Gender|Time_point,CIs=TRUE)

colnames(non_c_monocytes)[8] <- c("Non_classical_Monocytes")
summ_nonc <- summarySE(non_c_monocytes, measurevar="Non_classical_Monocytes", groupvars=c("Training","Time_point"))
pd <- position_dodge(0.4)
quartz()
ggplot(summ_nonc, aes(x=Time_point,y=Non_classical_Monocytes,linetype=Training,group=Training))+ 
  geom_errorbar(aes(ymin=`Non_classical_Monocytes`-se, ymax=`Non_classical_Monocytes`+se),colour="black", 
                width=0.1,size=0.7,stat="identity", position = pd) +
  geom_line(position=pd,size=1) +
  geom_point(position=pd,size=3,shape=21)+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_bw()+
  ggtitle("Mean-SE plot for GR expression on Non-classical Monocytes")+
  ylim(2000,3000)+
  theme(legend.justification=c(1,0),
        legend.position=c(1,0)) +
  ylab("MFI")+
  xlab("Time Point")
