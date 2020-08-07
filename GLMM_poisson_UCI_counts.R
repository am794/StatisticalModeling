##UCI normalized counts data: Poisson family
library(openxlsx)
library(lme4)
library("nlme")
library(phia)
library("ggplot2")
library("reshape")
library(dplyr)
library(plyr)
library(reshape2)
library("MASS")

################################### CD193 normalized percentages ##########################################
CD193_events <-  read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=3)
CD193_events <- as.data.frame(unclass(CD193_events))
cbc_193_203 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=5)
norm_CD193 <- (CD193_events$CD193posGRpos*cbc_193_203[,6]*100)/(CD193_events$Leukocytes)
norm_CD193_counts <- cbind(CD193_events,norm_CD193)

CD193_pos_mod1b <- glmer.nb(norm_CD193 ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                          Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Condition*Gender+
                          Training*Gender*Time_point+Condition*Gender*Time_point+(1|Subject)+
                          (1+Training|Subject),data=norm_CD193_counts)

CD193_pos_mod2b <- glmer.nb(norm_CD193 ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                          Training*Gender+Training*Condition+Gender*Condition+(1|Subject)+
                          (1+Training|Subject),data=norm_CD193_counts)

CD193_pos_mod1 <- glmmPQL(norm_CD193 ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                          Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Condition*Gender+
                          Training*Gender*Time_point+Condition*Gender*Time_point,random=~1+Training|Subject,family=quasipoisson,data = norm_CD193_counts)

CD193_pos_mod2 <- glmmPQL(norm_CD193 ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                          Training*Gender+Training*Condition+Gender*Condition,random=~1+Training|Subject,family=quasipoisson,data = norm_CD193_counts)

CD193_pos_mod3 <- glmmPQL(norm_CD193 ~Time_point*Training*Gender*Condition,random=~1|Subject/Training,family=poisson(link='log'),data = norm_CD193_counts)
#overdisp(CD193_pos_mod1)

#Post-hoc tests
cd193mfi_ls1_pos <- lsmeans(CD193_pos_mod2b,~Time_point|Training)
cd193mfi_test1_pos <- contrast(cd193mfi_ls1_pos,interaction = "trt.vs.ctrl1")

cd193mfi_ls2_pos <- lsmeans(CD193_pos_mod2b ,~Training|Time_point)
cd193mfi_test2_pos <- contrast(cd193mfi_ls2_pos,interaction = "trt.vs.ctrl1")

cd193mfi_ls3_pos <- lsmeans(CD193_pos_mod1b ,~Time_point|Training*Condition)
cd193mfi_test3_pos <- contrast(cd193mfi_ls3_pos,interaction = "trt.vs.ctrl1")

cd193mfi_ls4_pos <- lsmeans(CD193_pos_mod1b ,~Gender|Training*Time_point)
cd193mfi_test4_pos <- contrast(cd193mfi_ls4_pos,interaction = "trt.vs.ctrl1")

cd193mfi_ls5_pos <- lsmeans(CD193_pos_mod1b ,~Condition|Training*Time_point)
cd193mfi_test5_pos <- contrast(cd193mfi_ls5_pos,interaction = "trt.vs.ctrl1")

cd193mfi_ls6_pos <- lsmeans(CD193_pos_mod1b ,~ Training|Time_point*Condition)
cd193mfi_test6_pos <- contrast(cd193mfi_ls6_pos,interaction = "trt.vs.ctrl1")

quartz()
ggplot(data=cbind(norm_CD193_counts,pred=fitted(CD193_pos_mod1)))+
  facet_wrap(~Condition)+
  geom_point(aes(x=Time_point,y=pred,color=Training))+
  geom_boxplot(aes(x=Time_point,y=pred,colour=Training))+
  ylab("Normalized counts")+xlab("Time point")

quartz()
ggplot(norm_CD193_counts)+
  facet_wrap(~Condition)+
  geom_point(aes(x=Time_point,y=(norm_CD193),color=Training))+
  geom_boxplot(aes(x=Time_point,y=norm_CD193,colour=Training))+
  ylab("Normalized counts")+xlab("Time point")

quartz()
ggplot(data=cbind(norm_CD193_counts,pred=fitted(CD193_pos_mod1b,type="response")))+
  facet_wrap(~Subject)+geom_smooth(aes(x=Time_point,y=pred),method="lm",alpha=0.3)+
  geom_point(aes(x=Time_point,y=norm_CD193,colour=Training))+
  geom_line(aes(x=Time_point,y=pred,colour=Training,group=Training))+
  ylab("Normalized proportions")+xlab("Time point")+ggtitle("CD193 normalised percentages")


################################### CD203 normalized percentages ##########################################
CD203_events <-  read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=4)
cbc_193_203 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=5)
norm_CD203 <- (CD203_events$`CD203pGRp`*cbc_193_203[,5]*100)/(CD203_events$Leukocytes)
norm_CD203_counts <- cbind(CD203_events,norm_CD203)

CD203_pos_mod1b <- glmer.nb(norm_CD203 ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                          Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Condition*Gender+
                          Training*Gender*Time_point+Condition*Gender*Time_point+(1|Subject)+
                          (1+Training|Subject),data=norm_CD203_counts)

CD203_pos_mod2b <- glmer.nb(norm_CD203 ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                          Training*Gender+Training*Condition+Gender*Condition+(1|Subject)+
                          (1+Training|Subject),data=norm_CD203_counts)


CD203_pos_mod1 <- glmmPQL(norm_CD203 ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                            Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Condition*Gender+
                            Training*Gender*Time_point+Condition*Gender*Time_point,random=~1+Training|Subject,family=quasipoisson,data = norm_CD203_counts)

CD203_pos_mod2 <- glmmPQL(norm_CD203 ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                            Training*Gender+Training*Condition+Gender*Condition,random=~1+Training|Subject,family=quasipoisson,data = norm_CD203_counts)

#overdisp(CD193_pos_mod1)

#Post-hoc tests
cd203mfi_ls1_pos <- lsmeans(CD203_pos_mod2b,~Time_point|Training)
cd203mfi_test1_pos <- contrast(cd203mfi_ls1_pos,interaction = "trt.vs.ctrl1")

cd203mfi_ls2_pos <- lsmeans(CD203_pos_mod2b ,~Training|Time_point)
cd203mfi_test2_pos <- contrast(cd203mfi_ls2_pos,interaction = "trt.vs.ctrl1")

cd203mfi_ls3_pos <- lsmeans(CD203_pos_mod1b ,~Time_point|Training*Condition)
cd203mfi_test3_pos <- contrast(cd203mfi_ls3_pos,interaction = "trt.vs.ctrl1")

cd203mfi_ls4_pos <- lsmeans(CD203_pos_mod1b ,~Gender|Training*Time_point)
cd203mfi_test4_pos <- contrast(cd203mfi_ls4_pos,interaction = "trt.vs.ctrl1")

cd203mfi_ls5_pos <- lsmeans(CD203_pos_mod1b ,~Condition|Training*Time_point)
cd203mfi_test5_pos <- contrast(cd203mfi_ls5_pos,interaction = "trt.vs.ctrl1")

cd203mfi_ls6_pos <- lsmeans(CD203_pos_mod1b ,~ Training|Time_point*Condition)
cd203mfi_test6_pos <- contrast(cd203mfi_ls6_pos,interaction = "trt.vs.ctrl1")

quartz()
ggplot(data=cbind(norm_CD203_counts,pred=fitted(CD203_pos_mod2b,type="response")))+
  facet_wrap(~Subject)+geom_smooth(aes(x=Time_point,y=pred),method="lm",alpha=0.3)+
  geom_point(aes(x=Time_point,y=norm_CD203,colour=Training))+
  geom_line(aes(x=Time_point,y=pred,colour=Training,group=Training))+
  ylab("Normalized proportions")+xlab("Time point")+ggtitle("CD203 normalised percentages")

############# Read the events table ###############
CD3_event <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=1)
cbc_cd3 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=2)
CD3_events <- as.data.frame(unclass(CD3_event))

norm_CD3 <- (CD3_events$`CD3posGRpos`*cbc_cd3[,4]*100)/(CD3_events$Leukocytes)
norm_CD3_counts <- cbind(CD3_events,norm_CD3)

CD3_pos_mod1b <- glmer.nb(norm_CD3 ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                              Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Condition*Gender+
                              Training*Gender*Time_point+Condition*Gender*Time_point+(1|Subject)+
                              (1+Training|Subject),data=norm_CD3_counts)

CD3_pos_mod2b <- glmer.nb(norm_CD3 ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                              Training*Gender+Training*Condition+Gender*Condition+(1|Subject)+
                              (1+Training|Subject),data=norm_CD3_counts)

CD3_pos_mod1 <- glmmPQL(norm_CD3 ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                          Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Condition*Gender+
                          Training*Gender*Time_point+Condition*Gender*Time_point,random=~1|Subject,family=quasipoisson,data = norm_CD3_counts)

CD3_pos_mod2 <- glmmPQL(norm_CD3 ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                          Training*Gender+Training*Condition+Gender*Condition,random=~1+Training|Subject,family=quasipoisson,data = norm_CD3_counts)

cd3ct_ls1_pos <- lsmeans(CD3_pos_mod2b,~Time_point|Training)
cd3ct_test1_pos <- contrast(cd3ct_ls1_pos,interaction = "trt.vs.ctrl1")

cd3ct_ls2_pos <- lsmeans(CD3_pos_mod2b ,~Training|Time_point)
cd3ct_test2_pos <- contrast(cd3ct_ls2_pos,interaction = "trt.vs.ctrl1")

cd3ct_ls3_pos <- lsmeans(CD3_pos_mod1b ,~Time_point|Training*Condition)
cd3ct_test3_pos <- contrast(cd3ct_ls3_pos,interaction = "trt.vs.ctrl1")

cd3ct_ls4_pos <- lsmeans(CD3_pos_mod1b ,~Gender|Training*Time_point)
cd3ct_test4_pos <- contrast(cd3ct_ls4_pos,interaction = "trt.vs.ctrl1")

cd3ct_ls5_pos <- lsmeans(CD3_pos_mod1b ,~Condition|Training*Time_point)
cd3ct_test5_pos <- contrast(cd3ct_ls5_pos,interaction = "trt.vs.ctrl1")

cd3ct_ls6_pos <- lsmeans(CD3_pos_mod1b ,~ Training|Time_point*Condition)
cd3ct_test6_pos <- contrast(cd3ct_ls6_pos,interaction = "trt.vs.ctrl1")

contrast(cd3ct_ls3_pos,"consec",simple="each",combine = TRUE,adjust="",interaction="trt.vs.ctrl1")

quartz()
ggplot(data=cbind(norm_CD3_counts,pred=fitted(CD3_pos_mod2b,type="response")))+
  facet_wrap(~Subject)+geom_smooth(aes(x=Time_point,y=pred),method="lm",alpha=0.3)+
  geom_point(aes(x=Time_point,y=norm_CD3,colour=Training))+
  geom_line(aes(x=Time_point,y=pred,colour=Training,group=Training))+
  ylab("Normalized proportions")+xlab("Time point")+ggtitle("CD3 normalised percentages")

############### Combo
combo_events <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=6)
cbc_combo <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=7)
combo_events <- as.data.frame(unclass(combo_events))

###################### CD14 normalized #######################
cd14_events <- combo_events[,-c(10,11,13)]
cbc_combo <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=7)
norm_cd14 <- (cd14_events$`CD14posGRpos`*cbc_combo[,5]*100)/(cd14_events$Leukocytes)
norm_cd14_counts <- cbind(cd14_events,norm_cd14)

CD14_pos_mod1 <- glmmPQL(norm_cd14 ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                          Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Condition*Gender+
                          Training*Gender*Time_point+Condition*Gender*Time_point,random=~1|Subject,family=poisson,data = norm_cd14_counts)

CD14_pos_mod2 <- glmmPQL(norm_CD3 ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                          Training*Gender+Training*Condition+Gender*Condition,random=~1|Subject,family=quasipoisson,data = norm_cd14_counts)

CD14_pos_mod1b <- glmer.nb(norm_cd14 ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                            Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Condition*Gender+
                            Training*Gender*Time_point+Condition*Gender*Time_point+(1|Subject)+
                            (1+Training|Subject),data=norm_cd14_counts)

CD14_pos_mod2b <- glmer.nb(norm_cd14 ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                            Training*Gender+Training*Condition+Gender*Condition+(1|Subject)+
                            (1+Training|Subject),data=norm_cd14_counts)

cd14ct_ls1_pos <- lsmeans(CD14_pos_mod2b,~Time_point|Training)
cd14ct_test1_pos <- contrast(cd14ct_ls1_pos,interaction = "trt.vs.ctrl1")

cd14ct_ls2_pos <- lsmeans(CD14_pos_mod2b ,~Training|Time_point)
cd14ct_test2_pos <- contrast(cd14ct_ls2_pos,interaction = "trt.vs.ctrl1")

cd14ct_ls3_pos <- lsmeans(CD14_pos_mod1b ,~Time_point|Training*Condition)
cd14ct_test3_pos <- contrast(cd14ct_ls3_pos,interaction = "trt.vs.ctrl1")

cd14ct_ls4_pos <- lsmeans(CD14_pos_mod1b ,~Gender|Training*Time_point)
cd14ct_test4_pos <- contrast(cd14ct_ls4_pos,interaction = "trt.vs.ctrl1")

cd14ct_ls5_pos <- lsmeans(CD14_pos_mod1b ,~Condition|Training*Time_point)
cd14ct_test5_pos <- contrast(cd14ct_ls5_pos,interaction = "trt.vs.ctrl1")

cd14ct_ls6_pos <- lsmeans(CD14_pos_mod1b ,~ Training|Time_point*Condition)
cd14ct_test6_pos <- contrast(cd14ct_ls6_pos,interaction = "trt.vs.ctrl1")

quartz()
ggplot(data=cbind(norm_cd14_counts,pred=fitted(CD14_pos_mod1b,type="response")))+
  facet_wrap(~Subject)+geom_smooth(aes(x=Time_point,y=pred),method="lm",alpha=0.3)+
  geom_point(aes(x=Time_point,y=norm_cd14,colour=Training))+
  geom_line(aes(x=Time_point,y=pred,colour=Training,group=Training))+
  ylab("Normalized proportions")+xlab("Time point")+ggtitle("CD14 normalised percentages")

##############CD16Hi##########
cd16hi_events <- combo_events[,-c(11,12,13)]
cbc_combo <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=7)
norm_cd16hi <- (cd16hi_events$`CD16hiGRpos`*cbc_combo[,4]*100)/(cd16hi_events$Leukocytes)
norm_cd16hi_counts <- cbind(cd16hi_events,norm_cd16hi)

CD16hi_pos_mod1 <- glmmPQL(norm_cd16hi ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                           Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Condition*Gender+
                           Training*Gender*Time_point+Condition*Gender*Time_point,random=~1|Subject,family=quasipoisson,data = norm_cd16hi_counts)

CD16hi_pos_mod2 <- glmmPQL(norm_cd16hi ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                           Training*Gender+Training*Condition+Gender*Condition,random=~1|Subject,family=quasipoisson,data = norm_cd16hi_counts)

CD16hi_pos_mod1b <- glmer.nb(norm_cd16hi ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                             Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Condition*Gender+
                             Training*Gender*Time_point+Condition*Gender*Time_point+(1|Subject)+
                             (1+Training|Subject),data=norm_cd16hi_counts)

CD16hi_pos_mod2b <- glmer.nb(norm_cd16hi ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                             Training*Gender+Training*Condition+Gender*Condition+(1|Subject)+
                             (1+Training|Subject),data=norm_cd16hi_counts)

cd16hict_ls1_pos <- lsmeans(CD16hi_pos_mod2b,~Time_point|Training)
cd16hict_test1_pos <- contrast(cd16hict_ls1_pos,interaction = "trt.vs.ctrl1")

cd16hict_ls2_pos <- lsmeans(CD16hi_pos_mod2b ,~Training|Time_point)
cd16hict_test2_pos <- contrast(cd16hict_ls2_pos,interaction = "trt.vs.ctrl1")

cd16hict_ls3_pos <- lsmeans(CD16hi_pos_mod1b ,~Time_point|Training*Condition)
cd16hict_test3_pos <- contrast(cd16hict_ls3_pos,interaction = "trt.vs.ctrl1")

cd16hict_ls4_pos <- lsmeans(CD16hi_pos_mod1b ,~Gender|Training*Time_point)
cd16hict_test4_pos <- contrast(cd16hict_ls4_pos,interaction = "trt.vs.ctrl1")

cd16hict_ls5_pos <- lsmeans(CD16hi_pos_mod1b ,~Condition|Training*Time_point)
cd16hict_test5_pos <- contrast(cd16hict_ls5_pos,interaction = "trt.vs.ctrl1")

cd16hict_ls6_pos <- lsmeans(CD16hi_pos_mod1b ,~ Training|Time_point*Condition)
cd16hict_test6_pos <- contrast(cd16hict_ls6_pos,interaction = "trt.vs.ctrl1")

#Facets with predicted values
quartz()
ggplot(data=cbind(norm_cd16hi_counts,pred=fitted(CD16hi_pos_mod1b)))+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=norm_cd16hi,colour=Training))+
  geom_line(aes(x=Time_point,y=pred,colour=Training,group=Training))+
  ylab("Normalized proportions")+xlab("Time point")+ggtitle("CD16Hi normalised percentages")


