
############### Generalized linear mixed effects model for normalized percentages ##############

#load the packages
library(openxlsx)
library(lme4)
library("nlme")
library(phia)
library("ggplot")
library("reshape")
library(dplyr)
library(plyr)
library(emmeans)

################################### CD193 normalized percentages ##########################################
CD193_events <-  read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=3)
CD193_events <- as.data.frame(unclass(CD193_events))
cbc_193_203 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=5)

#model1 subject level variation model2 observational variation
cd193_mod1 <- glmer((`CD193posGRpos`*cbc_193_203[,6]/100)/(Leukocytes)~Time_point*Training*Gender*Condition+(1|Subject),weights = Leukocytes,family=binomial,data=CD193_events)
cd193_mod2 <- glmer((`CD193posGRpos`*cbc_193_203[,6]/100)/(Leukocytes)~Time_point*Training*Gender*Condition+(1+Training|Subject)+(1|Sample),weights = Leukocytes,family=binomial,data=CD193_events)
cd193_mod3 <- glmer((`CD193posGRpos`*cbc_193_203[,6]/100)/(Leukocytes)~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                      Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Condition*Gender+
                      Training*Gender*Time_point+Condition*Gender*Time_point+
                      (1+Training|Subject)+(1|Sample),weights = Leukocytes,family=binomial,data=CD193_events)
#cd193_mod3/mod4 is significant and mod2 with four way interaction is not significant
cd193_mod4 <- glmer((`CD193posGRpos`*cbc_193_203[,6]/100)/(Leukocytes)~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                      Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Condition*Gender+
                      Training*Gender*Time_point+Condition*Gender*Time_point+
                      (1+Training|Subject)+(1|Sample),weights = Leukocytes,family=binomial,data=CD193_events)

cd193_mod7 <- glmer((`CD193posGRpos`*cbc_193_203[,6]/100)/(Leukocytes)~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                      Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+
                      (1+Training|Subject)+(1|Sample),weights = Leukocytes,family=binomial,data=CD193_events)

cd193_mod6 <- glmer((`CD193posGRpos`*cbc_193_203[,6]/100)/(Leukocytes)~Time_point+Training+Gender+Condition+
                      (1|Subject)+(1|Sample),weights = Leukocytes,family=binomial,data=CD193_events)

#Post-hoc tests
#Test1: Pre vs post in baseline peak and recovery
cd193_ls1 <- lsmeans(cd193_mod3,~Training|Time_point)
cd193_test1 <- contrast(cd193_ls1,interaction = "trt.vs.ctrl1")

#Test2: Baseline vs Peak and recovery
cd193_ls2 <- emmeans(cd193_mod3,~Time_point|Condition*Training)
cd193_test2 <- contrast(cd193_ls2,alpha=0.05,interaction = "trt.vs.ctrl1")

#Test 3: Pre vs post in healthy/recovery and baseline/peak/recovery
cd193_ls3 <- lsmeans(cd193_mod3,~Training|Condition*Time_point)
cd193_test3 <- contrast(cd193_ls3,interaction = "trt.vs.ctrl1")

#Test 4: Healthy vs Asthmatic in pre/post and baseline/peak/recovery
cd193_ls4 <- lsmeans(cd193_mod3,~Condition|Training*Time_point)
cd193_test4 <- contrast(cd193_ls4,interaction = "trt.vs.ctrl1")

#Test 5: Male vs Female in pre/post and baseline/peak/recovery
cd193_ls5 <- lsmeans(cd193_mod3,~Gender|Training*Time_point)
cd193_test5 <- contrast(cd193_ls5,interaction = "trt.vs.ctrl1")

#Facets with predicted values line and actual value dot
quartz()
ggplot(CD193_events)+
  facet_wrap(~Condition+Gender)+
  geom_point(aes(x=Time_point,y=((CD193posGRpos*cbc_193_203[,6]/100)/(Leukocytes)),colour=Training))+
  geom_boxplot(data=cbind(CD193_events,pred=fitted(mod_perc1)),aes(x=Time_point,y=pred,colour=Training))+
  ylab("Normalized proportions")+xlab("Time point")

#Facets with predicted values line and actual value dot
quartz()
ggplot(CD193_events)+
  facet_wrap(~Time_point+Condition+Gender)+
  geom_point(aes(x=Training,y=((CD193posGRpos*cbc_193_203[,6]/100)/(Leukocytes)),colour=Training))+
  geom_boxplot(data=cbind(CD193_events,pred=fitted(mod_perc1)),aes(x=Training,y=pred,colour=Training))+
  ylab("Normalized proportions")+xlab("Time point")

#Facets with predicted values
quartz()
ggplot(data=cbind(CD193_events,pred=fitted(cd193_mod4)*100))+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=pred,colour=Training))+
  geom_line(aes(x=Time_point,y=pred,colour=Training,group=Training))+
  ylab("Normalized proportions")+xlab("Time point")+ggtitle("CD193 normalised percentages")

#Facets with actual values
quartz()
ggplot(CD193_events)+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=((CD193posGRpos*cbc_193_203[,6])*100/(Leukocytes)),colour=Training))+
  geom_line(data=CD193_events,aes(x=Time_point,y=((CD193posGRpos*cbc_193_203[,6]*100)/(Leukocytes)),colour=Training,group=Training))+
  ylab("Normalized percentages")+xlab("Time point")+ggtitle("CD193 normalized counts(thous/mcL)")

#Boxplot with actual values
quartz()
ggplot(CD193_events)+
  geom_boxplot(aes(x=Time_point,y=((CD193posGRpos*cbc_193_203[,6]/100)/(Leukocytes)),colour=Training))+
  facet_wrap(~Condition)

#Boxplot for condition with actual values
quartz()
ggplot(CD193_events)+
  geom_boxplot(aes(x=Time_point,y=((CD193posGRpos*cbc_193_203[,6])/(Leukocytes)),colour=Condition),width=0.3)+
  facet_wrap(~Training)+xlab("Time point")+ylab("GR+ CD3+ normalized percentages")+
  ggtitle("GR+ CD3+ normalized percentages across the time points and training")+
  scale_color_manual(values=c("#999999", "#E69F00"))+theme_test()

#Boxplot for condition with predicted values
quartz()
ggplot(CD193_events)+
  geom_boxplot(data=cbind(CD193_events,pred=fitted(cd193_mod4)),
               aes(x=Time_point,y=pred,colour=Condition),width=0.3)+
  facet_wrap(~Training)+xlab("Time point")+ylab("GR+ CD3+ normalized percentages")+
  ggtitle("GR+ CD3+ normalized percentages across the time points and training")+
  scale_color_manual(values=c("#999999", "#E69F00"))+theme_test()

quartz()
ggplot(CD193_events)+
  geom_boxplot(data=cbind(CD193_events,pred=fitted(cd193_mod2)),
               aes(x=Time_point,y=pred,colour=Training),width=0.3)+
  facet_wrap(~Condition+Gender)+xlab("Time point")+ylab("GR+ CD3+ normalized percentages")+
  ggtitle("GR+ CD3+ normalized percentages across the time points and training")+
  scale_color_manual(values=c("#999999", "#E69F00"))+theme_test()

#Visulaize the model
#library(effects)
eff <- allEffects(cd193_mod2)
eff.df <- as.data.frame(eff[[1]])
plot(eff,multiline=TRUE,confint=TRUE,ci.style="bars")

################################### CD203 normalized percentages ##########################################
CD203_events <-  read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=4)
cbc_193_203 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=5)

#model1 includes subject level variation and level 2 includes observational level variation
cd203_mod1 <- glmer((`CD203pGRp`*cbc_193_203[,5]/100)/(Leukocytes)~Time_point*Training*Gender*Condition+(1|Subject),weights = Leukocytes,family=binomial,data=CD203_events)
cd203_mod2 <- glmer((`CD203pGRp`*cbc_193_203[,5]/100)/(Leukocytes)~Time_point*Training*Gender*Condition+(1|Subject)+(1|Sample),weights = Leukocytes,family=binomial,data=CD203_events)
##mod 3 is significant i.e the four way interaction term is significant
cd203_mod3 <- glmer((`CD203pGRp`*cbc_193_203[,5]/100)/(Leukocytes)~Time_point*Training*Gender*Condition+(1+Training|Subject)+(1|Sample),weights = Leukocytes,family=binomial,data=CD203_events)
cd203_mod4 <- glmer((`CD203pGRp`*cbc_193_203[,5]/100)/(Leukocytes)~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                      Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Condition*Gender+
                      Training*Gender*Time_point+Condition*Gender*Time_point+
                      (1+Training|Subject)+(1|Sample),weights = Leukocytes,family=binomial,data=CD203_events)


#Post-hoc tests
#Test1: Pre vs post in baseline peak and recovery
cd203_ls1 <- lsmeans(cd203_mod3,~Training|Time_point)
cd203_test1 <- contrast(cd203_ls1,interaction = "trt.vs.ctrl1")

#Test2: Baseline vs Peak and recovery
cd203_ls2 <- lsmeans(cd203_mod3,~Time_point|Condition*Training)
cd203_test2 <- contrast(cd203_ls2,alpha=0.05,interaction = "trt.vs.ctrl1")

#Test 3: Pre vs post in healthy/recovery and baseline/peak/recovery
cd203_ls3 <- lsmeans(cd203_mod3,~Training|Condition*Time_point)
cd203_test3 <- contrast(cd203_ls3,interaction = "trt.vs.ctrl1")

#Test 4: Healthy vs Asthmatic in pre/post and baseline/peak/recovery
cd203_ls4 <- lsmeans(cd203_mod3,~Condition|Training*Time_point)
cd203_test4 <- contrast(cd203_ls4,interaction = "trt.vs.ctrl1")

#Test 5: Male vs Female in pre/post and baseline/peak/recovery
cd203_ls5 <- lsmeans(cd203_mod3,~Gender|Training*Time_point)
cd203_test5 <- contrast(cd203_ls5,interaction = "trt.vs.ctrl1")

quartz()
ggplot(CD203_events)+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=((`CD203pGRp`*cbc_193_203[,5]*100)/(Leukocytes)),colour=Training))+
  geom_line(data=CD203_events,aes(x=Time_point,y=((`CD203pGRp`*cbc_193_203[,5]*100)/(Leukocytes)),colour=Training,group=Training))+
  ylab("Normalized percentages")+xlab("Time point")+ggtitle("CD203 normalized counts(thous/mcL)")

############# Read the events table ###############
CD3_event <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=1)
cbc_cd3 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=2)
CD3_events <- as.data.frame(unclass(CD3_event))

combo_events <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=6)
cbc_combo <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=7)
combo_events <- as.data.frame(unclass(combo_events))

cd3_mod1 <- glmer((`CD3posGRpos`*cbc_cd3[,4]/100)/(Leukocytes)~Time_point*Training*Gender*Condition+(1|Subject),weights = Leukocytes,family = binomial,data = CD3_events)
cd3_mod2 <- glmer((`CD3posGRpos`*cbc_cd3[,4]/100)/(Leukocytes)~Time_point*Training*Gender*Condition+(1|Subject)+(1|Sample),weights = Leukocytes,family = binomial,data = CD3_events)
#cd3_mod3 is more significant i.e. four way interaction term
cd3_mod3 <- glmer((`CD3posGRpos`*cbc_cd3[,4]/100)/(Leukocytes)~Time_point*Training*Gender*Condition+(1+Training|Subject)+(1|Sample),weights = Leukocytes,family = binomial,data = CD3_events)
cd3_mod4 <- glmer((`CD3posGRpos`*cbc_cd3[,4]/100)/(Leukocytes)~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                      Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+
                      (1+Training|Subject)+(1|Sample),weights = Leukocytes,family=binomial,data=CD3_events)

#Post-hoc tests
#Test1: Pre vs post in baseline peak and recovery
cd3_ls1 <- lsmeans(cd3_mod3,~Training|Time_point)
cd3_test1 <- contrast(cd3_ls1,interaction = "trt.vs.ctrl1")

#Test2: Baseline vs Peak and recovery
cd3_ls2 <- lsmeans(cd3_mod3,~Time_point|Condition*Training)
cd3_test2 <- contrast(cd3_ls2,alpha=0.05,interaction = "trt.vs.ctrl1")

#Test 3: Pre vs post in healthy/recovery and baseline/peak/recovery
cd3_ls3 <- lsmeans(cd3_mod3,~Training|Condition*Time_point)
cd3_test3 <- contrast(cd3_ls3,interaction = "trt.vs.ctrl1")

#Test 4: Healthy vs Asthmatic in pre/post and baseline/peak/recovery
cd3_ls4 <- lsmeans(cd3_mod3,~Condition|Training*Time_point)
cd3_test4 <- contrast(cd3_ls4,interaction = "trt.vs.ctrl1")

#Test 5: Male vs Female in pre/post and baseline/peak/recovery
cd3_ls5 <- lsmeans(cd3_mod3,~Gender|Training*Time_point)
cd3_test5 <- contrast(cd3_ls5,interaction = "trt.vs.ctrl1")

quartz()
ggplot(CD3_events)+
  geom_boxplot(data=cbind(CD3_events,pred=fitted(cd3_mod3)),
               aes(x=Time_point,y=pred,colour=Condition),width=0.3)+
  facet_wrap(~Training)+xlab("Time point")+ylab("GR+ CD3+ normalized percentages")+
  ggtitle("GR+ CD3+ normalized percentages across the time points and training")+
  scale_color_manual(values=c("#999999", "#E69F00"))+theme_test()

quartz()
ggplot(CD3_events)+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=((`CD3posGRpos`*cbc_cd3[,4]*100)/(Leukocytes)),colour=Training))+
  geom_line(data=CD3_events,aes(x=Time_point,y=((`CD3posGRpos`*cbc_cd3[,4]*100)/(Leukocytes)),colour=Training,group=Training))+
  ylab("Normalized percentages(thous/mcL)")+xlab("Time point")+ggtitle("CD3 normalized counts(thous/mcL)")

############CD14###########
cd14_events <- combo_events[,-c(10,11,13)]
cbc_combo <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=7)

cd14_mod1 <- glmer((`CD14posGRpos`*cbc_combo[,5]/100)/(Leukocytes)~Time_point*Training*Gender*Condition+(1|Subject),weights = Leukocytes,family = binomial,data = cd14_events)
cd14_mod2 <- glmer((`CD14posGRpos`*cbc_combo[,5]/100)/(Leukocytes)~Time_point*Training*Gender*Condition+(1|Subject)+(1|Sample),weights = Leukocytes,family = binomial,data = cd14_events)
cd14_mod3 <- glmer((`CD14posGRpos`*cbc_combo[,5]/100)/(Leukocytes)~Time_point*Training*Gender*Condition+(1+Training|Subject)+(1|Sample),weights = Leukocytes,family = binomial,data = cd14_events)
cd14_mod4 <- glmer((`CD14posGRpos`*cbc_combo[,5]/100)/(Leukocytes)~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                    Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+
                    (1+Training|Subject)+(1|Sample),weights = Leukocytes,family=binomial,data=cd14_events)

#Post-hoc tests
#Test1: Pre vs post in baseline peak and recovery
cd14_ls1 <- lsmeans(cd14_mod4,~Training|Time_point)
cd14_test1 <- contrast(cd14_ls1,interaction = "trt.vs.ctrl1")

#Test2: Baseline vs Peak and recovery
cd14_ls2 <- lsmeans(cd14_mod4,~Time_point|Condition*Training)
cd14_test2 <- contrast(cd14_ls2,alpha=0.05,interaction = "trt.vs.ctrl1")

#Test 3: Pre vs post in healthy/recovery and baseline/peak/recovery
cd14_ls3 <- lsmeans(cd14_mod4,~Training|Condition*Time_point)
cd14_test3 <- contrast(cd14_ls3,interaction = "trt.vs.ctrl1")

#Test 4: Healthy vs Asthmatic in pre/post and baseline/peak/recovery
cd14_ls4 <- lsmeans(cd14_mod4,~Condition|Training*Time_point)
cd14_test4 <- contrast(cd14_ls4,interaction = "trt.vs.ctrl1")

#Test 5: Male vs Female in pre/post and baseline/peak/recovery
cd14_ls5 <- lsmeans(cd14_mod4,~Gender|Training*Time_point)
cd14_test5 <- contrast(cd14_ls5,interaction = "trt.vs.ctrl1")

quartz()
ggplot(cd14_events)+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=(`CD14posGRpos`*cbc_combo[,5]/100)/(Leukocytes),colour=Training),alpha=0.6)+
  scale_color_manual(values=c("steelblue","orange"),guide = guide_legend(reverse = TRUE))+scale_linetype_manual(values=c("twodash", "dotted"),guide = FALSE)+
  geom_line(data=CD3_events,aes(x=Time_point,y=((`CD14posGRpos`*cbc_combo[,5]/100)/(Leukocytes)),colour=Training,group=Training,linetype=Training),size=0.7)+
  ylab("Normalized counts(thous/mcL)")+xlab("Time point")+ggtitle("CD14 normalized counts across basline, peak and recovery")


##############CD16Hi##########
cd16hi_events <- combo_events[,-c(11,12,13)]
cbc_combo <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=7)
cd16hi_ev_mod1 <- glmer((`CD16hiGRpos`*cbc_combo[,4]/100)/(Leukocytes)~Time_point*Training*Gender*Condition+(1|Subject),weights = Leukocytes,family = binomial,data = cd16hi_events)
cd16hi_ev_mod2 <- glmer((`CD16hiGRpos`*cbc_combo[,4]/100)/(Leukocytes)~Time_point*Training*Gender*Condition+(1|Subject)+(1|Sample),weights = Leukocytes,family = binomial,data = cd16hi_events)
#mod3 is significant i.e. four way interaction
cd16hi_ev_mod3 <- glmer((`CD16hiGRpos`*cbc_combo[,4]/100)/(Leukocytes)~Time_point*Training*Gender*Condition+(1+Training|Subject)+(1|Sample),weights = Leukocytes,family = binomial,data = cd16hi_events)
cd16hi_ev_mod4 <- glmer((`CD16hiGRpos`*cbc_combo[,4]/100)/(Leukocytes)~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                     Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+
                     (1+Training|Subject)+(1|Sample),weights = Leukocytes,family=binomial,data=cd16hi_events)

pairs(cd16hi_ls2,by=NULL,adjust="Dunnett")
#Post-hoc tests
em1 <- emmeans(cd16hi_ev_mod3,~Time_point|Training*Condition*Gender)
contrast(em1,"consec",simple="each",combine = TRUE,adjust="mvt",interaction="trt.vs.ctrl1")
emmip(cd16hi_ev_mod3,Training~Time_point|Condition*Gender)
contrast(em1,interaction = "pairwise")
em1 <- emmeans(cd16hi_ev_mod3,~Time_point*Training*Condition)
contrast(em1,"consec",simple="each",combine = TRUE,adjust="mvt",interaction="trt.vs.ctrl1")
pairs(em1)
pairs(em1,interaction="trt.vs.ctrl1",adjust="mvt")

#Test1: Pre vs post in baseline peak and recovery
cd16hi_ls1 <- emmeans(cd16hi_ev_mod3,~Training|Time_point)
cd16hi_test1 <- contrast(cd16hi_ls1,interaction = "trt.vs.ctrl1")
ref_dfR <- as.data.frame(summary(cd16hi_ls1))
pd <- position_dodge(0.1)
quartz()
ggplot(ref_dfR,aes(Time_point,y=emmean,group=Training,colour=Training))+
  geom_errorbar(aes(ymin=emmean-SE,ymax=emmean+SE),width=0.1,position = pd)+
  geom_line(position = pd)+geom_point(position=pd)

#Test2: Baseline vs Peak and recovery
cd16hi_ls2 <- lsmeans(cd16hi_ev_mod3,~Time_point|Condition*Training)
cd16hi_test2 <- contrast(cd16hi_ls2,interaction = "trt.vs.ctrl1")

#Test 3: Pre vs post in healthy/recovery and baseline/peak/recovery
cd16hi_ls3 <- lsmeans(cd16hi_ev_mod3,~Training|Condition*Time_point)
cd16hi_test3 <- contrast(cd16hi_ls3,interaction = "trt.vs.ctrl1")

#Test 4: Healthy vs Asthmatic in pre/post and baseline/peak/recovery
cd16hi_ls4 <- lsmeans(cd16hi_ev_mod3,~Condition|Training*Time_point)
cd16hi_test4 <- contrast(cd16hi_ls4,interaction = "trt.vs.ctrl1")

#Test 5: Male vs Female in pre/post and baseline/peak/recovery
cd16hi_ls5 <- lsmeans(cd16hi_ev_mod3,~Gender|Training*Time_point)
cd16hi_test5 <- contrast(cd16hi_ls5,interaction = "trt.vs.ctrl1")

quartz()
ggplot(cd16hi_events)+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=((`CD16hiGRpos`*cbc_combo[,4]/100)/(Leukocytes)),colour=Training),alpha=0.6)+
  scale_color_manual(values=c("steelblue","orange"),guide = guide_legend(reverse = TRUE))+scale_linetype_manual(values=c("twodash", "dotted"),guide = FALSE)+
  geom_line(data=CD3_events,aes(x=Time_point,y=((`CD16hiGRpos`*cbc_combo[,4]/100)/(Leukocytes)),colour=Training,group=Training,linetype=Training),size=0.7)+
  ylab("Normalized counts(thous/mcL)")+xlab("Time point")+ggtitle("CD16hi normalized counts across basline, peak and recovery")

agg <- ddply(CD193_events,.(Condition,Training,Time_point,Gender),function(x) {c(mean=mean(x$CD193_percentage),sd=sd(x$CD193_percentage))})
agg$lower = agg$mean + agg$sd
agg$upper = agg$mean - agg$sd
quartz()
ggplot(agg,aes(x=Time_point,y=mean,colour=Condition))+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=0.3,position = pd)+
  geom_point(position = pd)+facet_grid(Gender~Training)+geom_line(position = pd)

#########CD16 subset 10:
flock_subset <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK_percentages.xlsx",sheet=1)
cbc_cd3_flock <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK_percentages.xlsx",sheet=2)
flock_subsets <- as.data.frame(unclass(flock_subset))
flock_perc <- flock_subsets$Pop10*flock_subsets$GRpos
cd16_pop10_mod1 <- glmer((flock_perc*cbc_cd3_flock[,4]/100)/(Leukocytes)~Time_point*Training*Gender*Condition+(1+Training|Subject)+(1|Sample),weights = Leukocytes,family = binomial,data = flock_subsets)
cd16_pop10_mod2 <- glmer((flock_perc*cbc_cd3_flock[,4]/100)/(Leukocytes)~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Gender*Condition+Time_point*Gender*Condition+Training*Time_point*Gender+
                          (1+Training|Subject)+(1|Sample),weights = Leukocytes,family=binomial,data=flock_subsets)



####### Fold change #########
cd14_percs <- cbind(cd14_events[,c(2,3,4)],(cd14_events$`CD14posGRpos`*cbc_combo[,4]/100)/(cd14_events$Leukocytes))
colnames(cd14_percs)[4] <- "cd14_percentage"
dcast(cd14_percs,Training*Subject~Time_point)

############################ Preliminary analysis
#Model1: Subject level random factors
model_cd3_1 <- glmer((`CD3posGRpos`*cbc_cd3[,4]/100)/(Leukocytes)~Time_point*Training*Gender*Condition+(1|Subject),weights = Leukocytes,family = binomial,data = CD3_events)
ls2 <- lsmeans(model_cd3_1,adjust="Tukey",specs="Time_point","Training")
t2 <- contrast(ls2)
t2 <- contrast(ls2,interaction = "pairwise")
mod1_means <- interactionMeans(model_cd3_1)

#Model2: Subject and observational level random factors
model_cd3_2 <- glmer((`CD3posGRpos`*cbc_cd3[,4]/100)/(Leukocytes)~Time_point*Training*Gender*Condition+(1|Subject)+(1|Sample),weights = Leukocytes,family = binomial,data = CD3_events)
ls3 <- lsmeans(model_cd3_2,~Time_point|Training,adjust="Tukey")
t3 <- contrast(ls3,interaction = 'trt.vs.ctrl1')
t2 <- contrast(ls2,interaction = "pairwise")
plot(mod1_means)
#m_nb <- glmer.nb((`CD3posGRpos`*cbc_cd3[,4])/(Leukocytes)~Time_point*Training*Gender*Condition+(1|Subject),data=CD3_events)
mod1 <- glmer((`CD193posGRpos`*cbc_193_203[,6]/100)/(Leukocytes)~Time_point*Training*Gender*Condition+(1|Subject)+(1|Sample),weights = Leukocytes,family=binomial,data=CD193_events)


#Overdispersion function
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2    
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")  
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

ls <- lsmeans(model_cd3_2,~Time_point|Condition*Training)
contrast(ls,interaction = "trt.vs.ctrl1",adjust="Tukey")
contrast(ls,interaction = "pairwise",adjust="Bonferroni")

CD3_events %>% mutate(pred_dist=fitted(model_cd3_2)) %>% ggplot(aes(x=Time_point,y=pred_dist,group=Subject,color=Subject))+geom_line(size=1)
quartz()
CD3_events %>% mutate(pred_perc=fitted(model_cd3_2)) %>% ggplot(aes(x=Time_point,y=pred_perc,color=Training))+
  geom_point()+geom_smooth(method="lm")+facet_wrap(~Subject)+theme_bw() +
  theme(panel.grid = element_blank())

# Facets with predicted values
quartz()
ggplot(CD3_events)+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=((CD3posGRpos*cbc_cd3[, 4]/100)/(Leukocytes)),colour=Training))+
  geom_line(data=cbind(CD3_events,pred=fitted(model_cd3_2)),aes(x=Time_point,y=pred,colour=Training,group=Training))+
  ylab("Normalized proportions")+xlab("Time point")

# Facets with actual values
quartz()
ggplot(CD3_events)+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=((CD3posGRpos*cbc_cd3[, 4])/(Leukocytes)),colour=Training))+
  geom_line(data=CD3_events,aes(x=Time_point,y=((CD3posGRpos*cbc_cd3[, 4])/(Leukocytes)),colour=Training,group=Training))+
  ylab("Normalized percentages")+xlab("Time point")

#Boxplot with actual values
quartz()
ggplot(CD3_events)+
  geom_boxplot(aes(x=Time_point,y=((CD3posGRpos*cbc_cd3[, 4]/100)/(Leukocytes)),colour=Training))+
  facet_wrap(~Condition)

#Boxplot for condition with actual values
quartz()
ggplot(CD3_events)+
  geom_boxplot(aes(x=Time_point,y=((CD3posGRpos*cbc_cd3[, 4])/(Leukocytes)),colour=Condition),width=0.3)+
  facet_wrap(~Training)+xlab("Time point")+ylab("GR+ CD3+ normalized percentages")+
  ggtitle("GR+ CD3+ normalized percentages across the time points and training")+
  scale_color_manual(values=c("#999999", "#E69F00"))+theme_test()

#Boxplot with predicted values
quartz()
ggplot(CD3_events)+
  geom_boxplot(data=cbind(CD3_events,pred=fitted(model_cd3_2)),
               aes(x=Time_point,y=pred,colour=Condition),width=0.3)+
  facet_wrap(~Training)+xlab("Time point")+ylab("GR+ CD3+ normalized percentages")+
  ggtitle("GR+ CD3+ normalized percentages across the time points and training")+
  scale_color_manual(values=c("#999999", "#E69F00"))+theme_test()

quartz()
ggplot(CD3_events)+
  geom_boxplot(aes(x=Time_point,y=((CD3posGRpos*cbc_cd3[, 4])/(Leukocytes)),colour=Condition),width=0.3)+
  facet_wrap(~Training)+xlab("Time point")+ylab("GR+ CD3+ normalized percentages")+
  ggtitle("GR+ CD3+ normalized percentages across the time points and training")+
  scale_color_manual(values=c("#999999", "#E69F00"))+theme_test()

quartz()
ggplot(CD3_events)+
  geom_boxplot(aes(x=Condition,y=((CD3posGRpos*cbc_cd3[,4])/(Leukocytes)),colour=Training),width=0.3)+
  facet_wrap(~Time_point)+xlab("Condition")+ylab("GR+ CD3+ normalized percentages")+
  ggtitle("GR+ CD3+ normalized percentages across the time points and training")+
  scale_color_manual(values=c("#999999", "#E69F00"))+theme_test()

quartz()
ggplot(CD3_events)+
  geom_boxplot(data=cbind(CD3_events,pred=fitted(model_cd3_2)),
               aes(x=Condition,y=pred,colour=Training),width=0.3)+
  facet_wrap(~Time_point)+xlab("Condition")+ylab("GR+ CD3+ normalized percentages")+
  ggtitle("GR+ CD3+ normalized percentages across the time points and training")+
  scale_color_manual(values=c("#999999", "#E69F00"))+theme_test()

ls2 <- lsmeans(model_cd3_2,~Time_point|Training)
contrast(ls2,interaction = "trt.vs.ctrl1",adjust="Bonferroni")
contrast(ls2,interaction = "pairwise",adjust="Bonferroni")


ls2 <- lsmeans(model_cd3_2,~Training|Time_point)
contrast(ls2,interaction = "trt.vs.ctrl1",adjust="Bonferroni")
contrast(ls2,interaction = "pairwise",adjust="Bonferroni")

quartz()
ggplot(CD3_events)+
  geom_boxplot(aes(x=Time_point,y=((CD3posGRpos*cbc_cd3[,4]/100)/(Leukocytes))))+
  facet_wrap(~Training)


ls3 <- lsmeans(model_cd3_2,~Training|Time_point*Condition)
contrast(ls3,interaction = "trt.vs.ctrl1",adjust="Bonferroni")
contrast(ls3,interaction = "pairwise",adjust="Bonferroni")


##MFI LMM
mfi_CD3 <- mfi_combo[,c(1,2,3,4,5,6,7,11)]
mfi_model <- lmer(GRpos_CD3pos~Time_point*Training*Gender*Condition+(1|Subject),data = mfi_CD3)
ls <- lsmeans(mfi_model,~Time_point|Condition*Training)
contrast(ls,interaction = "trt.vs.ctrl1",adjust="Bonferroni")
contrast(ls,interaction = "pairwise",adjust="Bonferroni")
quartz()
mfi_combo %>% mutate(pred_perc=fitted(mfi_model)) %>% ggplot(aes(x=Time_point,y=pred_perc,color=Training))+
  geom_point()+geom_smooth(method="lm")+facet_wrap(~Subject)+theme_bw() +
  theme(panel.grid = element_blank())


quartz()
ggplot(mfi_combo)+
  geom_boxplot(aes(x=Time_point,y=GRpos_CD3pos,colour=Training))+
  facet_wrap(~Condition)

quartz()
ggplot(mfi_combo)+
  geom_boxplot(aes(x=Time_point,y=GRpos_CD3pos,colour=Condition))+
  facet_wrap(~Training)
ls <- lsmeans(mfi_model,~Time_point|Training)
contrast(ls,interaction = "pairwise",adjust="Bonferroni")


> K <- cbind(0,diag(length(fixef(cd3_mod3))-1)) 
> rownames(K) <- names(fixef(cd3_mod3))[-1] 
> model_glht <- glht(cd3_mod3,linfct=K) 
> summary(model_glht)


