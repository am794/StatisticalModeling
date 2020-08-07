###############################
# UCI GLMM for FLOCK analysis #
###############################

#Load the packages
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

#Data
flock <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK_percentages_combo.xlsx",sheet = 1)
flock_cbc <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Oct11th/FLOCK_percentages_combo.xlsx",sheet = 2)

##NK cells
NKcells <- ((flock[,11]+flock[,16])*flock$GRpos*flock_cbc[,6]*10)/(flock$Leukocytes)
nkcells <- cbind(flock[,c(1:6)],NKcells)
nk_cells <- as.data.frame(unclass(nkcells))

#glmmTMB model
NK_mod1 <- glmmTMB(NKcells ~ Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Condition*Gender+
                Training*Gender*Time_point+Condition*Gender*Time_point+(1+Training|Subject)+
                  (1|Sample),family=poisson,data = nk_cells)

NK_mod2<- glmmTMB(NKcells ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                    Training*Gender+Training*Condition+Gender*Condition+(1+Training|Subject)+
                  (1|Sample),family=nbinom2,data = nk_cells)
#glmer neg binom model
mod1_NK <- glmer.nb(NKcells ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                  Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Condition*Gender+
                  Training*Gender*Time_point+Condition*Gender*Time_point+(1|Subject)+(1+Training|Subject),family= negative.binomial,data = nk_cells)

mod2_NK <- glmer.nb(NKcells ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                    Training*Gender+Training*Condition+Gender*Condition+(1|Subject)+(1+Training|Subject),data = nk_cells)

#Post-hoc tests
NK_ls1 <- lsmeans(mod2_NK,~Time_point|Training)
NK_test1 <- contrast(NK_ls1,interaction = "trt.vs.ctrl1")

NK_ls2 <- lsmeans(mod2_NK,~Training|Time_point)
NK_test2 <- contrast(NK_ls2,interaction = "trt.vs.ctrl1")

NK_ls3 <- lsmeans(mod1_NK,~Time_point|Training*Condition)
NK_test3 <- contrast(NK_ls3,interaction = "trt.vs.ctrl1")

NK_ls4 <- lsmeans(mod1_NK,~Gender|Training*Time_point)
NK_test4 <- contrast(NK_ls4,interaction = "trt.vs.ctrl1")

NK_ls5 <- lsmeans(mod1_NK,~Condition|Training*Time_point)
NK_test5 <- contrast(NK_ls5,interaction = "trt.vs.ctrl1")

NK_ls6 <- lsmeans(mod1_NK,~ Training|Time_point*Condition)
NK_test6 <- contrast(NK_ls6,interaction = "trt.vs.ctrl1")

#Post hoc and interaction plots
quartz()
ggplot(data=cbind(nk_cells,pred=fitted(mod1_NK,type="response")))+
  facet_wrap(~Subject)+geom_smooth(aes(x=Time_point,y=pred),method="lm",alpha=0.3)+
  geom_point(aes(x=Time_point,y=NKcells,colour=Training))+
  scale_color_manual(values=c( "#E69F00", "#56B4E9"))+
  geom_line(aes(x=Time_point,y=pred,colour=Training,group=Training))+
  ylab("Normalized counts")+xlab("Time point")+ggtitle("Line plots for NK cells")

summ_nk <- Rmisc::summarySE(nk_cells, measurevar="NKcells", groupvars=c("Training","Time_point"))
pd <- position_dodge(0.4)
quartz()
ggplot(summ_nk, aes(x=Time_point,y=NKcells,linetype=Training,group=Training))+ 
  geom_errorbar(aes(ymin=NKcells-se, ymax=NKcells+se),colour="black", 
                width=0.1,size=0.7,stat="identity", position = pd) +
  geom_line(position=pd,size=1) +
  geom_point(position=pd,size=3,shape=21)+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  ylim(c(0,1000))+theme_bw()+
  theme(legend.justification=c(1,0),
        legend.position=c(1,0)) +
  ylab("Normalized counts")+
  xlab("Time Point")+
  ggtitle("Mean-SE plot for normalized counts of NK cells")

quartz()
emmip(NK_mod1,Training~Condition|Time_point,CIs=TRUE)

quartz()
emmip(NK_mod1,Training~Gender|Time_point,CIs=TRUE)

quartz()
emmip(NK_mod1,Training~Time_point|Condition,CIs=TRUE)

quartz()
emmip(NK_mod1,Time_point~Training|Condition,CIs=TRUE)

quartz()
ggplot(nk_cells,aes(x=Time_point,y=NKcells,group=desc(Training),fill=Training))+
  geom_bar(stat='summary',fun.y='mean',position=position_dodge(0.9),width = 0.55)+
  stat_summary(fun.data=mean_se,geom="errorbar",
               color="gray50",width=0.15,position = position_dodge(0.9))+
  scale_fill_manual(values=c("steelblue","orange"),guide = guide_legend(reverse = TRUE))+
  labs(x="Time point",y="normalized count for NKcells")+
  ggtitle("Barplot for normalized counts for NK cells ")

#####CD3_flock
pop2_norm_flock <- flock[,8]*flock$GRpos*flock_cbc$absolute.lymphocyte.count*10/flock$Leukocytes
pop2_flock <- cbind(flock[,c(1:6)],pop2_norm_flock)
pop2flock <- as.data.frame(unclass(pop2_flock))

#glmmTMB model
pop2_mod1 <- glmmTMB(pop2_norm_flock ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                     Training*Gender+Training*Condition+Gender*Condition+Training*Condition*Time_point+Training*Condition*Gender+
                     Training*Gender*Time_point+Condition*Gender*Time_point+(1+Training|Subject)+
                     (1|Sample),family=poisson,data = pop2flock)

pop2_mod2<- glmmTMB(pop2_norm_flock ~Time_point+Training+Gender+Condition+Time_point*Training+Time_point*Gender+Time_point*Condition+
                    Training*Gender+Training*Condition+Gender*Condition+(1+Training|Subject)+
                    (1|Sample),family=poisson,data = pop2flock)

quartz()
ggplot(data=cbind(pop2flock,pred=fitted(pop2_mod1,type="response")))+
  facet_wrap(~Subject)+geom_smooth(aes(x=Time_point,y=pred),method="lm",alpha=0.3)+
  geom_point(aes(x=Time_point,y=NKcells,colour=Training))+
  scale_color_manual(values=c( "#E69F00", "#56B4E9"))+
  geom_line(aes(x=Time_point,y=pred,colour=Training,group=Training))+
  ylab("Normalized counts")+xlab("Time point")+ggtitle("Line plots for Population 2(FLOCK)")


summ_nk <- Rmisc::summarySE(pop2flock, measurevar="pop2_norm_flock", groupvars=c("Training","Time_point"))
pd <- position_dodge(0.4)
quartz()
ggplot(summ_nk, aes(x=Time_point,y=pop2_norm_flock,linetype=Training,group=Training))+ 
  geom_errorbar(aes(ymin=pop2_norm_flock-se, ymax=pop2_norm_flock+se),colour="black", 
                width=0.1,size=0.7,stat="identity", position = pd) +
  geom_line(position=pd,size=1) +
  geom_point(position=pd,size=3,shape=21)+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  ylim(c(0,2000))+theme_bw()+
  theme(legend.justification=c(1,0),
        legend.position=c(1,0)) +
  ylab("Normalized counts")+
  xlab("Time Point")+
  ggtitle("Mean-SE plot for normalized counts of Population#2 cells")

#Post-hoc tests
p2_ls1 <- lsmeans(pop2_mod2,~Time_point|Training)
p2_test1 <- contrast(p2_ls1,interaction = "trt.vs.ctrl1")

p2_ls2 <- lsmeans(pop2_mod2,~Training|Time_point)
p2_test2 <- contrast(p2_ls2,interaction = "trt.vs.ctrl1")

p2_ls3 <- lsmeans(pop2_mod1,~Time_point|Training*Condition)
p2_test3 <- contrast(p2_ls3,interaction = "trt.vs.ctrl1")

p2_ls4 <- lsmeans(pop2_mod1,~Gender|Training*Time_point)
p2_test4 <- contrast(p2_ls4,interaction = "trt.vs.ctrl1")

p2_ls5 <- lsmeans(pop2_mod1,~Condition|Training*Time_point)
p2_test5 <- contrast(p2_ls5,interaction = "trt.vs.ctrl1")

p2_ls6 <- lsmeans(pop2_mod1,~ Training|Time_point*Condition)
p2_test6 <- contrast(p2_ls6,interaction = "trt.vs.ctrl1")

quartz()
emmip(pop2_mod1,Training~Time_point|Condition,CIs=TRUE)
