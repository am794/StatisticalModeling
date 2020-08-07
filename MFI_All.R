############### Mixed-effects model for MFI ##############

library(openxlsx)
library(lme4)
library("nlme")
library(phia)
library("ggplot2")
library("reshape")
library(dplyr)
library(plyr)
library(reshape2)

############### Read the MFI table ################
mfi <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/MFI_all_panels.xlsx",sheet=5)
mfi_all <- as.data.frame(unclass(mfi))
mfi_wide <- dcast(mfi_all,Subject+Gender+Condition+Time_point+Training~Population,mean,value.var = "MFI")
#reshape(mfi_all,direction="wide",idvar = "Subject",timevar=Population)
# Fixed and random terms without interaction terms
All_mod1 <- lmer(MFI~Population+Condition+Gender+Time_point+Training+(1+Training|Subject),data=mfi_all)

# Two way interaction terms
All_mod2 <- lmer(MFI~Population+Condition+Gender+Time_point+Training+Population*Condition+
               Population*Gender+Population*Time_point+Population*Training+
               Condition*Training+Condition*Time_point+Condition*Gender+Gender*Time_point+
               Gender*Training+Time_point*Training+(1+Training|Subject),data=mfi_all)

# Three way interaction terms
All_mod3 <- update(mod2,.~ Population*Condition*Gender+Population*Condition*Time_point+
                 Population*Condition*Training+Population*Gender*Time_point+Population*Gender*Training+
                 Population*Time_point*Training+Condition*Gender*Time_point+Condition*Gender*Training+
                 Condition*Gender*Time_point+Gender*Training*Time_point+(1+Training|Subject))

###Post hoc tests
#Test 1: Pre-training vs post-training in baseline, peak, recovery
em1 <- emmeans(All_mod3,~Population|Time_point*Training)
contrast(em1,adjust="Tukey",interaction = "pairwise")

#Test 2: Baseline vs Peak/Recovery
em2 <- emmeans(All_mod3,~Time_point|Population*Training)
contrast(em2,interaction = "trt.vs.ctrl",adjust = "Tukey")

#Test 3: 
contrast(em1,"consec",simple="each",combine = TRUE,adjust="mvt",interaction="trt.vs.ctrl1")
em3 <- emmeans(mod3,~Training|Population+Time_point+Condition+Gender)
contrast(em3,interaction = "trt.vs.ctrl",adjust = "Tukey")


quartz()
ggplot(data=mfi_all)+geom_boxplot(aes(x=Time_point,y=MFI,fill=Training,colour=Training),width=0.4)+
  facet_wrap(~Population*Gender)

quartz()
ggplot(mfi_all)+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=MFI,colour=Training))+
  geom_line(data=mfi_all,aes(x=Time_point,y=MFI,colour=Training,group=Training))+
  ylab("MFI")+xlab("Time point")+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

quartz()
ggplot(mfi_all)+
  facet_wrap(.~Population+Condition+Gender)+
  geom_boxplot(data=cbind(mfi_all,pred=fitted(mod3)),aes(x=Time_point,y=pred,colour=Training),alpha=0.5)+
  ylab("MFI of GR+")+xlab("Time point")+ggtitle("Boxplot of MFI of GR+ across cell populations")+
scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme(strip.text = element_text(size=7,margin = margin(.06, 0, .06, 0, "cm")),
        strip.background = element_rect( fill = "#858585", color = NA ),    
        panel.background = element_rect( fill = "#efefef", color = NA ),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line( color = "#b2b2b2" ),
        panel.spacing.x = unit( 1, "cm" ),
        panel.spacing.y = unit( 0.5, "cm" ))

quartz()
ggplot(mfi_all)+
  facet_wrap(.~Population+Training)+
  geom_boxplot(data=cbind(mfi_all,pred=fitted(mod3)),aes(x=Time_point,y=pred),alpha=0.5)+
  ylab("MFI of GR+")+xlab("Time point")+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  ggtitle("Boxplot of MFI of GR+ across cell populations")+
  theme(strip.text = element_text(size=7,margin = margin(.06, 0, .06, 0, "cm")),
        strip.background = element_rect( fill = "#858585", color = NA ),    
        panel.background = element_rect( fill = "#efefef", color = NA ),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line( color = "#b2b2b2" ),
        panel.spacing.x = unit( 1, "cm" ),
        panel.spacing.y = unit( 0.5, "cm" ))

quartz()
emmip(mod3,Training+Condition~Time_point|Population+Gender)
monocyte_types <- mfi_all[with(mfi_all, Population %in%  c("Classical_Monocytes","Non-classical_monocytes","Intermediate_Monocytes")),]
monocyte_types <- as.data.frame(unclass(monocyte_types))
quartz()
ggplot(monocyte_types,aes(x=Population,y=MFI))+
  geom_bar(stat='summary',fun.y='mean',position=position_dodge(0.9),width = 0.4)+
  stat_summary(fun.data=mean_se,geom="errorbar",
               color="gray50",width=0.15,position = position_dodge(0.9))+
  scale_fill_manual(values=c("steelblue","orange"),guide = guide_legend(reverse = TRUE))+
  labs(x="Time point",y="normalized count for CD193(thousand/mcL)")

######### Bar plots for each time point figure1 from the paper #############
cell_types_baseline <-  mfi_all[with(mfi_all, Population %in%  c("GR+_CD193pos","GR+_CD203pos","GR+_CD14pos","GR+_CD16Hi","GR+_CD16pos","GR+_CD3pos") & Time_point %in% "baseline"),]
cell_types_baseline <- as.data.frame(unclass(cell_types_baseline))
quartz()
ggplot(cell_types_baseline,aes(x=Population,y=MFI))+
  geom_bar(stat='summary',fun.y='mean',position=position_dodge(0.9),width = 0.5)+
  stat_summary(fun.data=mean_se,geom="errorbar",
               color="gray50",width=0.15,position = position_dodge(0.9))+
  scale_fill_manual(values=c("steelblue","orange"),guide = guide_legend(reverse = TRUE))+
  labs(x="Cell Population",y="GR Expression(MFI)")+ggtitle("GR expression in Baseline")

cell_types_peak <-  mfi_all[with(mfi_all, Population %in%  c("GR+_CD193pos","GR+_CD203pos","GR+_CD14pos","GR+_CD16Hi","GR+_CD16pos","GR+_CD3pos") & Time_point %in% "peak"),]
cell_types_peak <- as.data.frame(unclass(cell_types_peak))
quartz()
ggplot(cell_types_peak,aes(x=Population,y=MFI))+
  geom_bar(stat='summary',fun.y='mean',position=position_dodge(0.9),width = 0.5)+
  stat_summary(fun.data=mean_se,geom="errorbar",
               color="gray50",width=0.15,position = position_dodge(0.9))+
  scale_fill_manual(values=c("steelblue","orange"),guide = guide_legend(reverse = TRUE))+
  labs(x="Cell Population",y="GR Expression(MFI)")+ggtitle("GR expression in peak")

cell_types_recovery <-  mfi_all[with(mfi_all, Population %in%  c("GR+_CD193pos","GR+_CD203pos","GR+_CD14pos","GR+_CD16Hi","GR+_CD16pos","GR+_CD3pos") & Time_point %in% "recovery"),]
cell_types_recovery <- as.data.frame(unclass(cell_types_recovery))
quartz()
ggplot(cell_types_recovery,aes(x=Population,y=MFI))+
  geom_bar(stat='summary',fun.y='mean',position=position_dodge(0.9),width = 0.5)+
  stat_summary(fun.data=mean_se,geom="errorbar",
               color="gray50",width=0.15,position = position_dodge(0.9))+
  scale_fill_manual(values=c("steelblue","orange"),guide = guide_legend(reverse = TRUE))+
  labs(x="Cell Population",y="GR Expression(MFI)")+ggtitle("GR expression in recovery")

