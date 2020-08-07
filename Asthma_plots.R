######## Plots ########
library(ggplot2)

##########Subject wise normalized counts(actual values)
quartz()
ggplot(CD193_events)+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=((CD193posGRpos*cbc_193_203[,6])/(Leukocytes)),colour=Training),alpha=0.6)+
  scale_color_manual(values=c("steelblue","orange"),guide = guide_legend(reverse = TRUE))+scale_linetype_manual(values=c("twodash", "dotted"),guide = FALSE)+
  geom_line(data=CD193_events,aes(x=Time_point,y=((CD193posGRpos*cbc_193_203[,6])/(Leukocytes)),colour=Training,group=Training,linetype=Training),size=0.7)+
  ylab("Normalized counts(thous/mcL)")+xlab("Time point")+ggtitle("CD193 normalized counts across basline, peak and recovery")

quartz()
ggplot(CD203_events)+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=((`CD203pGRp`*cbc_193_203[,5])/(Leukocytes)),colour=Training),alpha=0.6)+
  scale_color_manual(values=c("steelblue","orange"),guide = guide_legend(reverse = TRUE))+scale_linetype_manual(values=c("twodash", "dotted"),guide = FALSE)+
  geom_line(data=CD203_events,aes(x=Time_point,y=((`CD203pGRp`*cbc_193_203[,5])/(Leukocytes)),colour=Training,group=Training,linetype=Training),size=0.7)+
  ylab("Normalized counts(thous/mcL)")+xlab("Time point")+ggtitle("CD203 normalized counts across basline, peak and recovery")

quartz()
ggplot(CD3_events)+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=((`CD3posGRpos`*cbc_cd3[,4])/(Leukocytes)),colour=Training),alpha=0.6)+
  scale_color_manual(values=c("steelblue","orange"),guide = guide_legend(reverse = TRUE))+scale_linetype_manual(values=c("twodash", "dotted"),guide = FALSE)+
  geom_line(data=CD3_events,aes(x=Time_point,y=((`CD3posGRpos`*cbc_cd3[,4])/(Leukocytes)),colour=Training,group=Training,linetype=Training),size=0.7)+
  ylab("Normalized counts(thous/mcL)")+xlab("Time point")+ggtitle("CD3 normalized counts across basline, peak and recovery")

quartz()
ggplot(cd14_events)+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=(`CD14posGRpos`*cbc_combo[,5]/100)/(Leukocytes),colour=Training),alpha=0.6)+
  scale_color_manual(values=c("steelblue","orange"),guide = guide_legend(reverse = TRUE))+scale_linetype_manual(values=c("twodash", "dotted"),guide = FALSE)+
  geom_line(data=CD3_events,aes(x=Time_point,y=((`CD14posGRpos`*cbc_combo[,5]/100)/(Leukocytes)),colour=Training,group=Training,linetype=Training),size=0.7)+
  ylab("Normalized counts(thous/mcL)")+xlab("Time point")+ggtitle("CD14 normalized counts across basline, peak and recovery")

quartz()
ggplot(cd16hi_events)+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=((`CD16hiGRpos`*cbc_combo[,4]/100)/(Leukocytes)),colour=Training),alpha=0.6)+
  scale_color_manual(values=c("steelblue","orange"),guide = guide_legend(reverse = TRUE))+scale_linetype_manual(values=c("twodash", "dotted"),guide = FALSE)+
  geom_line(data=CD3_events,aes(x=Time_point,y=((`CD16hiGRpos`*cbc_combo[,4]/100)/(Leukocytes)),colour=Training,group=Training,linetype=Training),size=0.7)+
  ylab("Normalized counts(thous/mcL)")+xlab("Time point")+ggtitle("CD16hi normalized counts across basline, peak and recovery")

##### CD193 Test 1:
CD193_events <-  read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=3)
CD193_events <- as.data.frame(unclass(CD193_events))
CD193_events <- cbind(CD193_events,(CD193_events$`CD193posGRpos`*cbc_193_203[,6])/(CD193_events$Leukocytes))
colnames(CD193_events)[11] <- "CD193_percentage"
cbc_193_203 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=5)
#levels(CD193_events$Training) <- c("Pre","Post")
quartz()
ggplot(CD193_events,aes(x=Time_point,y=CD193_percentage,group=desc(Training),fill=Training))+
  geom_bar(stat='summary',fun.y='mean',position=position_dodge(0.9),width = 0.55)+
  stat_summary(fun.data=mean_se,geom="errorbar",
               color="gray50",width=0.15,position = position_dodge(0.9))+
                 scale_fill_manual(values=c("steelblue","orange"),guide = guide_legend(reverse = TRUE))+
                 labs(x="Time point",y="normalized count for CD193(thousand/mcL)")

##### CD203
CD203_events <-  read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=4)
cbc_193_203 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=5)
CD203_events <- as.data.frame(unclass(CD203_events))
CD203_events <- cbind(CD203_events,(CD203_events$`CD203pGRp`*cbc_193_203[,5])/(CD203_events$Leukocytes))
colnames(CD203_events)[13] <- "CD203_percentage"
#levels(CD203_events$Training) <- c("Pre","Post")

quartz()
ggplot(CD203_events,aes(x=Time_point,y=CD203_percentage,group=desc(Training),fill=Training))+
  geom_bar(stat='summary',fun.y='mean',position=position_dodge(0.9),width = 0.55)+
  stat_summary(fun.data=mean_se,geom="errorbar",
               color="gray50",width=0.15,position = position_dodge(0.9))+
  scale_fill_manual(values=c("steelblue","orange"),guide = guide_legend(reverse = TRUE))+
  labs(x="Time point",y="normalized count for CD203 (thousand/mcL)")


##### CD3
CD3_event <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=1)
cbc_cd3 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=2)
CD3_events <- as.data.frame(unclass(CD3_event))
CD3_events <- cbind(CD3_events,(CD3_events$`CD3posGRpos`*cbc_cd3[,4])/(CD3_events$Leukocytes))
colnames(CD3_events)[14] <- "CD3_percentage"
ddply(CD3_events,~Training+Time_point,summarise,mean=mean(CD3_percentage),sd=sd(CD3_percentage))
quartz()
ggplot(CD3_events,aes(x=Time_point,y=CD3_percentage,group=desc(Training),fill=Training))+
  geom_bar(stat='summary',fun.y='mean',position=position_dodge(0.9),width = 0.55)+
  stat_summary(fun.data=mean_se,geom="errorbar",
               color="gray50",width=0.15,position = position_dodge(0.9))+
  scale_fill_manual(values=c("steelblue","orange"),guide = guide_legend(reverse = TRUE))+
  labs(x="Time point",y="normalized count for CD3 (thousand/mcL)")

##### CD14
cd14_events <- combo_events[,-c(10,11,13)]
cbc_combo <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=7)
CD14_events <- as.data.frame(unclass(cd14_events))
CD14_events <- cbind(CD14_events,(CD14_events$`CD14posGRpos`*cbc_combo[,5])/(CD14_events$Leukocytes))
colnames(CD14_events)[11] <- "CD14_percentage"

quartz()
ggplot(CD14_events,aes(x=Time_point,y=CD14_percentage,group=desc(Training),fill=Training))+
  geom_bar(stat='summary',fun.y='mean',position=position_dodge(0.9),width = 0.55)+
  stat_summary(fun.data=mean_se,geom="errorbar",
               color="gray50",width=0.15,position = position_dodge(0.9))+
  scale_fill_manual(values=c("steelblue","orange"),guide = guide_legend(reverse = TRUE))+
  labs(x="Time point",y="normalized count for CD14 (thousand/mcL)")

#Predicted values barplot
quartz()
ggplot(data=cbind(CD14_events,pred=fitted(cd14_mod4)*100),aes(x=Time_point,y=pred*100,group=desc(Training),fill=Training))+
  geom_bar(stat='summary',fun.y='mean',position=position_dodge(0.9),width = 0.55)+
  stat_summary(fun.data=mean_se,geom="errorbar",
               color="gray50",width=0.15,position = position_dodge(0.9))+
  scale_fill_manual(values=c("steelblue","orange"),guide = guide_legend(reverse = TRUE))+
  labs(x="Time point",y="normalized count for CD14 (thousand/mcL)")

quartz()
ggplot(CD14_events)+
  facet_wrap(~Subject)+
  geom_point(aes(x=Time_point,y=GRpos_CD193pos,colour=Training))+
  geom_line(data=cbind(mfi_cd193,pred=fitted(model7)),aes(x=Time_point,y=pred,colour=Training,group=Training))+
  ylab("MFI")+xlab("Time point")+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))

##### CD16high
cd16hi_events <- combo_events[,-c(11,12,13)]
cbc_combo <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/Events_all_populations.xlsx",sheet=7)
CD16hi_events <- as.data.frame(unclass(cd16hi_events))
CD16hi_events <- cbind(CD16hi_events,(CD16hi_events$`CD16hiGRpos`*cbc_combo[,4])/(CD16hi_events$Leukocytes))
colnames(CD16hi_events)[11] <- "CD16_percentage"
ddply(CD16hi_events,~Training+Time_point,summarise,mean=mean(CD16_percentage),sd=sd(CD16_percentage))
quartz()
ggplot(CD16hi_events,aes(x=Time_point,y=CD16_percentage,group=desc(Training),fill=Training))+
  geom_bar(stat='summary',fun.y='mean',position=position_dodge(0.9),width = 0.55)+
  stat_summary(fun.data=mean_se,geom="errorbar",
               color="gray50",width=0.15,position = position_dodge(0.9))+
  scale_fill_manual(values=c("steelblue","orange"),guide = guide_legend(reverse = TRUE))+
  labs(x="Time point",y="normalized count for CD16hi (thousand/mcL)")


####Test4: Healthy vs Asthmatic
CD193_events$Training <- factor(CD193_events$Training,levels = c("Pre","Post"))
quartz()
ggplot(CD193_events)+
  facet_wrap(~Training)+
  geom_boxplot(data=cbind(CD193_events,pred=fitted(cd193_mod3)),alpha=0.55,
               aes(x=Time_point,y=pred,fill=Condition,colour=Condition),width=0.35)+
xlab("Time point")+ylab("GR+ CD193+ normalized counts")+
  ggtitle("GR+ CD193+ Healthy vs Asthmatic normalized counts across the time points and training")+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  scale_colour_manual(values=c("#999999", "#E69F00"))

CD203_events$Training <- factor(CD203_events$Training,levels = c("Pre","Post"))
quartz()
ggplot(CD203_events)+
  facet_wrap(~Training)+
  geom_boxplot(data=cbind(CD203_events,pred=fitted(cd203_mod3)),alpha=0.55,
               aes(x=Time_point,y=pred,fill=Condition,colour=Condition),width=0.35)+
  xlab("Time point")+ylab("GR+ CD203+ normalized counts")+
  ggtitle("GR+ CD203+ Healthy vs Asthmatic normalized counts across the time points and training")+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  scale_colour_manual(values=c("#999999", "#E69F00"))

CD203_events$Training <- factor(CD203_events$Training,levels = c("Pre","Post"))
quartz()
ggplot(CD203_events)+facet_wrap(~Training)+
  geom_boxplot(aes(x=Time_point,y=((`CD203pGRp`*cbc_193_203[,5]/100)/(Leukocytes)),fill=Condition,colour=Condition),width=0.35,alpha=0.55)+
  xlab("Time point")+ylab("GR+ CD203+ normalized counts")+
  ggtitle("GR+ CD203+ Healthy vs Asthmatic normalized counts across the time points and training")+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  scale_colour_manual(values=c("#999999", "#E69F00"))


CD14_events$Training <- factor(CD14_events$Training,levels = c("Pre","Post"))
quartz()
ggplot(CD14_events)+
  facet_wrap(~Training)+
  geom_boxplot(data=cbind(CD14_events,pred=fitted(cd14_mod4)),alpha=0.55,
               aes(x=Time_point,y=pred,fill=Condition,colour=Condition),width=0.35)+
  xlab("Time point")+ylab("GR+ CD14+ normalized counts")+
  ggtitle("GR+ CD14+ Healthy vs Asthmatic normalized counts across the time points and training")+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  scale_colour_manual(values=c("#999999", "#E69F00"))

CD3_events$Training <- factor(CD3_events$Training,levels = c("Pre","Post"))
quartz()
ggplot(CD3_events)+
  geom_boxplot(data=cbind(CD3_events,pred=fitted(cd3_mod3)),alpha=0.55,
               aes(x=Time_point,y=pred,fill=Condition,colour=Condition),width=0.35)+
  facet_wrap(~Training)+xlab("Time point")+ylab("GR+ CD3+ normalized counts")+
  ggtitle("GR+ CD3+ Healthy vs Asthmatic normalized counts across the time points and training")+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  scale_colour_manual(values=c("#999999", "#E69F00"))

CD3_events$Training <- factor(CD3_events$Training,levels = c("Pre","Post"))
quartz()
ggplot(CD3_events)+
  geom_boxplot(data=cbind(CD3_events,pred=fitted(cd3_mod3)),alpha=0.55,
               aes(x=Time_point,y=pred,fill=Condition,colour=Condition),width=0.35)+
  facet_wrap(~Training)+xlab("Time point")+ylab("GR+ CD3+ normalized counts")+
  ggtitle("GR+ CD3+ Healthy vs Asthmatic normalized counts across the time points and training")+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  scale_colour_manual(values=c("#999999", "#E69F00"))

cd16hi_events$Training <- factor(cd16hi_events$Training,levels = c("Pre","Post"))
quartz()
ggplot(cd16hi_events)+
  geom_boxplot(data=cbind(cd16hi_events,pred=fitted(cd16hi_ev_mod3)),alpha=0.55,
               aes(x=Time_point,y=pred,fill=Condition,colour=Condition),width=0.35)+
  facet_wrap(~Training)+xlab("Time point")+ylab("GR+ CD16hi normalized counts")+
  ggtitle("GR+ CD16hi Healthy vs Asthmatic normalized counts across the time points and training")+
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  scale_colour_manual(values=c("#999999", "#E69F00"))


