
CD4_data <- read.xlsx("/Users/amandava/Documents/Rochester_AsPIRES/AsPIRES_data/redo/AsPIRES_CD4_CD8_redo.xlsx",sheet = 4,colNames = TRUE,cols=c(1:14))
CD8_data <- read.xlsx("/Users/amandava/Documents/Rochester_AsPIRES/AsPIRES_data/redo/AsPIRES_CD4_CD8_redo.xlsx",sheet = 4,colNames = TRUE,cols=c(1:5,15))
comp_data <- read.xlsx("/Users/amandava/Documents/Rochester_AsPIRES/AsPIRES_data/redo/AsPIRES_CD4_CD8_redo.xlsx",sheet = 4,colNames = TRUE)

########### Interaction Plots ##########
Mild <- comp_data[which(comp_data$Condition == "Mild"),]
Severe <- comp_data[which(comp_data$Condition == "Severe"),]

i=15
quartz()
ggplot(data=Mild,aes(x=Visit,y=Mild[,i],group=SubjectID,color=SubjectID)) +
  geom_line() + geom_point(size=0.5) +
  facet_wrap( ~ Stimulation,scales = "free_y")+aes(color=SubjectID)+
  theme_bw()+ggtitle(paste0(colnames(Mild)[i]," for Mild "))+ylab("Population percentage")

quartz()
ggplot(data=Severe,aes(x=Visit,y=Severe[,i],group=SubjectID,color=factor(SubjectID))) +
  geom_line() + geom_point(size=0.5) +
  facet_wrap( ~ Stimulation,scales = "free_y")+aes(color=SubjectID)+
  theme_bw()+ggtitle(paste0(colnames(Severe)[i]," for Severe "))+ylab("Population percentage")

############# Boxplots ##############

# Test 5
t5_G5_12 <- comp_data %>% dplyr::filter(Stimulation %in% c("G5") & Visit %in% c("Visit11","Visit12"))
t5_M10_12 <- comp_data %>% dplyr::filter(Stimulation %in% c("M10") & Visit %in% c("Visit11","Visit12"))

quartz()
ggplot(data = t5_G5_12,aes(x=Visit,y=t5_G5_12[,i],fill=Visit))+geom_boxplot()+
  geom_line(aes(group=SubjectID,color=SubjectID))+
  geom_jitter(size=1.7,width=0.2)+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("G5 Popuation percentage for ",colnames(t5_G5_12)[i]," between Visit 11 and Visit 12"))

quartz()
ggplot(data = t5_M10_12,aes(x=Visit,y=t5_M10_12[,i],fill=Visit))+geom_boxplot()+
  geom_jitter(size=1.7)+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("M10 Popuation percentage for ",colnames(t5_M10_12)[i]," between Visit 11 and Visit 12"))

t5_G5_13 <- comp_data %>% dplyr::filter(Stimulation %in% c("G5") & Visit %in% c("Visit11","Visit13"))
t5_M10_13 <- comp_data %>% dplyr::filter(Stimulation %in% c("M10") & Visit %in% c("Visit11","Visit13"))
quartz()
ggplot(data = t5_G5_13,aes(x=Visit,y=t5_G5_13[,i],fill=Visit))+geom_boxplot()+
  geom_jitter(size=1.7)+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("G5 Popuation percentage for ",colnames(t5_G5_13)[i]," between Visit 11 and Visit 13"))
  
quartz()
ggplot(data = t5_M10_13,aes(x=Visit,y=t5_M10_13[,i],fill=Visit))+geom_boxplot()+
  geom_jitter(size=1.7)+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("M10 Popuation percentage for ",colnames(t5_M10_13)[i]," between Visit 11 and Visit 13"))


# Test 4
a1 <- comp_data %>% dplyr::filter(Visit %in% c("Visit11")  & Stimulation %in% c("DMSO","G5"))
t4_V11_G5 <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]
a1 <- comp_data %>% dplyr::filter(Visit %in% c("Visit11")  & Stimulation %in% c("DMSO","M10"))
t4_V11_M10 <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]

quartz()
ggplot(data = t4_V11_G5,aes(x=Stimulation,y=t4_V11_G5[,i],fill=Stimulation))+geom_boxplot(width=0.2)+
  geom_line(aes(group=SubjectID,color=SubjectID))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 11 Popuation percentage for ",colnames(t4_V11_G5)[i]," between DMSO and G5"))

quartz()
ggplot(data = t4_V11_M10,aes(x=Stimulation,y=t4_V11_M10[,i],fill=Stimulation))+geom_boxplot()+
  geom_line(aes(group=SubjectID,color=SubjectID))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 11 Popuation percentage for ",colnames(t4_V11_M10)[i]," between DMSO and M10"))

a1 <- comp_data %>% dplyr::filter(Visit %in% c("Visit12")  & Stimulation %in% c("DMSO","G5"))
t4_V12_G5 <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]
a1 <- comp_data %>% dplyr::filter(Visit %in% c("Visit12")  & Stimulation %in% c("DMSO","M10"))
t4_V12_M10 <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]

quartz()
ggplot(data = t4_V12_G5,aes(x=Stimulation,y=t4_V12_G5[,i],fill=Stimulation))+geom_boxplot()+
  geom_line(aes(group=SubjectID,color=SubjectID))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 12 Popuation percentage for ",colnames(t4_V12_G5)[i]," between DMSO and G5"))

quartz()
ggplot(data = t4_V12_M10,aes(x=Stimulation,y=t4_V12_M10[,i],fill= Stimulation))+geom_boxplot()+
  geom_line(aes(group=SubjectID,color=SubjectID))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 12 Popuation percentage for ",colnames(t4_V12_M10)[i]," between DMSO and M10"))

a1 <- comp_data %>% dplyr::filter(Visit %in% c("Visit13")  & Stimulation %in% c("DMSO","G5"))
t4_V13_G5 <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]
a1 <- comp_data %>% dplyr::filter(Visit %in% c("Visit13")  & Stimulation %in% c("DMSO","M10"))
t4_V13_M10 <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]

quartz()
ggplot(data = t4_V13_G5,aes(x=Stimulation,y=t4_V13_G5[,i],fill=Stimulation))+geom_boxplot()+
  geom_line(aes(group=SubjectID,color=SubjectID))+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 13 Popuation percentage for ",colnames(t4_V11_G5)[i]," between DMSO and G5"))

quartz()
ggplot(data = t4_V13_M10,aes(x=Stimulation,y=t4_V13_M10[,i],fill=Stimulation))+geom_boxplot()+
  geom_line(aes(group=SubjectID,color=SubjectID))+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 13 Popuation percentage for ",colnames(t4_V11_M10)[i]," between DMSO and M10"))

  
#Test 3
a1 <- comp_data %>% dplyr::filter(Stimulation %in% c("DMSO") & Condition %in% c("Mild") & Visit %in% c("Visit11","Visit12"))
m_DMSO_12 <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]
a1 <- comp_data %>% dplyr::filter(Stimulation %in% c("F1") & Condition %in% c("Mild") & Visit %in% c("Visit11","Visit12"))
m_F1_12 <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]

quartz()
ggplot(data = m_DMSO_12,aes(x=Visit,y=m_DMSO_12[,i],fill=Visit))+geom_boxplot()+
  geom_jitter(size=1.7)+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("DMSO Mild Popuation percentage for ",colnames(m_DMSO_12)[i]," between Visit 11 and Visit 12"))

quartz()
ggplot(data = m_F1_12,aes(x=Visit,y=m_F1_12[,i],fill=Visit))+geom_boxplot()+
  geom_jitter(size=1.7)+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("F1 Mild Popuation percentage for ",colnames(m_F1_12)[i]," between Visit 11 and Visit 12"))


a1 <- comp_data %>% dplyr::filter(Stimulation %in% c("DMSO") & Condition %in% c("Mild") & Visit %in% c("Visit11","Visit13"))
m_DMSO_13 <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]
a1 <- comp_data %>% dplyr::filter(Stimulation %in% c("F1") & Condition %in% c("Mild") & Visit %in% c("Visit11","Visit13"))
m_F1_13 <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]

quartz()
ggplot(data = m_DMSO_13,aes(x=Visit,y=m_DMSO_13[,i],fill=Visit))+geom_boxplot()+
  geom_jitter(size=1.7)+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("DMSO Mild Popuation percentage for ",colnames(m_DMSO_13)[i]," between Visit 11 and Visit 13"))

quartz()
ggplot(data = m_F1_13,aes(x=Visit,y=m_F1_13[,i],fill=Visit))+geom_boxplot()+
  geom_jitter(size=1.7)+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("F1 Mild Popuation percentage for ",colnames(m_F1_13)[i]," between Visit 11 and Visit 13"))

a1 <- comp_data %>% dplyr::filter(Stimulation %in% c("F1") & Condition %in% c("Severe") & Visit %in% c("Visit11","Visit13"))
s_F1_13 <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]

quartz()
ggplot(data = s_F1_13,aes(x=Visit,y=s_F1_13[,i],fill=Visit))+geom_boxplot()+
  geom_jitter(size=1.7)+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("F1 Severe Popuation percentage for ",colnames(s_F1_13)[i]," between Visit 11 and Visit 13"))

#Test 2
a1 <- comp_data %>% dplyr::filter(Visit %in% c("Visit11") & Condition %in% c("Mild")  & Stimulation %in% c("DMSO","F1"))
V11_DMSO_F1 <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]
a1 <- comp_data %>% dplyr::filter(Visit %in% c("Visit11") & Condition %in% c("Mild")  & Stimulation %in% c("DMSO","F5"))
V11_DMSO_F5 <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]
a1 <- comp_data %>% dplyr::filter(Visit %in% c("Visit11") & Condition %in% c("Mild")  & Stimulation %in% c("DMSO","NP1"))
V11_DMSO_NP1 <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]
a1_np5 <- comp_data %>% dplyr::filter(Visit %in% c("Visit11") & Condition %in% c("Mild")  & Stimulation %in% c("DMSO","NP5"))
V11_DMSO_NP1 <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]

quartz()
ggplot(data = a1_np5,aes(x=Stimulation,y=a1_np5[,i],fill=Stimulation))+geom_boxplot(width=0.35)+ 
  geom_line(aes(group=SubjectID,color=SubjectID))+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 11 Population percentage for ",colnames(V11_DMSO_F1)[i]," between DMSO and NP5 Mild"))

quartz()
ggplot(data = V11_DMSO_F1,aes(x=Stimulation,y=V11_DMSO_F1[,i],fill=Stimulation))+geom_boxplot(width=0.35)+ 
  geom_line(aes(group=SubjectID,color=SubjectID))+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 11 Population percentage for ",colnames(V11_DMSO_F1)[i]," between DMSO and F1 Mild"))

quartz()
ggplot(data = V11_DMSO_NP1,aes(x=Stimulation,y=V11_DMSO_NP1[,i],fill=Stimulation))+geom_boxplot()+
  geom_line(aes(group=SubjectID,color=SubjectID))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 11 Popuation percentage for ",colnames(V11_DMSO_NP1)[i]," between DMSO and NP1 Mild"))

a1 <- comp_data %>% dplyr::filter(Visit %in% c("Visit12") & Condition %in% c("Mild")  & Stimulation %in% c("DMSO","F1"))
V12_DMSO_F1 <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]
a1_np1_v12 <- comp_data %>% dplyr::filter(Visit %in% c("Visit12") & Condition %in% c("Mild")  & Stimulation %in% c("DMSO","NP1"))


quartz()
ggplot(data = a1_np1_v12,aes(x=Stimulation,y=a1_np1_v12[,i],fill=Stimulation))+geom_boxplot(width=0.35)+
  geom_line(aes(group=SubjectID,color=SubjectID))+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 12 Population percentage for ",colnames(a1_np1_v12)[i]," between DMSO and NP1 Milid"))

quartz()
ggplot(data = V12_DMSO_F1,aes(x=Stimulation,y=V12_DMSO_F1[,i],fill=Stimulation))+geom_boxplot(width=0.35)+
  geom_line(aes(group=SubjectID,color=SubjectID))+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 12 Population percentage for ",colnames(V12_DMSO_F1)[i]," between DMSO and F1 Milid"))

a1 <- comp_data %>% dplyr::filter(Visit %in% c("Visit13") & Condition %in% c("Mild")  & Stimulation %in% c("DMSO","F1"))
V13_DMSO_F1 <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]
a1_f5_v13 <- comp_data %>% dplyr::filter(Visit %in% c("Visit13") & Condition %in% c("Mild")  & Stimulation %in% c("DMSO","F5"))
a1_np1_v13 <- comp_data %>% dplyr::filter(Visit %in% c("Visit13") & Condition %in% c("Mild")  & Stimulation %in% c("DMSO","NP1"))

quartz()
ggplot(data = a1_f5_v13,aes(x=Stimulation,y=a1_f5_v13[,i],fill=Stimulation))+geom_boxplot(width=0.35)+
  geom_line(aes(group=SubjectID,color=SubjectID))+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 13 Population percentage for ",colnames(a1_f5_v13)[i]," between DMSO and F5 Mild"))

quartz()
ggplot(data = a1_np1_v13,aes(x=Stimulation,y=a1_np1_v13[,i],fill=Stimulation))+geom_boxplot(width=0.35)+
  geom_line(aes(group=SubjectID,color=SubjectID))+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 13 Population percentage for ",colnames(a1_np1_v13)[i]," between DMSO and NP1 Mild"))

quartz()
ggplot(data = V13_DMSO_F1,aes(x=Stimulation,y=V13_DMSO_F1[,i],fill=Stimulation))+geom_boxplot()+
  geom_line(aes(group=SubjectID,color=SubjectID))+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 13 Popuation percentage for ",colnames(V13_DMSO_F1)[i]," between DMSO and F1 Mild"))

a1 <- comp_data %>% dplyr::filter(Visit %in% c("Visit11") & Condition %in% c("Severe")  & Stimulation %in% c("DMSO","F1"))
V11_DMSO_F1_s <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]
a1 <- comp_data %>% dplyr::filter(Visit %in% c("Visit11") & Condition %in% c("Severe")  & Stimulation %in% c("DMSO","F5"))
V11_DMSO_F5_s <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]
a1 <- comp_data %>% dplyr::filter(Visit %in% c("Visit11") & Condition %in% c("Severe")  & Stimulation %in% c("DMSO","NP1"))
V11_DMSO_NP1_s <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]

quartz()
ggplot(data = V11_DMSO_F1_s,aes(x=Stimulation,y=V11_DMSO_F1_s[,i],fill=Stimulation))+geom_boxplot()+
  geom_line(aes(group=SubjectID,color=SubjectID))+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 11 Popuation percentage for ",colnames(V11_DMSO_F1_s)[i]," between DMSO and F1 Severe"))

quartz()
ggplot(data = V11_DMSO_F5_s,aes(x=Stimulation,y=V11_DMSO_F5_s[,i],fill=Stimulation))+geom_boxplot(width=0.35)+
  geom_line(aes(group=SubjectID,color=SubjectID))+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 11 Population percentage for ",colnames(V11_DMSO_F5_s)[i]," between DMSO and F5 Severe"))

quartz()
ggplot(data = V11_DMSO_NP1_s,aes(x=Stimulation,y=V11_DMSO_NP1_s[,i],fill=Stimulation))+geom_boxplot(width=0.35)+
  geom_line(aes(group=SubjectID,color=SubjectID))+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 11 Population percentage for ",colnames(V11_DMSO_NP1_s)[i]," between DMSO and NP1 Severe"))


a1 <- comp_data %>% dplyr::filter(Visit %in% c("Visit12") & Condition %in% c("Severe")  & Stimulation %in% c("DMSO","NP1"))
V12_DMSO_NP1_s <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]

quartz()
ggplot(data = V12_DMSO_NP1_s,aes(x=Stimulation,y=V12_DMSO_NP1_s[,i],fill=Stimulation))+geom_boxplot(width=0.35)+
  geom_line(aes(group=SubjectID,color=SubjectID))+geom_jitter(aes(color=SubjectID))+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 12 Population percentage for ",colnames(V12_DMSO_NP1_s)[i]," between DMSO and NP1 Severe"))



a1 <- comp_data %>% dplyr::filter(Visit %in% c("Visit13") & Condition %in% c("Severe")  & Stimulation %in% c("DMSO","NP1"))
V13_DMSO_NP1_s <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]
a1 <- comp_data %>% dplyr::filter(Visit %in% c("Visit13") & Condition %in% c("Severe")  & Stimulation %in% c("DMSO","NP5"))
V13_DMSO_NP5_s <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]
a1 <- comp_data %>% dplyr::filter(Visit %in% c("Visit13") & Condition %in% c("Severe")  & Stimulation %in% c("DMSO","F5"))
V13_DMSO_F5_s <- a1[a1$SubjectID %in% names(which(table(a1$SubjectID)==2)),]

quartz()
ggplot(data = V13_DMSO_NP1_s,aes(x=Stimulation,y=V13_DMSO_NP1_s[,i],fill=Stimulation))+geom_boxplot(width=0.35)+
  geom_line(aes(group=SubjectID,color=SubjectID))+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 13 Population percentage for ",colnames(V12_DMSO_NP1_s)[i]," between DMSO and NP1 Severe"))

quartz()
ggplot(data = V13_DMSO_NP5_s,aes(x=Stimulation,y=V13_DMSO_NP5_s[,i],fill=Stimulation))+geom_boxplot(width=0.35)+
  geom_line(aes(group=SubjectID,color=SubjectID))+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 13 Population percentage for ",colnames(V13_DMSO_NP5_s)[i]," between DMSO and NP5 Severe"))

quartz()
ggplot(data = V13_DMSO_F5_s,aes(x=Stimulation,y=V13_DMSO_F5_s[,i],fill=Stimulation))+geom_boxplot(width=0.35)+
  geom_line(aes(group=SubjectID,color=SubjectID))+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 13 Population percentage for ",colnames(V13_DMSO_F5_s)[i]," between DMSO and F5 Severe"))

### Test1
F1_12 <- comp_data %>% dplyr::filter(Visit %in% c("Visit12") & Stimulation %in% c("F1"))
DMSO_12 <- comp_data %>% dplyr::filter(Visit %in% c("Visit12") & Stimulation %in% c("DMSO"))

quartz()
ggplot(data = DMSO_12,aes(x=Condition,y=DMSO_12[,i],fill=Condition))+geom_boxplot()+
  geom_jitter(size=1.7)+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 12 DMSO Popuation percentage for ",colnames(DMSO_12)[i]," between Mild and Severe"))
quartz()
ggplot(data = F1_12,aes(x=Condition,y=F1_12[,i],fill=Condition))+geom_boxplot()+
  geom_jitter(size=1.7)+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_grey()+ylab("Population percentage")+
  ggtitle(paste0("Visit 12 F1 Popuation percentage for ",colnames(F1_12)[i]," between Mild and Severe"))
