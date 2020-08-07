#################
#Plots
#################
library("gdata")
library("reshape2")
library("ggplot2")
library("plyr")
library("lattice")
setwd("/Users/amandava/Desktop/command_line_uci/uci_plots")
#mfi <- read.xls("/Users/amandava/Desktop/command_line_uci/uci_plots/MFI_new.xlsx",header=TRUE)
mfi <- read.xls("/Users/amandava/Desktop/command_line_uci/test_new/UCI_stats/MFI_stats.xlsx",header=TRUE,sheet=1)
row.names(mfi) <- mfi[1:23,1]
mfi <- mfi[,2:9]

####Fig2
my_mean <- apply(mfi[,1:7], 2, mean, na.rm=TRUE)
mfi_mean <- rbind(mfi[,1:7],my_mean)
se <- std.error(mfi[,1:7],na.rm = TRUE)
quartz()
bp = barplot(my_mean,names.arg = FALSE,ylim = c(0,3000),
ylab="GR Expression-MFI",main="GR expression by leukocyte subtype")#,names.arg = colnames(mfi)[1:7],beside = TRUE,cex.names = 0.75, ylab = "MFI")
text(x=bp, srt = 45,y = par("usr")[3] - 1,
     adj = 1, labels =colnames(mfi)[1:7] , xpd = TRUE,cex=1)
segments(bp, my_mean-se,bp,my_mean+se,lwd=1.5)
arrows(bp,my_mean-se,bp,my_mean+se,lwd=1.5,angle = 90,code=3,length=0.05)

####Fig3
####Fig3
mfi_sub <- mfi[,c(1,2,4,5)]
mfi_new<- data.frame(matrix(ncol = 4, nrow = 23))
mfi_new[,1] <- mfi_sub$T.cells*df2$T.cells
mfi_new[,2] <- mfi_sub$Monocytes*df2$Monocytes
mfi_new[,3] <- mfi_sub$NK.cells*df2$NK.cells
mfi_new[,4] <- mfi_sub$NKT.cells*df2$NKT.cells
mean_Samp<- data.frame(matrix(ncol = 4, nrow = 23))
df3 <- df2[,c(1,2,4,5)]
mean_Samp<-c()
for(i in 1:23){
  mean_Samp[i] <- weighted.mean(mfi_sub[i,],df3[i,],na.rm=TRUE)
}
std.error(mean_Samp)
mfi_pbmc <- as.data.frame(cbind(mfi_sub,mean_Samp))
mfi_pbmc <- as.data.frame(cbind(mfi_pbmc,mfi$Gender))
mfi_pbmc <- as.data.frame(cbind(mfi_pbmc,se))
quartz()
ggplot(data = mfi_pbmc,aes(x=mfi_pbmc$`mfi$Gender`, y=mean_Samp))+ 
  geom_dotplot(position="identity",binaxis = 'y',stackdir = 'center',dotsize = 0.3)+xlab("Gender")+ylab("MFI")+
  stat_summary(fun.data = mean_sdl,fun.args = list(mult=1),geom = "errorbar",
               color="blue",width=0.1)+stat_summary(fun.y=mean,geom = "point",color="blue")+ 
  scale_x_discrete(limits=c("Male","Female"))


#####Fig4
levels <- c("Male","Female")
melted$Gender <- levels
melted <- melt(mfi,na.rm = TRUE)
data_summ_2 <- ddply(melted, c("variable","Gender"), summarize,n=length(value),
                     sd = sd(value),
                     se = sd/sqrt(n),
                     mean=mean(value))
data_summ_2$Gender <- levels
quartz()
#ggplot(data_summ_2,aes(x=variable,y=mean,fill=Gender))+geom_bar(position = "dodge")
barchart(data_summ_2$mean~data_summ_2$variable,groups=data_summ_2$Gender,xlab = list(label = "Cell type", 
         fontsize = 20),ylab = list(label = "MFI", fontsize = 20), 
         auto.key = TRUE)

list(space = "inside", rectangles = TRUE, points = FALSE)), ll = data_summ$llim,
         ul = data_summ$ulim)
#, panel.segments(as.numeric(x), ll,as.numeric(x), ul,col = 'black', lwd = 1))


###Fig5
#mono_1 <- read.xls("/Users/amandava/Desktop/command_line_uci/uci_plots/monocyte_subtypes.xlsx",header=TRUE)
mono <- read.xls("/Users/amandava/Desktop/command_line_uci//test_new/UCI_stats/MFI_stats.xlsx",header=TRUE,sheet=2)
row.names(mono) <- mono[,1]
mono <- mono[,2:5]
melted_mono <- melt(mono,na.rm = TRUE)
data_summ_2 <- ddply(melted_mono, c("variable","Gender"), summarize,n=length(value),
                     sd = sd(value),
                     se = sd/sqrt(n),
                     mean=mean(value))
quartz()
barchart(data_summ_2$mean ~ data_summ_2$variable,groups=data_summ_2$Gender,xlab = list(label = "Cell type", 
                                                                                       fontsize = 20),ylab = list(label = "MFI", fontsize = 20),
         ylim=c(0,2900),
         auto.key = list(space = "inside", rectangles = TRUE, points = FALSE), ll = data_summ$llim,
         ul = data_summ$ulim)
#, panel.segments(as.numeric(x), ll,as.numeric(x), ul,col = 'black', lwd = 1))


#wilcoxon rank sum test
#significant difference in GR expression between male and female subjects
for(i in 1:7){
  wil_test <- wilcox.test(mfi[,1]~Gender, data=mfi, alternative="two.sided")
  p_vals[i] <- wil_test$p.value
}
names(p_vals) <- colnames(mfi)[1:7]
wil_t <- wilcox.test(T.cells ~ Gender ,data=mfi,alternative="two.sided")

#Proportions
props <- read.xls("/Users/amandava/Desktop/command_line_uci/uci_plots/proportions_all_panels.xlsx",header=TRUE,sheet=5)
props <- as.data.frame(t(props))
colnames(props) <- c("T.cells", "Monocytes", "Granulocytes", "NK.cells", "NKT.cells",
                     "Eosinophils", "CD203+ basophils", "Classical" ,"Intermediate", "Non-classical")
props <- as.data.frame(props[2:24,])
df2 <- data.frame(sapply(props, function(x) as.numeric(as.character(x))))
df3 <- df2[,c(1,2,4,5)]
mel <- melt(df3)
quartz()
ggplot(data = mel)+geom_boxplot(aes(x=variable,y=value,fill=variable))+
  ggtitle("cell population proportions")+theme_bw()+
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
  text = element_text(size = 12, family = "Tahoma"),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 11)) +scale_fill_brewer(palette = "Accent")  

Gender <- mfi$Gender
df4<-cbind(df3,Gender)
melt_2 <- melt(df4)
quartz()
ggplot(data = melt_2)+geom_boxplot(aes(x=variable,y=value,fill=Gender))+
  ggtitle("cell population proportions")+theme_bw()+
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 11)) +scale_fill_brewer(palette = "Accent")  


############################################################################################################



####Fig2
melted_2 <- melt(mfi)
melted <- na.omit(melted)
cdata <- ddply(melted, c("Gender","variable"), summarise,
               N    = length(melted$value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N))
quartz()
bp = barplot(cdata$mean, beside=TRUE, las=1, ylim=c(0,3000),
             cex.names=0.75,ylab="MFI",axes=TRUE,legend.text = TRUE,
             args.legend=list(title="Gender",legend=c("Female","Male"),x="topright",cex=0.7),col = c("red","blue"))
text(x=bp, srt = 45,y = par("usr")[3] - 1,
     adj = 1, labels =colnames(mfi)[1:7] , xpd = TRUE,cex=1)
segments(bp, cdata$mean-cdata$se*2,bp,cdata$mean+cdata$se*2,lwd=1.5)
arrows(bp,cdata$mean-cdata$se*2,bp,cdata$mean+cdata$se*2,lwd=1.5,angle = 90,code=3,length=0.05)


######2
melted <- melt(mfi,na.rm = TRUE)
data_summ <- ddply(melted, c("variable"), summarize,n=length(value),
               sd = sd(value),
               se = sd/sqrt(n),
               mean=mean(value),na.rm=TRUE)
data_summ$ulim <- data_summ$mean + data_summ$se
data_summ$llim <- data_summ$mean - data_summ$se
quartz()
barchart(data_summ$mean ~ data_summ$variable,xlab = list(label = "Cell type", 
                                                         fontsize = 20),
         ylab = list(label = "MFI", fontsize = 20), ll = data_summ$llim,
         ul = data_summ$ulim), panel.segments(as.numeric(x), ll,  
                                            as.numeric(x), ul,
                                            col = 'black', lwd = 1))

data_summ_2 <- ddply(melted, c("variable","Gender"), summarize,n=length(value),
                   sd = sd(value),
                   se = sd/sqrt(n),
                   mean=mean(value))
quartz()
barchart(data_summ_2$mean ~ data_summ_2$variable,groups=data_summ_2$Gender,xlab = list(label = "Cell type", 
        fontsize = 20),ylab = list(label = "MFI", fontsize = 20),
         auto.key = list(space = "inside", rectangles = TRUE, points = FALSE), ll = data_summ$llim,
         ul = data_summ$ulim), panel.segments(as.numeric(x), ll,  
                                             as.numeric(x), ul,
                                             col = 'black', lwd = 1))
##fig2


  
