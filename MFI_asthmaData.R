############UCI MFI###############
#Combo
setwd("/Users/amandava/Desktop/UCI_Config_Reports/Combo/Combo_old_correct/Gated_Combo/")
files.list_gated <- dir("/Users/amandava/Desktop/UCI_Config_Reports/Combo/Combo_old_correct/Gated_Combo/",pattern = "\\.fcs")
GR_pos <- list(list())
for(i in 1:142){
  mfi <- read.table(paste0("/Users/amandava/Desktop/UCI_Config_Reports/Combo/Combo_old_correct/Gated_Combo/",files.list_gated[i],"/DUnSup_pop_MFI.txt"),
                       header=F,sep = "\t",row.names=1)
  GR_pos[[i]] <- mfi[,4]
}
GRpositive <- t(do.call(cbind, lapply(GR_pos, data.frame)))
row.names(GRpositive) <- files.list_gated
colnames(GRpositive) <- c("Leukocytes","GR+","CD16hiGR","CD16+GR+","CD14+GR+",
                          "CD14+GR_Backgated","Classical_Monocytes","Classical_Backgated",
                          "CD14+CD16+","Intermediate_Monocytes","Non-classicalMonocytes",
                          "Intermediate_Backgated","Non-classical_Backgated","CD3+GR+")
write.xlsx(GRpositive,"/Users/amandava/Desktop/UCI_Config_Reports/Combo/Combo_old_correct/MFI_GRpositive_combo.xlsx",quote=F,sep="\t",row.names=T)

#CD193
setwd("/Users/amandava/Desktop/UCI_Config_Reports/CD193/CD193_results/Gated_CD193/")
files.list_gated <- dir("/Users/amandava/Desktop/UCI_Config_Reports/CD193/CD193_results/Gated_CD193/",pattern = "\\.fcs")
GR_pos_cd193 <- list(list())
for(i in 1:145){
  mfi <- read.table(paste0("/Users/amandava/Desktop/UCI_Config_Reports/CD193/CD193_results/Gated_CD193/",files.list_gated[i],"/DUnSup_pop_MFI.txt"),
                    header=F,sep = "\t",row.names=1)
  GR_pos_cd193[[i]] <- mfi[,3]
}
GRpositive_cd193 <- t(do.call(cbind, lapply(GR_pos_cd193, data.frame)))
row.names(GRpositive_cd193) <- files.list_gated
colnames(GRpositive_cd193) <- c("Leukocytes","GR+","GR+CD193+","GR+CD193+_Backgated")
write.xlsx(GRpositive_cd193,"/Users/amandava/Desktop/UCI_Config_Reports/CD193/CD193_results/MFI_GRpositive_CD193.xlsx",quote=F,sep="\t",row.names=T)

#CD203
setwd("/Users/amandava/Desktop/UCI_Config_Reports/CD203/CD203_results/Gated_CD203/")
files.list_gated <- dir("/Users/amandava/Desktop/UCI_Config_Reports/CD203/CD203_results/Gated_CD203/",pattern = "\\.fcs")
GR_pos_cd203 <- list(list())
for(i in 1:141){
  mfi <- read.table(paste0("/Users/amandava/Desktop/UCI_Config_Reports/CD203/CD203_results/Gated_CD203/",files.list_gated[i],"/DUnSup_pop_MFI.txt"),
                    header=F,sep = "\t",row.names=1)
  GR_pos_cd203[[i]] <- mfi[,3]
}
GRpositive_cd203 <- t(do.call(cbind, lapply(GR_pos_cd203, data.frame)))
row.names(GRpositive_cd203) <- files.list_gated
colnames(GRpositive_cd203) <- c("Leukocytes","GR+","GR+CD203+_bisecting","CD203+_GR+_backgated",
                                "GR+CD203+_Clustering","GR+CD203+_clust_backgating")
write.xlsx(GRpositive_cd203,"/Users/amandava/Desktop/UCI_Config_Reports/CD203/CD203_results//MFI_GRpositive_CD203.xlsx",quote=F,sep="\t",row.names=T)


CD193 <- read.xlsx("/Users/amandava/Desktop/UCI_Config_Reports/MFI_all_panels_10Oct.xlsx",sheet=1)
filt_cd193 <- CD193[,c(2,4,5,6)]
CD193_pre <- filt_cd193[which(filt_cd193$Training == "Pre"),c(1,2,4)]
CD193_post <- filt_cd193[which(filt_cd193$Training == "Post"),c(1,2,4)]
quartz()
ggplot(data=CD193_pre,aes(x=Time.Point,y=`pop4_GR+CD193+_GR`,group=Subject,color=Subject))+
  geom_line(aes(linetype=Subject),size=0.8)+geom_point(aes(shape=Subject),size=1.5)+
  ggtitle("Median MFI of GR+CD193+ population for Pre training")+
  xlab("Time Point")+ylab("Median MFI of GR+CD193+ population")+
  scale_linetype_manual(values = c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5))+
  scale_shape_manual(values = c(0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,0))+
  ylim(c(1000,4096))

quartz()
ggplot(data=CD193_post,aes(x=Time.Point,y=`pop4_GR+CD193+_GR`,group=Subject,color=Subject))+
  geom_line(aes(linetype=Subject),size=0.8)+geom_point(aes(shape=Subject),size=1.5)+
  ggtitle("Median MFI of GR+CD193+ population for Post training")+
  xlab("Time Point")+ylab("Median MFI of GR+CD193+ population")+
  scale_linetype_manual(values = c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5))+
  scale_shape_manual(values = c(0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,0))+
  ylim(c(1000,4096))


