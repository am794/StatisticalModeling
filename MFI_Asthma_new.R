#######UCI MFI######

#CD193
library("openxlsx")
files.list <- dir("/Volumes/Samsung_T5/UCI/Asthma/Oct11th/CD193_results/Gated/",pattern="\\.fcs")
MFI_CD193 <- list()
for (i in 1:length(files.list)) {
  MFI_table <- read.table(paste0("/Volumes/Samsung_T5/UCI/Asthma/Oct11th/CD193_results/Gated/",
                                   files.list[i],"/DUnSup_pop_MFI.txt"),header=F,row.names = 1)
  MFI_CD193[[i]] <- MFI_table[,3]
}
MFI_GR_CD193pos <- t(do.call(cbind,lapply(MFI_CD193,data.frame)))
rownames(MFI_GR_CD193pos) <- noquote(files.list)
colnames(MFI_GR_CD193pos) <- c("Leukocytes","GRpos","GRpos_CD193pos")
write.xlsx(MFI_GR_CD193pos,"/Volumes/Samsung_T5/UCI/Asthma/Oct11th/CD193_results/Gated/CD193_MFI.xlsx",quote=F,sep="\t",row.names=T)


#Combo
library("openxlsx")
files.list <- dir("/Volumes/Samsung_T5/UCI/Asthma/Oct11th/Combo_results/Gated/",pattern="\\.fcs")
MFI_Combo <- list()
for (i in 1:length(files.list)) {
  MFI_table <- read.table(paste0("/Volumes/Samsung_T5/UCI/Asthma/Oct11th/Combo_results/Gated/",
                                 files.list[i],"/DUnSup_pop_MFI.txt"),header=F,row.names = 1)
  MFI_Combo[[i]] <- MFI_table[,4]
}
MFI_GR_Combo <- t(do.call(cbind,lapply(MFI_Combo,data.frame)))
rownames(MFI_GR_Combo) <- noquote(files.list)
colnames(MFI_GR_Combo) <- c("Leukocytes","GRpos","CD16hiGR","CD16+GR+","CD14+GR+","CD14+GR_Backgated",
                            "Classical_Monocytes","Classical_Backgated","CD14+CD16+","Intermediate_Monocytes",
                            "Non-classical_Monocytes","Intermediate_Backgated","Non-classical_Backgated",
                            "CD3+GR+","CD3+GR+_Backgated")
write.xlsx(MFI_GR_Combo,"/Volumes/Samsung_T5/UCI/Asthma/Oct11th/Combo_results/Gated/Combo_MFI.xlsx",quote=F,sep="\t",row.names=T)

#CD203
library(openxlsx)
files.list <- dir("/Volumes/Samsung_T5/UCI/Asthma/Oct11th/CD203_results/Gated/",pattern="\\.fcs")
MFI_CD203 <- list()
for (i in 1:length(files.list)) {
  MFI_table <- read.table(paste0("/Volumes/Samsung_T5/UCI/Asthma/Oct11th/CD203_results/Gated/",
                                 files.list[i],"/DUnSup_pop_MFI.txt"),header=F,row.names = 1)
  MFI_CD203[[i]] <- MFI_table[,3]
}
MFI_GR_CD203 <- t(do.call(cbind,lapply(MFI_CD203,data.frame)))
rownames(MFI_GR_CD203) <- noquote(files.list)
colnames(MFI_GR_CD203) <- c("Leukocytes","GRpos","CD203_bisecting","CD203+_GR+","CD203+_GR+_Backgated")
write.xlsx(MFI_GR_CD203,"/Volumes/Samsung_T5/UCI/Asthma/Oct11th/CD203_results/Gated/CD203_MFI.xlsx",quote=F,sep="\t",row.names=T)

##############################################
# Normalizing the proportions w.r.t CBC data #
##############################################


