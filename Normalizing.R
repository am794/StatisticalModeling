####################################################################
# Normalising the data across the panels based on the median values
####################################################################

#Calculate the median looping through the directories and files
setwd("/Users/amandava/Desktop/command_line_uci/")
dir_list <- c("Combo", "CD193", "CD203", "CD3_CD56")
for (i in dir_list){
files <-  list.files(path=paste0(i,"/Preprocessed_TXT"),pattern = "*.txt", full.names=T, recursive=FALSE)
median_values <- lapply(files,function(x){
t <- read.table(x,header=TRUE)
median_t <- median(t[,3])
})
median <- median(unlist(median_values))
assign(paste("median_",i,sep=""),median)
}
median_panels <- c(median_Combo, median_CD3_CD56, median_CD193, median_CD203)
ref_median <- min(median_panels)
diff <- as.data.frame(abs(median_panels-ref_median))
rownames(diff) <- dir_list

for(i in dir_list){
  dirs <- list.dirs(path=paste0(i,"/Gated_TXT"), full.names=T,recursive=FALSE)
  for(j in dirs){
    file_new <- list.files(path = )
  }
  file <- read.table()
}


ex <- read.table("/Users/amandava/Desktop/command_line_uci/CD3_CD56/Preprocessed_TXT/1012_gc_CD3apc_CD56pecy5_2014-09-20.001.txt",header=TRUE)
norm_ex <- ex
for(i in names(ex)){
  norm_ex[,i] <- ex[,i] - diff[[1]][4] 
}
#out <- function(
#write.table(out,"",sep="\t")