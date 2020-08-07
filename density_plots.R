library(ggplot2)
library(reshape2)
combo_1013 <- read.table("/Users/amandava/Desktop/command_line_uci/Combo/693/694/1013_gc_CD3apc_CD16viob_CD14percp_2014-09-25.001/flock_results_all.txt",header=TRUE)
cd56_1013 <- read.table("/Users/amandava/Desktop/command_line_uci/CD3_CD56/708/709/1013_gc_CD3apc_CD56pecy5_2014-09-25.001/flock_results_all.txt",header=TRUE)
cd193_1013 <- read.table("/Users/amandava/Desktop/command_line_uci/CD193/705/706/1013_gc_CD193apc_2014-09-25.001/flock_results_all.txt", header=TRUE)
cd203_1013 <- read.table("/Users/amandava/Desktop/command_line_uci/CD203/702/703/1013_gc_CD203apc_2014-09-25.001/flock_results_all.txt", header=TRUE)
combo_1013_new <- combo_1013[combo_1013$pop1 == "0",]
cd56_1013_new <- cd56_1013[cd56_1013$pop1 == "0",]
cd193_1013_new <- cd193_1013[cd193_1013$pop1 == "0",]
cd203_1013_new <- cd203_1013[cd203_1013$pop1 == "0",]

quartz()
par(mfrow=c(2,2))

dens_combo <- density(as.numeric(combo_1013_new[,4]))
plot(dens_combo, main = "combo_1013", xlab = "GR-FITC")
polygon(dens_combo , col = "red", border = "gray")

dens_cd56 <- density(as.numeric(cd56_1013_new[,4]))
plot(dens_cd56, main = "cd56_1013", xlab = "GR-FITC")
polygon(dens_cd56 , col = "red", border = "gray")

dens_cd193 <- density(as.numeric(cd193_1013_new[,4]))
plot(dens_cd193, main = "cd193_1013", xlab = "GR-FITC")
polygon(dens_cd193 , col = "red", border = "gray")

dens_cd203 <- density(as.numeric(cd203_1013_new[,4]))
plot(dens_cd203, main = "cd203_1013" , xlab = "GR-FITC")
polygon(dens_cd203 , col = "red", border = "gray")

#frequency plots
quartz()
h<-hist(as.numeric(combo_1013_new[,4]),col="red",main = "", xlab="GR-FITC")
xfit <- seq(min(combo_1013_new[,4]), max(combo_1013_new[,4]))
yfit<-dnorm(xfit,mean=mean(as.numeric(combo_1013_new[,4])),sd=sd(combo_1013_new[,4])) 
yfit <- yfit*diff(h$mids[1:2])*length(combo_1013_new[,4])
lines(xfit, yfit, col="blue", lwd=2)

##########################################################################
cd56_list <- c("1012_gc_CD3apc_CD56pecy5_2014-09-20.001","1013_gc_CD3apc_CD56pecy5_2014-09-25.001",
               "1014_gc_CD3apc_CD56pecy5_2014-09-27.001","1015_GC+CD3apc+CD56peyc5_2014-09-29.001",
               "1016_GC+CD3apc+CD56peyc5_2014-09-29.002","1017_GC+CD3apc+CD56peyc5_2014-10-13.001",
               "1018_GC+CD3apc+CD56peyc5_2014-10-20.001","1019_GC+CD3apc+CD56peyc5_2014-10-20.001",
               "1020_GC+CD3apc+CD56peyc5_2014-11-21.001","1021_GC+CD3apc+CD56peyc5_2014-11-21.001",
               "1022_GC+CD3apc+CD56peyc5_2014-11-21.001","1023_GC+CD3apc+CD56peyc5_2015-02-03.001",
               "1024_GC+CD3apc+CD56peyc5_2015-02-06.001","1025_GC+CD3apc+CD56peyc5_2015-02-06.001",
               "1029_GC+CD3apc+CD56peyc5_2015-03-03.001","1030_GC+CD3apc+CD56peyc5_2015-03-17.001",
               "1031_GC+CD3apc+CD56peyc5_2015-04-14.001","1032_GC+CD3apc+CD56peyc5_2015-04-14.001",
               "1033_GC+CD3apc+CD56peyc5_2015-04-16.001","1035_GC+CD3apc+CD56peyc5_2015-06-09.001",
               "1036_GC+CD3apc+CD56peyc5_2015-06-16.001","1037_GC+CD3apc+CD56peyc5_2015-06-16.001",
               "1038_GC+CD3apc+CD56peyc5_2015-06-16.001")

quartz(width = 15, height = 12)
nf<-layout(matrix(1:24,4,6,byrow=TRUE), respect= TRUE)
par(oma=c(1,1,1,1),mar = c(2,2,2,2))
layout.show(24)
for(i in cd56_list){
  a<-paste0("/Users/amandava/Desktop/command_line_uci/CD3_CD56/745/746/",i)
  setwd(a)
  sample <- strsplit(i, "_")
  cd56 <- read.table("flock_results_all.txt",header=TRUE)
  cd56_new <- cd56[cd56$pop1 == "0",]
  dens_cd56 <- density(as.numeric(cd56_new[,4]))
  plot(dens_cd56,main = paste0("cd56_",sample[[1]][1]), xlab = "GR-FITC")
  polygon(dens_cd56 , col = "red", border = "black")
}
dev.off()

#frequency plots
quartz(width = 15, height = 12)
nf<-layout(matrix(1:24,4,6,byrow=TRUE), respect= TRUE)
par(oma=c(1,1,1,1),mar = c(2,2,2,2))
layout.show(24)
for(i in cd56_list){
  a<-paste0("/Users/amandava/Desktop/command_line_uci/CD3_CD56/745/746/",i)
  setwd(a)
  sample <- strsplit(i, "_")
  cd56 <- read.table("flock_results_all.txt",header=TRUE)
  cd56_new <- cd56[cd56$pop1 == "0",]
  h<-hist(cd56_new[,4],col="red", main = paste0("cd56_",sample[[1]][1]), xlab="GR-FITC")
  xfit <- seq(min(cd56_new[,4]), max(cd56_new[,4]))
  yfit<-dnorm(xfit,mean=mean(as.numeric(cd56_new[,4])),sd=sd(cd56_new[,4])) 
  yfit <- yfit*diff(h$mids[1:2])*length(cd56_new[,4])
  lines(xfit, yfit, col="blue", lwd=2)
}

##################################################################
combo_list <- c("1012_gc_CD3apc_CD16viob_CD14percp_2014-09-20.001","1013_gc_CD3apc_CD16viob_CD14percp_2014-09-25.001",
                "1014_gc_CD3apc_CD16viob_CD14percp_2014-09-27.001","1015_gcr_CD3apc_cd16viob_cd14perc_2014-09-29.001",
                "1016_gcr_CD3apc_cd16viob_cd14perc_2014-09-29.002","1017_gcr_CD3apc_cd16viob_cd14perc_2014-10-13.001",
                "1018_gcr_CD3apc_cd16viob_cd14perc_2014-10-20.001","1019_gcr_CD3apc_cd16viob_cd14perc_2014-10-20.001",
                "1020_gcr_CD3apc_cd16viob_cd14perc_2014-11-21.001","1021_gcr_CD3apc_cd16viob_cd14perc_2014-11-21.001",
                "1022_gcr_CD3apc_cd16viob_cd14perc_2014-11-21.001","1023_gcr_CD3apc_cd16viob_cd14perc_2015-02-03.001",
                "1024_gcr_CD3apc_cd16viob_cd14perc_2015-02-06.001","1025_gcr_CD3apc_cd16viob_cd14perc_2015-02-06.001",
                "1029_gcr_CD3apc_cd16viob_cd14perc_2015-03-03.001","1030_gcr_CD3apc_cd16viob_cd14perc_2015-03-17.001",
                "1031_gcr_CD3apc_cd16viob_cd14perc_2015-04-14.001","1032_gcr_CD3apc_cd16viob_cd14perc_2015-04-14.001",
                "1033_gcr_CD3apc_cd16viob_cd14perc_2015-04-16.001","1035_gcr_CD3apc_cd16viob_cd14perc_2015-06-09.001",
                "1036_gcr_CD3apc_cd16viob_cd14perc_2015-06-16.001","1037_gcr_CD3apc_cd16viob_cd14perc_2015-06-16.001",
                "1038_gcr_CD3apc_cd16viob_cd14perc_2015-06-16.001")

quartz(title="Combo_GR-FITC",width = 15, height = 12)
nf<-layout(matrix(1:24,4,6,byrow=TRUE), respect= TRUE)
par(oma=c(1,1,1,1),mar = c(2,2,2,2))
layout.show(24)
for(i in combo_list){
  a<-paste0("/Users/amandava/Desktop/command_line_uci/Combo/742/743/",i)
  setwd(a)
  sample <- strsplit(i, "_")
  combo <- read.table("flock_results_all.txt",header=TRUE)
  combo_new <- combo[combo$pop1 == "0",]
  dens_combo <- density(as.numeric(combo_new[,4]))
  plot(dens_combo,main = paste0("combo_",sample[[1]][1]), xlab = "GR-FITC")
  polygon(dens_combo , col = "red", border = "black")
}
dev.off()


#frequency plots
quartz(width = 15, height = 12)
nf<-layout(matrix(1:24,4,6,byrow=TRUE), respect= TRUE)
par(oma=c(1,1,1,1),mar = c(2,2,2,2))
layout.show(24)
for(i in combo_list){
  a<-paste0("/Users/amandava/Desktop/command_line_uci/Combo/742/743/",i)
  setwd(a)
  sample <- strsplit(i, "_")
  combo <- read.table("flock_results_all.txt",header=TRUE)
  combo_new <- cd56[cd56$pop1 == "0",]
  h<-hist(combo_new[,4],col="red", main = paste0("combo_",sample[[1]][1]), xlab="GR-FITC")
  xfit <- seq(min(combo_new[,4]), max(combo_new[,4]))
  yfit<-dnorm(xfit,mean=mean(as.numeric(combo_new[,4])),sd=sd(combo_new[,4])) 
  yfit <- yfit*diff(h$mids[1:2])*length(combo_new[,4])
  lines(xfit, yfit, col="blue", lwd=2)
}


###################################################################
cd203_list <- c("1012_gc_CD203apc_2014-09-20.001","1013_gc_CD203apc_2014-09-25.001",
"1014_gc_CD203apc_2014-09-27.001","1015_GC_CD203_apc_2014-09-29.001",
"1016_GC_CD203_apc_2014-09-29.001","1017_GC_CD203_apc_2014-10-13.001",
"1018_GC_CD203_apc_2014-10-20.001","1019_GC_CD203_apc_2014-10-20.001",
"1020_GC_CD203_apc_2014-11-21.001","1021_GC_CD203_apc_2014-11-21.001",
"1022_GC_CD203_apc_2014-11-21.001","1023_GC_CD203_apc_2015-02-03.001",
"1024_GC_CD203_apc_2015-02-06.001","1025_GC_CD203_apc_2015-02-06.001",
"1029_GC_CD203_apc_2015-03-03.001","1030_GC_CD203_apc_2015-03-17.001",
"1031_GC_CD203_apc_2015-04-14.001","1032_GC_CD203_apc_2015-04-14.001",
"1033_GC_CD203_apc_2015-04-16.001","1035_GC_CD203_apc_2015-06-09.001",
"1036_GC_CD203_apc_2015-06-16.001","1037_GC_CD203_apc_2015-06-16.001",
"1038_GC_CD203_apc_2015-06-16.001")

quartz(title="Cd203_GR-FITC",width = 15, height = 12)
nf<-layout(matrix(1:24,4,6,byrow=TRUE), respect= TRUE)
par(oma=c(1,1,1,1),mar = c(2,2,2,2))
layout.show(24)
for(i in cd203_list){
  a<-paste0("/Users/amandava/Desktop/command_line_uci/CD203/748/749/",i)
  setwd(a)
  sample <- strsplit(i, "_")
  cd203 <- read.table("flock_results_all.txt",header=TRUE)
  cd203_new <- cd203[cd203$pop1 == "0",]
  dens_cd203 <- density(as.numeric(cd203_new[,4]))
  plot(cd203_new[,4],main = paste0("cd203_",sample[[1]][1]), xlab = "GR-FITC")
  polygon(cd203_new[,4] , col = "red", border = "black")
}

#frequency plots
quartz(width = 15, height = 12)
nf<-layout(matrix(1:24,4,6,byrow=TRUE), respect= TRUE)
par(oma=c(1,1,1,1),mar = c(2,2,2,2))
layout.show(24)
for(i in cd203_list){
  a<-paste0("/Users/amandava/Desktop/command_line_uci/CD203/748/749/",i)
  setwd(a)
  sample <- strsplit(i, "_")
  cd203 <- read.table("flock_results_all.txt",header=TRUE)
  cd203_new <- cd203[cd203$pop1 == "0",]
  h<-hist(cd203_new[,4],col="red", main = paste0("cd203_",sample[[1]][1]), xlab="GR-FITC")
  xfit <- seq(min(cd203_new[,4]), max(cd203_new[,4]))
  yfit<-dnorm(xfit,mean=mean(as.numeric(cd203_new[,4])),sd=sd(cd203_new[,4])) 
  yfit <- yfit*diff(h$mids[1:2])*length(cd203_new[,4])
  lines(xfit, yfit, col="blue", lwd=2)
}

#########################################################################################
cd193_list <- c("1012_gc_CD193apc_2014-09-20.001","1013_gc_CD193apc_2014-09-25.001",
"1014_gc_CD193apc_2014-09-27.001","1015_GC+CD193_apc_2014-09-29.001",
"1016_GC+CD193_apc_2014-09-29.002","1017_GC+CD193_apc_2014-10-13.001",
"1018_GC+CD193_apc_2014-10-20.001","1019_GC+CD193_apc_2014-10-20.001",
"1020_GC+CD193_apc_2014-11-21.001","1021_GC+CD193_apc_2014-11-21.001",
"1022_GC+CD193_apc_2014-11-21.001","1023_GC+CD193_apc_2015-02-03.001",
"1024_GC+CD193_apc_2015-02-06.001","1025_GC+CD193APC_2015-02-06.001",
"1029_GC+CD193_apc_2015-03-03.001","1030_GC+CD193_apc_2015-03-17.001",
"1031_GC+CD193_apc_2015-04-14.001","1032_GC+CD193_apc_2015-04-14.001",
"1033_GC+CD193_apc_2015-04-16.001","1035_GC+CD193_apc_2015-06-09.001",
"1036_GC+CD193_apc_2015-06-16.001","1037_GC+CD193_apc_2015-06-16.001",
"1038_GC+CD193_apc_2015-06-16.001")

quartz(title="Cd193_GR-FITC",width = 15, height = 12)
nf<-layout(matrix(1:24,4,6,byrow=TRUE), respect= TRUE)
par(oma=c(1,1,1,1),mar = c(2,2,2,2))
layout.show(24)
for(i in cd193_list){
  a<-paste0("/Users/amandava/Desktop/command_line_uci/CD193/A/751/752/",i)
  setwd(a)
  sample <- strsplit(i, "_")
  cd193 <- read.table("flock_results_all.txt",header=TRUE)
  cd193_new <- cd193[cd193$pop1 == "0",]
  dens_cd193 <- density(as.numeric(cd193_new[,4]))
  plot(dens_cd193,main = paste0("cd193_",sample[[1]][1]), xlab = "GR-FITC")
  polygon(dens_cd193 , col = "red", border = "black")
}


