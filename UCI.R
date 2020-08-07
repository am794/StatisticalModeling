source('~/Desktop/amandava/FCSTrans_FCS3.1_R3.0.r')
source("https://bioconductor.org/biocLite.R")
biocLite("flowViz")
require("flowViz")
require("ggplot2")
require("reshape2")
require("marray")
require("flowCore")
require("lattice")
read_combo <- read.FCS("/Users/amandava/Desktop/UCI_Data_29April2017/HealthyAdultsGR/Combo/1013_gc_CD3apc_CD16viob_CD14percp_2014-09-25.001.fcs",emptyValue = FALSE)
read_cd193 <- read.FCS("/Users/amandava/Desktop/UCI_Data_29April2017/HealthyAdultsGR/CD193/1013_gc_CD193apc_2014-09-25.001.fcs", emptyValue = FALSE)
read_cd203 <- read.FCS("/Users/amandava/Desktop/UCI_Data_29April2017/HealthyAdultsGR/CD203/1013_gc_CD203apc_2014-09-25.001.fcs", emptyValue = FALSE)
read_un <- read.FCS("/Users/amandava/Desktop/UCI_Data_29April2017/HealthyAdultsGR/UN/1013_UN_2014-09-25.001.fcs",emptyValue = FALSE)
fcs_combo <- as.data.frame(convertfcs(read_combo)[,2:13])
fcs_cd193 <- as.data.frame(convertfcs(read_cd193)[,2:13])
fcs_cd203 <- as.data.frame(convertfcs(read_cd203)[,2:13])
fcs_un <- as.data.frame(convertfcs(read_un)[,2:13])
names_fcs <- c("FSC-A","SSC-A","VioBlue-A","FITC-A","PE-A","PI/PE-Cy5.5-A","PE-Cy7-A","APC-A","APC-Cy-7-A")
rownames(fcs_combo[3:12]) <- names_fcs
write.table(fcs_combo,file = "/Users/amandava/Desktop/scatterplots/combo.txt",quote=F,sep = "\t", row.names = FALSE)
write.table(fcs_cd193,file = "/Users/amandava/Desktop/scatterplots/cd193.txt",quote=F,sep = "\t", row.names = FALSE)
write.table(fcs_cd203,file = "/Users/amandava/Desktop/scatterplots/cd203.txt",quote=F,sep = "\t", row.names = FALSE)
write.table(fcs_un,file = "/Users/amandava/Desktop/scatterplots/un.txt",quote=F,sep = "\t", row.names = FALSE)
png("Plot_combo.png", width=4.5, height = 4, units = 'in', res=300)
quartz()
splom(~fcs_combo[,4:12]|fcs_combo$`FL1-A`,data=fcs_combo,pch=10,cex=0.01,cex.labels=0.1,title="combo")
dev.off()

png("Plot_cd193.png", width=4.5, height = 4, units = 'in', res=300)
splom(~fcs_cd193[,4:12],pch=10,cex=0.01)
dev.off()

png("Plot_cd203.png", width=4.5, height = 4, units = 'in', res=300)
splom(~fcs_cd203[,4:12],pch=10,cex=0.01)
dev.off()
png("Plot_un.png", width=4.5, height = 4, units = 'in', res=300)
splom(~fcs_un[,4:12],pch=10,cex=0.01)
dev.off()


for(i in 1:2){
  for(j in 1:2)
{
  quartz()
  splom(~fcs_combo| fcs_combo$`FSC-A` ,pch=10,cex=0.01)
    }
}
#######################################################################################
?fcs_files <- c("fcs_combo","fcs_cd193","fcs_cd203","fcs_un")


########################################################################################
for(i in 4:12)
  quartz()
pairs(fcs_combo[,4:12],panel = points,pch=10,cex.labels = 0.1)

quartz()
pairs(fcs_cd193[,4:5],panel = points,pch=16,cex.labels = 0.001, labels = colnames(fcs_cd193))
xyplot(`SSC-A` ~ `FSC-A` | read_combo@exprs)
xyplot(fcs_cd203$`FSC-A`~fcs_cd203$`SSC-A`, smooth=FALSE)
ggplot(fcs_files[1], aes(x=fcs_files[1][,4], y = fcs_combo[,5]))+ geom_point(size=0.5)
ggplot(fcs_combo, aes(x=fcs_combo[,4],y=fcs_combo[,5]))+ geom_point(size=0.5)
quartz()
ggpairs(fcs_combo[,4:12])