#########################
# GLMM for AsPIRES data #
#########################
library("lme4")
library("reshape2")
library("reshape")
library("ggplot2")
library("openxlsx")
library("emmeans")
library("dplyr")
library("plyr")
library("tidyr")
library("NMF")

genmix_dat <- read.xlsx("/Users/amandava/Documents/Rochester_AsPIRES/AsPIRES_data/redo/events_all.xlsx",sheet = 3,colNames = TRUE)
genmix_data <- genmix_dat[!(genmix_dat$Stimulation == "G10" | genmix_dat$Stimulation == "M20"),]
genmix_data$Stimulation <- as.factor(genmix_data$Stimulation)
genmix_data$Condition <- as.factor(genmix_data$Condition)
genmix_data$Visit <- as.factor(genmix_data$Visit)
genmix_data$SubjectID <- as.factor(genmix_data$SubjectID)
genmix_d <- genmix_data[which(genmix_data$Stimulation != "F1" & genmix_data$Stimulation != "NP1" & genmix_data$Stimulation != "SEB"),-3]
genmix_d$Stimulation <- factor(genmix_d$Stimulation,levels=c("DMSO","F5","G5","M10","NP5"))

Test_3 <- list(list())
for (i in 1:ncol(genmix_d)) {
  if(i>11 & i< 22){
    genmod_1.f <- genmix_d[,i]/genmix_d[,10] ~ Stimulation * Visit + (1|SubjectID)
    genfit_1.f <- lmer(genmod_1.f,weights=genmix_d[,10],data = genmix_d)
    isSingular(genfit_1.f) #FALSE
    ls_3 <- emmeans(genfit_1.f,~Stimulation| Visit)
    Test_3[[i-11]] <- contrast(ls_3,interaction = "trt.vs.ctrl",adjust="Dunnett")
    #quartz()
    #plot_emm <- emmip(genfit_1.f,~Stimulation|Visit,type="response",CI=T)+theme_gray()+ggtitle(paste0("Mean and 95% CI for ",colnames(genmix_d)[i]))
    #print(plot_emm)
    p <- genmix_d[,i]/genmix_d[,10]
    #quartz()
    #plot_box <- ggplot(data=genmix_d)+geom_boxplot(aes(x=Stimulation,y=p,color=Stimulation))+theme_gray()+
      #ylab("Cell proportion")+facet_wrap(~Visit)+ggtitle(paste0("Boxplot for ",colnames(genmix_d)[i]))
    #print(plot_box)
  }else if(i>21){
    genmod_1.f <- genmix_d[,i]/genmix_d[,11] ~ Stimulation * Visit + (1|SubjectID)
    genfit_1.f <- lmer(genmod_1.f,weights=genmix_d[,11],data = genmix_d)
    ls_3 <- emmeans(genfit_1.f,~Stimulation| Visit)
    Test_3[[i-11]] <- contrast(ls_3,interaction = "trt.vs.ctrl",adjust="Dunnett")
    #quartz()
    #plot_emm <- emmip(genfit_1.f,~Stimulation|Visit,type="response",CI=T)+theme_gray()+ggtitle(paste0("Mean and 95% CI for ",colnames(genmix_d)[i]))
    #print(plot_emm)
    p <- genmix_d[,i]/genmix_d[,10]
    #quartz()
    #plot_box <- ggplot(data=genmix_d)+geom_boxplot(aes(x=Stimulation,y=p,color=Stimulation))+theme_gray()+
      #ylab("Cell proportion")+facet_wrap(~Visit)+ggtitle(paste0("Boxplot for ",colnames(genmix_d)[i]))
    #print(plot_box)
  }
}

test3_pvals <- as.data.frame(Test_3)[,c(1,2,7,14,21,28,35,42,49,56,63,70,77,84)]
colnames(test3_pvals) <- c("Contrast", "Visit", colnames(genmix_d)[12:23])
#test3_pvals <- test3_pvals[,-8]
p_adj <- apply(test3_pvals[,3:14],1,p.adjust,method="BH")
#rownames(p_adj_fin) <- test3_pvals[,1]
p_adj_f <- cbind(as.data.frame(test3_pvals[,2]),t(p_adj))
colnames(p_adj_f)[1] <- "Visit"
p_adj_f <- cbind(as.data.frame(test3_pvals[,1]),p_adj_f)
colnames(p_adj_f)[1] <- "Contrast"
melted_p <- melt(data=p_adj_f,id.vars = c("Contrast","Visit"))
colnames(melted_p)[4] <- "pvalue"


p_adj_f[,3:14] <- -log10(p_adj_f[,3:14])
pp <- p_adj_f
colnames(pp)[3:14] <- c("CD4T/CD69+TNFa+","CD4T/CD69+IL2+","CD4T/CD69+IL4+",
                        "CD4T/CD69+IL10+","CD4T/CD69+IL17+","CD4_CD69_Ki67",
                        "CD4T/CD69+IFNg+","CD4T/CD25+FoxP3+_bi",
                        "CD4T/CD25+FoxP3+", "CD4T/CD25+FoxP3_IL10+","CD8T/CD69+TNFa+",     
                        "CD8T/IFNg+")
pp <- pp[order(pp$Contrast),]
pp_n <- pp[,c(3:9,11:14)]
#write.xlsx(pp,"/Users/amandava/Documents/Rochester_AsPIRES/AsPIRES_data/AsPIRES_neglog10_feb26.xlsx",rownames=T)
#rownames(pp_n) <- pp[,1]
#pp[,1] <- paste0(rep(c(1,2,3),each=4),pp[,1])
anno <- pp[,1:2]
annotation <- (anno %>% separate(Contrast,c("Contrast","DMSO")))[,c(1:3)]
t <- t(pp_n) > 1.30103
t[t=="TRUE"] <- "*"
t[t=="FALSE"] <- " "

quartz()
aheatmap(t(pp_n),annCol=list(Stimulation=annotation$Contrast,Visit=annotation$Visit),Colv=NA,
         Rowv= NA,cellwidth = 16, cellheight = 18, fontsize = 8,labCol = NA,color = 1L,
         breaks=seq(0,5,0.1),main = "-log10 BH corrected p-values",txt=t)

quartz()
aheatmap(t(pp_n),annCol=list(Stimulation=annotation$Contrast,Visit=annotation$Visit),Colv=NA,
         Rowv= NA,cellwidth = 16, cellheight = 18, fontsize = 8,labCol = NA,color = "-RdYlBu2",
         breaks=seq(0,5,0.1),main = "-log10 BH corrected p-values",txt=t)

