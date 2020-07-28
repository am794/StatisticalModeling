################################
# GLMM for flow cytometry data #
################################
library("lme4")
library("reshape2")
library("reshape")
library("ggplot2")
library("openxlsx")
library("emmeans")
library(dplyr)
library(NMF)
library(tidyr)
mixed_dat <- read.xlsx("/Users/amandava/Documents/Rochester_AsPIRES/AsPIRES_data/redo/AsPIRES_severity_socres.xlsx",sheet = 5,colNames = TRUE)
mixed_data <- mixed_dat[!(mixed_dat$Stimulation == "G10" | mixed_dat$Stimulation == "M20"),]
mixed_data$Stimulation <- as.factor(mixed_data$Stimulation)
mixed_data$Condition <- as.factor(mixed_data$Condition)
mixed_data$Visit <- as.factor(mixed_data$Visit)

genmix_dat <- read.xlsx("/Users/amandava/Documents/Rochester_AsPIRES/AsPIRES_data/redo/events_all.xlsx",sheet = 3,colNames = TRUE)
genmix_data <- genmix_dat[!(genmix_dat$Stimulation == "G10" | genmix_dat$Stimulation == "M20"),]
genmix_data$Stimulation <- as.factor(genmix_data$Stimulation)
genmix_data$Condition <- as.factor(genmix_data$Condition)
genmix_data$Visit <- as.factor(genmix_data$Visit)
genmix_data$SubjectID <- as.factor(genmix_data$SubjectID)
genmix_d <- genmix_data[which(genmix_data$Stimulation != "F1" & genmix_data$Stimulation != "NP1" & genmix_data$Stimulation != "SEB"),-3]
#genmix_d <- filter(genmix_data,Stimulation != "F1" & Stimulation != "NP1" & Stimulation != "SEB")
genmix_d$Stimulation <- factor(genmix_d$Stimulation,levels=c("DMSO","F5","G5","M10","NP5"))

##Fold change heatmap
fc_mean <- c()
for (i in 1:ncol(genmix_d)) {
  if(i>11 & i< 22){
    sub_data <- cbind(genmix_d[,1:5],genmix_d[,i]/genmix_d[,10])
    colnames(sub_data)[6] <- colnames(genmix_d)[i]
    wide <- dcast(sub_data,SubjectID~Visit+Stimulation)
    foldchange <- mutate(wide,DMSO_F5.Visit1=((wide[,3]/wide[,2])),DMSO_G5.Visit1=((wide[,4]/wide[,2])),
                         DMSO_M10.Visit1=((wide[,5]/wide[,3])),DMSO_NP5.Visit1=((wide[,6]/wide[,3])),
                         DMSO_F5.Visit2=((wide[,8]/wide[,7])),DMSO_G5.Visit2=((wide[,9]/wide[,7])),
                         DMSO_M10.Visit2=((wide[,10]/wide[,7])),DMSO_NP5.Visit2=((wide[,11]/wide[,7])),
                         DMSO_F5.Visit3=((wide[,13]/wide[,12])),DMSO_G5.Visit3=((wide[,14]/wide[,12])),
                         DMSO_M10.Visit3=((wide[,15]/wide[,12])),DMSO_NP5.Visit3=((wide[,16]/wide[,12])))[,c(1,17:28)]
    fc <- sapply(foldchange, function(x) replace(x, is.infinite(x),NA))
    foldchange_mean <- colMeans(fc[,-c(1)],na.rm = TRUE)
    fc_mean <- rbind(fc_mean,foldchange_mean)
  }else if(i>21){
    sub_data <- cbind(genmix_d[,1:5],genmix_d[,i]/genmix_d[,11])
    colnames(sub_data)[6] <- colnames(genmix_d)[i]
    wide <- dcast(sub_data,SubjectID~Visit+Stimulation)
    foldchange <- mutate(wide,DMSO_F5.Visit1=((wide[,3]/wide[,2])),DMSO_G5.Visit1=((wide[,4]/wide[,2])),
                         DMSO_M10.Visit1=((wide[,5]/wide[,3])),DMSO_NP5.Visit1=((wide[,6]/wide[,3])),
                         DMSO_F5.Visit2=((wide[,8]/wide[,7])),DMSO_G5.Visit2=((wide[,9]/wide[,7])),
                         DMSO_M10.Visit2=((wide[,10]/wide[,7])),DMSO_NP5.Visit2=((wide[,11]/wide[,7])),
                         DMSO_F5.Visit3=((wide[,13]/wide[,12])),DMSO_G5.Visit3=((wide[,14]/wide[,12])),
                         DMSO_M10.Visit3=((wide[,15]/wide[,12])),DMSO_NP5.Visit3=((wide[,16]/wide[,12])))[,c(1,17:28)]
    fc <- sapply(foldchange, function(x) replace(x, is.infinite(x),NA))
    foldchange_mean <- colMeans(fc[,-c(1)],na.rm = TRUE)
    fc_mean <- rbind(fc_mean,foldchange_mean)
  }
}
rownames(fc_mean) <- c("CD4T/CD69+TNFa+","CD4T/CD69+IL2+","CD4T/CD69+IL4+",
                       "CD4T/CD69+IL10+","CD4T/CD69+IL17+","CD4/CD69+Ki67+","CD4T/CD69+IFNg+","CD4T/CD25+FoxP3+_bi",
                       "CD4T/CD25+FoxP3+", "CD4T/CD25+FoxP3_IL10+","CD8T/CD69+TNFa+",
                       "CD8T/IFNg+")
mlt <- melt(fc_mean[c(2,3,5,7,10,12),])
colnames(mlt) <- c("Var1","Var2","value")
mlt_1 <- separate(data=mlt,col=Var2,into=c("Stimulation","Visit"),".Visit")
mlt_1 <- dcast(mlt_1,Stimulation+Visit~Var1)[,c(1,2,6,7,5,4,3,8)]
mlt_1 <- mlt_1[order(mlt_1$Stimulation),]
mlt1_n <- mlt_1[,c(3:8)]


######## ANNO ##########
Test_3 <- list(list())
for (i in 1:ncol(genmix_d)) {
  if(i>11 & i< 22){
    genmod_1.f <- genmix_d[,i]/genmix_d[,10] ~ Stimulation * Visit + (1|SubjectID)+(1|Sample)
    genfit_1.f <- glmer(genmod_1.f,weights=genmix_d[,10],family="binomial",data = genmix_d)
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
    genmod_1.f <- genmix_d[,i]/genmix_d[,11] ~ Stimulation * Visit + (1|SubjectID)+(1|Sample)
    genfit_1.f <- glmer(genmod_1.f,weights=genmix_d[,11],family="binomial",data = genmix_d)
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
test3_pvals <- test3_pvals[,-8]
p_adj <- apply(test3_pvals[,3:13],1,p.adjust,method="BH")
#rownames(p_adj_fin) <- test3_pvals[,1]
p_adj_f <- cbind(as.data.frame(test3_pvals[,2]),t(p_adj))
colnames(p_adj_f)[1] <- "Visit"
p_adj_f <- cbind(as.data.frame(test3_pvals[,1]),p_adj_f)
colnames(p_adj_f)[1] <- "Contrast"
melted_p <- melt(data=p_adj_f,id.vars = c("Contrast","Visit"))
colnames(melted_p)[4] <- "pvalue"

p_adj_f[,3:13] <- -log10(p_adj_f[,3:13])
pp <- p_adj_f
colnames(pp)[3:13] <- c("CD4T/CD69+TNFa+","CD4T/CD69+IL2+","CD4T/CD69+IL4+",
                        "CD4T/CD69+IL10+","CD4T/CD69+IL17+","CD4T/CD69+IFNg+","CD4T/CD25+FoxP3+_bi",
                        "CD4T/CD25+FoxP3+", "CD4T/CD25+FoxP3_IL10+","CD8T/CD69+TNFa+",
                        "CD8T/IFNg+")
pp <- pp[order(pp$Contrast),]
#pp_n <- pp[,c(3:8,10,11,13)]
pp_n <- pp[,c(4,5,7,8,11,13)]
#write.xlsx(pp,"/Users/amandava/Documents/Rochester_AsPIRES/AsPIRES_data/AsPIRES_neglog10_feb26.xlsx",rownames=T)
#rownames(pp_n) <- pp[,1]
#pp[,1] <- paste0(rep(c(1,2,3),each=4),pp[,1])
anno <- pp[,1:2]
annotation <- (anno %>% separate(Contrast,c("Contrast","DMSO")))[,c(1:3)]
annotation$Contrast <- c(rep("F",3),rep("G",3),rep("M",3),rep("NP",3))
t <- t(pp_n) > 1.30103
t[t=="TRUE"] <- "*"
t[t=="FALSE"] <- " "

############ HeatMap ###########
col_fun = colorRamp2(c(0,1,8), c("blue", "white", "red"))
col_fun(seq(-3, 3))
quartz()
aheatmap(t(mlt1_n),annCol=list(Stimulation=annotation$Contrast,Visit=annotation$Visit),Colv=NA,
         Rowv= NA,cellwidth = 16, cellheight = 18, fontsize = 8,labCol = NA,
         color = '-RdYlBu:100',
         breaks=NA,main = "Mean foldchange in comparison to DMSO",txt=t)
text(x=0.95,y=0.3,"* Significantly different",cex=0.6)


############ New heatmap ##############
ha <- HeatmapAnnotation(df=data.frame(Stimulation=annotation$Contrast,Visit=annotation$Visit))
Heatmap(t(mlt1_n),cluster_columns = FALSE, cluster_rows = FALSE,col = colorRamp2(c(0, 1,2, 8), c("darkblue","green", "yellow","red")),top_annotation = ha,show_column_names = FALSE)

quartz()
hp <- Heatmap(t(mlt1_n),cluster_columns = FALSE, name='',cluster_rows = F,
              col = colorRamp2(c(0, 1, 2,8), c("blue2","cyan","gold","firebrick3")),
              heatmap_width = unit(12, "cm"),heatmap_height = unit(5, "cm"),
              top_annotation = ha,show_column_names = FALSE,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
  grid.text(t[i, j], x, y)
})
draw(hp,annotation_legend_list=Legend(pch = "*", type = "points", labels = "Significantly different",size=0.1))

quartz()
hp <- Heatmap(t(mlt1_n),cluster_columns = FALSE, name='',cluster_rows = F,
              col = colorRamp2(c(0, 1,8), c("blue2","gold","firebrick3")),
              heatmap_width = unit(12, "cm"),heatmap_height = unit(5, "cm"),
              top_annotation = ha,show_column_names = FALSE,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(t[i, j], x, y)
              })
draw(hp,annotation_legend_list=Legend(pch = "*", type = "points", labels = "Significantly different",size=0.1))


quartz()
hp <- Heatmap(t(mlt1_n),cluster_columns = FALSE, name='',cluster_rows = F,
              col = colorRamp2(c(0,8), c("blue2","firebrick3")),
              heatmap_width = unit(12, "cm"),heatmap_height = unit(5, "cm"),
              top_annotation = ha,show_column_names = FALSE,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(t[i, j], x, y)
              })
draw(hp,annotation_legend_list=Legend(pch = "*", type = "points", labels = "Significantly different",size=0.1))
col_function = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))
ha_new <- HeatmapAnnotation(df=data.frame(Visit=annotation$Visit),col = list(Visit = c("Visit 1" =  "midnightblue", "Visit 2" = "brown4", "Visit 3" = "darkgreen")))
#Heatmap(t(mlt1_n),cluster_columns = FALSE, cluster_rows = FALSE,col = colorRamp2(c(0, 1,2, 8), c("darkblue","green", "yellow","red")),top_annotation = ha_new,show_column_names = FALSE)

quartz()
hp <- Heatmap(t(mlt1_n),cluster_columns = FALSE, name='',cluster_rows = F,
              col = colorRamp2(c(0, 1, 2,8), c("blue2","cyan","gold","firebrick3")),
              heatmap_width = unit(12, "cm"),heatmap_height = unit(5, "cm"),
              top_annotation = ha_new,show_column_names = FALSE,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(t[i, j], x, y)
              })
draw(hp,annotation_legend_list=Legend(pch = "*", type = "points", labels = "Significantly different",size=0.1))
