################################
# GLMM for flow Cytometry data #
################################
library("lme4")
library("reshape2")
library("reshape")
library("ggplot2")
library("openxlsx")
mixed_dat <- read.xlsx("/Users/amandava/Documents/Rochester_AsPIRES/AsPIRES_data/redo/AsPIRES_severity_socres.xlsx",sheet = 5,colNames = TRUE)
mixed_data <- mixed_dat[!(mixed_dat$Stimulation == "G10" | mixed_dat$Stimulation == "M20"),]
mixed_data$Stimulation <- as.factor(mixed_data$Stimulation)
mixed_data$Condition <- as.factor(mixed_data$Condition)
mixed_data$Visit <- as.factor(mixed_data$Visit)

##############Linear Models with interaction################
mod_1 <- CD4_CD69_IL2 ~ Condition + Stimulation + Visit + (1|SubjectID)
fit_1 <- lmer(mod_1,data=mixed_data,REML=FALSE)
summary(fit_1)

mod_2 <- CD4_CD69_IL2 ~ Condition + Stimulation + Visit + (1|SubjectID) + Stimulation*Visit
fit_2 <- lmer(mod_2,data=mixed_data,REML=FALSE)
summary(fit_1)

mod_1.a <- CD4_CD69_IL2 ~ Condition + Stimulation + Visit + (1|SubjectID) + Condition*Stimulation + Condition*Visit + Stimulation*Visit
fit_1.a <- lmer(mod_1.a,data=mixed_data,REML=FALSE)
summary(fit_1.a)

mod_1.b <- CD4_CD69_IL2 ~ Condition * Stimulation * Visit +  (1+Stimulation|SubjectID)
fit_1.b <- lmer(mod_1.b,data=mixed_data)
summary(fit_1.b)

mod_1.c <- CD4_CD69_IL2 ~ Severity_score * Stimulation * Visit +  (1|SubjectID)
fit_1.c <- lmer(mod_1.c,data=mixed_data)
summary(fit_1.c)

quartz();
plot(residuals(r))
plot(fit_1.c) #Residuals vs fitted
qqnorm(residuals(fit_1.c),ylim=c(-2,2)) #Q-Q plot
qqline(residuals(fit_1.c))
plot(residuals(fit_1.a)) #Residuals
plot(fit_1.a,which=c(1:3,5))
spreadLevelPlot(CD4_CD69_IL2 ~ Condition * Stimulation * Visit,data=mixed_data)
shapiro.test(residuals(fit_1.b))

mod_1_v2 <- CD4T_CD69_IL2 ~ Condition + Stimulation + Visit
fit_1_v2 <- lme(CD4_CD69_TNFa ~ Condition + Stimulation + Visit,random = ~1|SubjectID,data=mixed_data)
summary(fit_1_v2)
mod_1.a_v2 <- CD4_CD69_TNFa ~ Condition + Stimulation + Visit + Condition:Stimulation:Visit
fit_1.a_v2 <- lme(CD4_CD69_TNFa ~ Condition * Stimulation * Visit,random = ~1|SubjectID,data=mixed_data,weights=varPower())
summary(fit_1.a_v2)

############ GLMM ###################

#Count data
genmix_dat <- read.xlsx("/Users/amandava/Documents/Rochester_AsPIRES/AsPIRES_data/redo/events_all.xlsx",sheet = 3,colNames = TRUE)
genmix_data <- genmix_dat[!(genmix_dat$Stimulation == "G10" | genmix_dat$Stimulation == "M20"),]
genmix_data$Stimulation <- as.factor(genmix_data$Stimulation)
genmix_data$Condition <- as.factor(genmix_data$Condition)
genmix_data$Visit <- as.factor(genmix_data$Visit)
genmix_data$SubjectID <- as.factor(genmix_data$SubjectID)

genmod_1.a <- genmix_data[,12]/genmix_data[,10] ~ Condition + Stimulation + Visit +  (1|SubjectID)
genfit_1.a <- glmer(genmod_1.a,weights=genmix_data[,10],family=binomial,data = genmix_data)

genmod_1.b <- genmix_data[,12]/genmix_data[,10] ~ Condition + Stimulation + Visit +
             Condition * Stimulation + Stimulation * Visit + Condition * Visit + (1|SubjectID) +   (1+Stimulation | SubjectID)
genfit_1.b <- glmer(genmod_1.b,weights=genmix_data[,10],family=binomial,data = genmix_data)

genmod_1.c <- genmix_data[,12]/genmix_data[,10] ~ Condition * Stimulation * Visit +  (1|SubjectID)
genfit_1.c <- glmer(genmod_1.c,weights=genmix_data[,10],family=binomial,data = genmix_data)

genmod_1.d <- genmix_data[,12]/genmix_data[,10] ~ Condition * Stimulation * Visit +  (1|SubjectID) + (1+Stimulation | SubjectID)
genfit_1.d <- glmer(genmod_1.d,weights=genmix_data[,10],family=binomial,data = genmix_data)

##model e: each sample receives a unique level of a random effect that models the extra varaition present in the data
genmod_1.e <- genmix_data[,17]/genmix_data[,10] ~ Condition * Stimulation * Visit +  (1|SubjectID) + (1|Sample)
genfit_1.e <- glmer(genmod_1.e,weights=genmix_data[,10],family=binomial,data = genmix_data)




############ overdispersion function #############
overdisp_fun <- function(model) {
  ## number of variance parameters in
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
## do residual deviance/df from the output
#####################################################

ls <- emmeans::lsmeans(genfit_1.e,~Stimulation| Condition*Visit,adjust="Dunnett")
contrast(ls,interaction = "trt.vs.ctrl")

ls <- emmeans::lsmeans(genfit_1.e,~Visit|Condition* Stimulation,adjust="none")
contrast(ls,interaction = "pairwise")

ls <- emmeans::lsmeans(genfit_1.e,~Condition | Stimulation* Visit,adjust=none)
contrast(ls,interaction = "pairwise")

###################### Hash and keys #######################
key <- c(12,13,14,15,16,17,18,20,21,23)
value <- c(10,10,10,10,10,10,10,10,10,11)
h <- hash()
h <- hash(key,value)

Test_1 <- list(list())
Test_2 <- list(list())
Test_3 <- list(list())

genmod_1.a <- genmix_data[,13]/genmix_data[,11] ~ Condition + Stimulation + Visit +  (1|SubjectID) + (1|Sample)
genfit_1.a <- glmer(genmod_1.a,weights=genmix_data[,11],family=binomial,data = genmix_data)

genmod_1.b <- genmix_data[,16]/genmix_data[,11] ~ Condition + Stimulation + Visit +  (1|SubjectID) + (1|Sample) + Stimulation*Visit
genfit_1.b <- glmer(genmod_1.b,weights=genmix_data[,11],family=binomial,data = genmix_data)

genmod_1.e <- genmix_data[,13]/genmix_data[,11] ~ Condition * Stimulation * Visit +  (1|SubjectID) + (1|Sample)
genfit_1.e <- glmer(genmod_1.e,weights=genmix_data[,11],family=binomial,data = genmix_data)

ls_3 <- emmeans::lsmeans(genfit_1.e,~Stimulation|Condition * Visit,adjust="mvt")
Test_3[[10]] <- contrast(ls_3,interaction = "trt.vs.ctrl")

ls_2 <- emmeans::lsmeans(genfit_1.e,~Visit|Condition* Stimulation,adjust="Tukey")
Test_2[[10]] <- contrast(ls_2,interaction = "pairwise")

ls_1 <- emmeans::lsmeans(genfit_1.e,~Condition | Stimulation* Visit,adjust="Tukey")
Test_1[[10]] <- contrast(ls_1,interaction = "pairwise")

library("openxlsx")
library("plyr")
library("dplyr")
library("cowplot")
test1_pvals <- read.xlsx("/Users/amandava/Documents/Rochester_AsPIRES/AsPIRES_data/redo/Test1.xlsx")
p_adj_1 <- apply(test1_pvals[,4:13],1,p.adjust,method="BH")
t_1 <- cbind(test1_pvals[,1:3],t(p_adj_1))

test2_pvals <- read.xlsx("/Users/amandava/Documents/Rochester_AsPIRES/AsPIRES_data/redo/Test2.xlsx",rowNames = TRUE)
p_adj_2 <- apply(test2_pvals[,4:13],1,p.adjust,method="BH")
t_2 <- cbind(test2_pvals[,1:3],t(p_adj_2))

test3_pvals <- read.xlsx("/Users/amandava/Documents/Rochester_AsPIRES/AsPIRES_data/redo/Test3.xlsx",rowNames = FALSE,sheet=2)
p_adj_3 <- apply(test3_pvals[,4:13],1,p.adjust,method="BH")
t_3 <- cbind(test3_pvals[,1:3],t(p_adj_3))

t_3 <- t_3[which(t_3[,1]!="DMSO vs SEB"),]
plot(-log(t_3[,4],base=10),type = "h")
t_3_log <- sapply(t_3[,4:13],function(x){-log(x,base=10)})
t_3_log <- cbind(t_3[,1:3],t_3_log)
melt_t3 <- melt(t_3_log)

M11 <- melt_t3 %>% dplyr::filter(Condition %in% c("Mild") & Visit %in% c("Visit11"))
M12 <- melt_t3 %>% dplyr::filter(Condition %in% c("Mild") & Visit %in% c("Visit12"))
M12 <- M12[which(M12[,1]!="DMSO vs NP5"),]
M13 <- melt_t3 %>% dplyr::filter(Condition %in% c("Mild") & Visit %in% c("Visit13"))

S11 <- melt_t3 %>% dplyr::filter(Condition %in% c("Severe") & Visit %in% c("Visit11"))
S12 <- melt_t3 %>% dplyr::filter(Condition %in% c("Severe") & Visit %in% c("Visit12"))
S13 <- melt_t3 %>% dplyr::filter(Condition %in% c("Severe") & Visit %in% c("Visit13"))

quartz()
ggplot(M11,aes(y=value,x=Stimulation_pairwise,group=variable,fill=variable))+geom_bar(stat = "identity",position="dodge",width = 0.9)+
  geom_hline(yintercept= -log(0.05,base=10),linetype="dotted")+ xlab("Stimulation") + ylab("-log10 Pvalue") + labs(fill="Cell Population")+
  ggtitle("Mild Visit11")

  plots_list[[2]] <- ggplot(M12,aes(y=value,x=Stimulation_pairwise,group=variable,fill=variable))+geom_bar(stat = "identity",position="dodge",width = 0.6) +
    geom_hline(yintercept= -log(0.05,base=10),linetype="dotted") + xlab("Stimulation") + ylab("-log10 Pvalue") + labs(fill="Cell Population")+
    ggtitle("BH corrected P-values for DMSO vs stimulations for Mild, Visit12")+ coord_flip() + theme_grey() +
    theme(axis.text=element_text(size=10),axis.title=element_text(size=10,face="bold"))+
    theme(legend.text=element_text(size=8),legend.title = element_text(size=9,face="bold"))+
    theme(plot.title = element_text(size=11,face="bold"))

#####################
plots_list <- list()
plots_list[[1]] <- ggplot(M11,aes(y=value,x=Stimulation_pairwise,group=variable,fill=variable))+geom_bar(stat = "identity",position="dodge",width = 0.8) +
  geom_hline(yintercept= -log(0.05,base=10),linetype="dotted") + ylab(expression(paste("-Log" ["10"], " BH corrected p-values"))) + labs(fill="Cell Population")+
  ggtitle("Mild-Visit11")+ coord_flip() + theme_grey()+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=8),axis.title.y=element_blank(),axis.ticks.y=element_blank())+
  theme(legend.text=element_text(size=8),legend.title = element_text(size=8,face="bold"),legend.key.size = unit(0.5, "cm"))+
  theme(plot.title = element_text(size=11,face="bold",hjust=0.5))+scale_fill_discrete(guide=guide_legend(reverse=T))+ylim(c(0,3))


plots_list[[2]] <- ggplot(M12,aes(y=value,x=Stimulation_pairwise,group=variable,fill=variable))+geom_bar(stat = "identity",position="dodge",width = 0.4) +
  geom_hline(yintercept= -log(0.05,base=10),linetype="dotted") + ylab(expression(paste("-Log" ["10"], " BH corrected p-values"))) + labs(fill="Cell Population")+
  ggtitle("Mild-Visit12")+ coord_flip() + theme_grey()+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=8),axis.title.y=element_blank(),axis.ticks.y=element_blank())+
  theme(legend.text=element_text(size=8),legend.title = element_text(size=8,face="bold"),legend.key.size = unit(0.5, "cm"))+
  theme(plot.title = element_text(size=11,face="bold",hjust=0.5))+scale_fill_discrete(guide=guide_legend(reverse=T))+ylim(c(0,3))

plots_list[[3]] <- ggplot(M13,aes(y=value,x=Stimulation_pairwise,group=variable,fill=variable))+geom_bar(stat = "identity",position="dodge",width = 0.8) +
  geom_hline(yintercept= -log(0.05,base=10),linetype="dotted") + ylab(expression(paste("-Log" ["10"], " BH corrected p-values"))) + labs(fill="Cell Population")+
  ggtitle("Mild-Visit13")+ coord_flip() + theme_grey()+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=8),axis.title.y=element_blank(),axis.ticks.y=element_blank())+
  theme(legend.text=element_text(size=8),legend.title = element_text(size=8,face="bold"),legend.key.size = unit(0.5, "cm"))+
  theme(plot.title = element_text(size=11,face="bold",hjust=0.5))+scale_fill_discrete(guide=guide_legend(reverse=T))+ylim(c(0,3))

plots_list[[4]] <- ggplot(S11,aes(y=value,x=Stimulation_pairwise,group=variable,fill=variable))+geom_bar(stat = "identity",position="dodge",width = 0.8) +
  geom_hline(yintercept= -log(0.05,base=10),linetype="dotted") + ylab(expression(paste("-Log" ["10"], " BH corrected p-values"))) + labs(fill="Cell Population")+
  ggtitle("Severe-Visit11")+ coord_flip() + theme_grey()+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=8),axis.title.y=element_blank(),axis.ticks.y=element_blank())+
  theme(legend.text=element_text(size=8),legend.title = element_text(size=8,face="bold"),legend.key.size = unit(0.5, "cm"))+
  theme(plot.title = element_text(size=11,face="bold",hjust=0.5))+scale_fill_discrete(guide=guide_legend(reverse=T))+ylim(c(0,3))

plots_list[[5]] <- ggplot(S12,aes(y=value,x=Stimulation_pairwise,group=variable,fill=variable))+geom_bar(stat = "identity",position="dodge",width = 0.8) +
  geom_hline(yintercept= -log(0.05,base=10),linetype="dotted") + ylab(expression(paste("-Log" ["10"], " BH corrected p-values"))) + labs(fill="Cell Population")+
  ggtitle("Severe-Visit12")+ coord_flip() + theme_grey()+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=8),axis.title.y=element_blank(),axis.ticks.y=element_blank())+
  theme(legend.text=element_text(size=8),legend.title = element_text(size=8,face="bold"),legend.key.size = unit(0.5, "cm"))+
  theme(plot.title = element_text(size=11,face="bold",hjust=0.5))+scale_fill_discrete(guide=guide_legend(reverse=T))+ylim(c(0,3))

plots_list[[6]] <- ggplot(S13,aes(y=value,x=Stimulation_pairwise,group=variable,fill=variable))+geom_bar(stat = "identity",position="dodge",width = 0.8) +
  geom_hline(yintercept= -log(0.05,base=10),linetype="dotted") + ylab(expression(paste("-Log" ["10"], " BH corrected p-values"))) + labs(fill="Cell Population")+
  ggtitle("Severe-Visit13")+ coord_flip() + theme_grey()+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=8),axis.title.y=element_blank(),axis.ticks.y=element_blank())+
  theme(legend.text=element_text(size=8),legend.title = element_text(size=8,face="bold"),legend.key.size = unit(0.5, "cm"))+
  theme(plot.title = element_text(size=11,face="bold",hjust=0.5))+scale_fill_discrete(guide=guide_legend(reverse=T))+ylim(c(0,3))

quartz()
plot_grid(plots_list[[1]],plots_list[[2]],plots_list[[3]],plots_list[[4]],plots_list[[5]],plots_list[[6]],label_size = 0.2)

##HeatMap
m11_hm <- M11[,c(1,4,5)] %>% dcast(variable~Stimulation_pairwise,value.var="value")
rownames(m11_hm) <- m11_hm[,1]
m11_hm <- m11_hm[,-1]
m11_hm <- as.matrix(m11_hm)
heatmap.2(as.matrix(m11_hm),col = bluered(100), trace = "none", Rowv = F,
                 Colv = F, scale = "row",
                 mar = c(10, 20))

mild <- rbind(M11,M12,M13)
m123 <- mild[,c(1,3,4,5)] %>% dcast(variable+Visit~Stimulation_pairwise,value.var = "value",fun.aggregate = NULL)
mild_123 <- m123[,c(1,2,4,5,6,8)]

#####Redo the models by merging the mild and severe groups
##model f: without sample level to account for overdispersion
genmix_d <- genmix_data[which(genmix_data$Stimulation != "F1" & genmix_data$Stimulation != "NP1" & genmix_data$Stimulation != "SEB"),-3]
#genmix_d <- filter(genmix_data,Stimulation != "F1" & Stimulation != "NP1" & Stimulation != "SEB")
genmix_d$Stimulation <- factor(genmix_d$Stimulation,levels=c("DMSO","F5","G5","M10","NP5"))
#genmod_1.f <- genmix_d[,12]/genmix_d[,10] ~ Stimulation * Visit + (1|SubjectID)+(1|Sample)
#genfit_1.f <- glmer(genmod_1.f,weights=genmix_d[,10],family="binomial",data = genmix_d)

Test_3 <- list(list())
for (i in 1:ncol(genmix_d)) {
  if(i>11 & i< 22){
    genmod_1.f <- genmix_d[,i]/genmix_d[,10] ~ Stimulation * Visit + (1|SubjectID)+(1|Sample)
    genfit_1.f <- glmer(genmod_1.f,weights=genmix_d[,10],family="binomial",data = genmix_d)
    ls_3 <- emmeans(genfit_1.f,~Stimulation| Visit)
    Test_3[[i-11]] <- contrast(ls_3,interaction = "trt.vs.ctrl",adjust="Dunnett")
    quartz()
    plot_emm <- emmip(genfit_1.f,~Stimulation|Visit,type="response",CI=T)+theme_gray()+ggtitle(paste0("Mean and 95% CI for ",colnames(genmix_d)[i]))
    print(plot_emm)
    p <- genmix_d[,i]/genmix_d[,10]
    quartz()
    plot_box <- ggplot(data=genmix_d)+geom_boxplot(aes(x=Stimulation,y=p,color=Stimulation))+theme_gray()+
      ylab("Cell proportion")+facet_wrap(~Visit)+ggtitle(paste0("Boxplot for ",colnames(genmix_d)[i]))
    print(plot_box)
  }else if(i>21){
    genmod_1.f <- genmix_d[,i]/genmix_d[,11] ~ Stimulation * Visit + (1|SubjectID)+(1|Sample)
    genfit_1.f <- glmer(genmod_1.f,weights=genmix_d[,11],family="binomial",data = genmix_d)
    ls_3 <- emmeans(genfit_1.f,~Stimulation| Visit)
    Test_3[[i-11]] <- contrast(ls_3,interaction = "trt.vs.ctrl",adjust="Dunnett")
    quartz()
    plot_emm <- emmip(genfit_1.f,~Stimulation|Visit,type="response",CI=T)+theme_gray()+ggtitle(paste0("Mean and 95% CI for ",colnames(genmix_d)[i]))
    print(plot_emm)
    p <- genmix_d[,i]/genmix_d[,10]
    quartz()
    plot_box <- ggplot(data=genmix_d)+geom_boxplot(aes(x=Stimulation,y=p,color=Stimulation))+theme_gray()+
      ylab("Cell proportion")+facet_wrap(~Visit)+ggtitle(paste0("Boxplot for ",colnames(genmix_d)[i]))
    print(plot_box)
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

quartz()
ggplot(data = melted_p, mapping = aes(x = Contrast,y = variable,
                                                       fill = pvalue))+
  facet_grid(~Visit, scales = "free_x", space = "free_x",drop=T)+geom_tile()+
scale_fill_gradient2('pvalue', limits=c(0,1), breaks = c(0,0.05,0.3,0.5,0.75,0.9, 1),
                     low = "#d53e4f",mid="#f46d43",high = "#fdae61",midpoint = 0.5)+
  theme_bw()+
  theme(strip.placement = "outside",plot.title = element_text(hjust = 0.4),
        axis.title.y = element_blank(), # Remove y-axis title
        strip.text = element_text(face="bold", size=9),
        strip.background = element_blank(),
        axis.text.x = element_text(size=7,angle=25,face="bold"), axis.text.y = element_text(size=7,face="bold"),
        legend.text=element_text(size=6,face="bold"),legend.key.width=grid::unit(0.3,"cm"),legend.title = element_text(size=8,
                                                                                 face="bold"),
        panel.border = element_blank(),plot.background=element_blank(),panel.spacing=unit(0, "lines"))+
  ggtitle(label = "DMSO vs other stimulations across visits")+
  scale_y_discrete(limits = rev(levels(as.factor(melted_p$variable))))+xlab("Stimulations")+ylab("Cell-population")

p_adj_f[,3:13] <- -log10(p_adj_f[,3:13])
pp <- p_adj_f
colnames(pp)[3:13] <- c("CD4T/CD69+TNFa+","CD4T/CD69+IL2+","CD4T/CD69+IL4+",
                        "CD4T/CD69+IL10+","CD4T/CD69+IL17+","CD4T/CD69+IFNg+","CD4T/CD25+FoxP3+_bi",
                        "CD4T/CD25+FoxP3+", "CD4T/CD25+FoxP3_IL10+","CD8T/CD69+TNFa+",
                        "CD8T/IFNg+"  )
pp <- pp[order(pp$Contrast),]
pp_n <- pp[,c(3:8,10,11,13)]
#write.xlsx(pp,"/Users/amandava/Documents/Rochester_AsPIRES/AsPIRES_data/AsPIRES_neglog10_feb26.xlsx",rownames=T)
#rownames(pp_n) <- pp[,1]
#pp[,1] <- paste0(rep(c(1,2,3),each=4),pp[,1])
anno <- pp[,1:2]
annotation <- (anno %>% separate(Contrast,c("Contrast","DMSO")))[,c(1:3)]
t <- t(pp_n) > 1.30103
t[t=="TRUE"] <- "*"
t[t=="FALSE"] <- " "
#t[8,4] <- t[5,4] <- t[3,1] <- "**"

quartz()
aheatmap(t(pp_n),annCol=list(Stimulation=annotation$Contrast,Visit=annotation$Visit),Colv=NA,
         Rowv= NA,cellwidth = 16, cellheight = 18, fontsize = 8,labCol = NA,color = "-RdYlBu2",
         breaks=seq(0,5,0.1),main = "-log10 BH corrected p-values",txt=t,annLegend=c("**: Significant increase from DMSO","*: Significant difference from DMSO"))
quartz()
aheatmap(t(pp_n),annCol=list(Stimulation=annotation$Contrast,Visit=annotation$Visit),Colv=NA,
         Rowv= NA,cellwidth = 16, cellheight = 18, fontsize = 8,labCol = NA,color = 1L,
         breaks=seq(0,5,0.1),main = "-log10 BH corrected p-values",txt=t)

##Fold change heatmap
for (i in 1:ncol(genmix_d)) {
  if(i>11 & i< 22){
    sub_data <- cbind(genmix_d[,1:5],genmix_d[,i]/genmix_d[,10])
    colnames(sub_data)[6] <- colnames(genmix_d)[i]
    wide <- dcast(sub_data,SubjectID+Visit~Stimulation)
    foldchange <-

    genfit_1.f <- glmer(genmod_1.f,weights=genmix_d[,10],family="binomial",data = genmix_d)
    ls_3 <- emmeans(genfit_1.f,~Stimulation| Visit)
    Test_3[[i-11]] <- contrast(ls_3,interaction = "trt.vs.ctrl",adjust="Dunnett")
  }else if(i>21){
    genmod_1.f <- genmix_d[,i]/genmix_d[,11] ~ Stimulation * Visit + (1+Visit|SubjectID)+(1|Sample)
    genfit_1.f <- glmer(genmod_1.f,weights=genmix_d[,11],family="binomial",data = genmix_d)
    ls_3 <- emmeans(genfit_1.f,~Stimulation| Visit)
    Test_3[[i-11]] <- contrast(ls_3,interaction = "trt.vs.ctrl",adjust="Dunnett")
  }
}

color = colorRampPalette(c("yellow","lightgreen","lightblue","steelblue"))(50)
