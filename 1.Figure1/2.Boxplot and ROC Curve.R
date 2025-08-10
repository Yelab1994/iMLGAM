
library(survival)
library(survminer)
library(pROC)
###Read the risk scores of various Cohort
bas_learn<- readRDS("riskscore.rds")
names(bas_learn) <- c("Training Set","Internal Validation Set","External Validation Set")
###Receiver Operating Characteristic curve plot
for (i in names(bas_learn)) {   dat <- bas_learn[[i]]
dat <- dat[!is.na(dat$res),]
roc_obj <- roc(dat$res, dat$riskscore,smooth=F) 
roc_coords <- coords(roc_obj)
ROC_1 <- ggplot(data = roc_coords, aes(x = 1-specificity, y = sensitivity)) +
  geom_ribbon(aes(ymin = 0, ymax =  sensitivity), fill = "#029978", alpha = 0.2)+
  geom_line(color = "#029978",lwd=1)  +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black", size = 1)+
  labs(title = "ROC Curve",
       x = "1 - Specificity",
       y = "Sensitivity") +
  annotate("text",x=0.7,y=0.3,label=paste0("AUC = ",round(as.numeric(roc_obj $auc),3)),
           colour = "black", size = 5, fontface = 2)+
  theme(
    axis.text = element_text(face = "plain", colour = "black", size = 12, hjust = 0.5),
    axis.line = element_blank(),
    plot.title = element_text(face = "bold", colour = "black", size = 16, hjust = 0.5),
    title = element_text(face = "bold", colour = "black", size = 13),
    panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid"),
    panel.background = element_rect(fill = "#F1F6FC"),
    panel.grid.major = element_line(color = "#CFD3D6", size = .5, linetype = "dotdash"),
    legend.position = "none"
  )
ggsave(paste0(i,"_roc.pdf"),units = "mm",width = 80,height =80)
dat$res <- ifelse(dat$res==1,"Non-res","Res")
boxplot_1 <- ggplot(data = dat, aes(x = res, y = riskscore, fill = res)) +     ##需要修改
  scale_fill_manual(values = c("#C31820", "#106C61")) +
  geom_boxplot(notch = FALSE, outlier.size = 1, color = "black", lwd = 0.8, alpha = 1,width=0.8,outlier.colour = NA) +
  theme_classic() +
  ylab('riskScore') +
  xlab("") +
  ggtitle(i) +
  theme(
    axis.text = element_text(face = "plain", colour = "black", size = 12, hjust = 0.5),
    axis.line = element_blank(),
    plot.title = element_text(face = "bold", colour = "black", size = 16, hjust = 0.5),
    title = element_text(face = "bold", colour = "black", size = 13),
    panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid"),
    panel.background = element_rect(fill = "#F1F6FC"),
    panel.grid.major = element_line(color = "#CFD3D6", size = .5, linetype = "dotdash"),
    legend.position = "none"
  ) +
  stat_compare_means(label.x =1.5,aes(group=res),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
ggsave(paste0(i,"_boxplot.pdf"),units = "mm",width = 80,height = 80)
}
