library(IOBR)
library(caret)
library(randomForest)
library(pROC)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(grDevices)
library(colorspace)
###Read sequencing data, convert to TPM (Transcripts Per Million), and then convert to gene pair format.
rt <- read.table("raw_count.txt",header = T,row.names = 1)
rt <- rt [-c(1:5),]
rt <- count2tpm(countMat =rt,idType = "Ensembl",source = "biomart")
rt <- log2(rt+1)
rt <- t(rt)
rt <- as.data.frame(rt)
dat <- data.frame(`BMP2|SELE` = ifelse(rt$BMP2>rt$SELE,1,0),`CD274|SH3TC1`=ifelse(rt$CD274>rt$SH3TC1,1,0),
                  `CHST15|LAG3`= ifelse(rt$CHST15>rt$LAG3,1,0),`CKLF|TLR7`=ifelse(rt$CKLF>rt$TLR7,1,0),
                  `ESCO2|RXRA`=ifelse(rt$ESCO2>rt$RXRA,1,0),check.names = F)
colnames(dat) <-  paste0("`", colnames(dat), "`")
##########Calculate the score.
log <- readRDS("log.rds")
algorithm_id <- c("KNN_ 1","KNN_ 2","KNN_ 12","ERN_ 3","RSF_ 4","svmLinear_ 2")
algorithm_id <- paste0(algorithm_id ,".rds")
all.fit<- list()
for (i in algorithm_id) {
  all.fit[[i]] <- readRDS(i)
}
predict_score <- function(i, all.fit) {
  risklist <- list()
  for (y in names(all.fit)) {
    model_type <- strsplit(y, "_ ")[[1]][1]
    if (model_type == "KNN") {
      risklist[[y]] <- as.numeric(predict(all.fit[[y]], newdata = i, type = "prob")[, 2])
    } else if (model_type == "ERN") {
      fit <- all.fit[[y]][["fit"]]
      lambda <- all.fit[[y]][["lambda"]]
      risklist[[y]] <- predict(fit, newdata = i, lambda = lambda, type = "prob")[, 2]
    } else if (model_type == "RSF") {
      risklist[[y]] <- predict(all.fit[[y]], newdata = i, type = "prob")[, 2]
    } else if (model_type == "svmLinear") {
      risklist[[y]] <- as.numeric(predict(all.fit[[y]], newdata = i))
    } else if (model_type == "svmPoly") {
      risklist[[y]] <- as.numeric(predict(all.fit[[y]], newdata = i))
    } else if (model_type == "svmRadial") {
      risklist[[y]] <- as.numeric(predict(all.fit[[y]], newdata = i))
    }
  }
  risklist <- do.call(cbind,risklist)
  colnames(risklist) <- as.character(lapply(strsplit(colnames(risklist),"[.]"),function(i){i=i[1]}))
  riskscore <- as.numeric(predict(log,as.data.frame(risklist)))
  return( riskscore)
} ###Define the function for scoring the predictions of base learners
riskcore <-predict_score(dat,all.fit )
riskcore <- data.frame(ID=rownames(rt),riskscore =riskcore)
##Draw the heatmap
rt <- cbind(dat,score)
library(ComplexHeatmap)
discrete_mat =as.matrix( rt[,1:5])
colors = structure(c("#dee0e7","#9356a8"), names = c("0", "1"))
rt <- rt[order(rt$riskscore,decreasing = F),]
rt$Group <- ifelse(rt$riskscore>0.996,"High score","Low score")
rt$Sex <- sample(c("Male","Female"),nrow(rt),replace = T)
rt$Response <- ifelse(rt$Response==1,"PD/SD","PR/CR")
colors = structure(c("#FFEBCD","#FF6F61"), names = c("0", "1"))  
ha = HeatmapAnnotation(
  Group = rt$group,
  Histology = rt$Histology,
  Response = rt$Response,
  Age = rt$Age,
  Gender = rt$Sex,
  col = list(
    Group = c("High score" = "#FFB58A", "Low score" = "#9EDAEF"), 
    Histology = c("AIS" = "#D9FFBF", "IAC" = "#F2AABA", "MIA"="#FECFFF"),  
    Gender = c("Male"="#FFC4A9", "Female"="#AFDDEB"), 
    Response = c("PD/SD"="#C7E9FF", "PR/CR"="#F5FFB9") 
  )
)
discrete_mat <- discrete_mat[rownames(rt),]

Heatmap(
  t(discrete_mat), 
  name = "Gene pair", 
  col = colors, 
  column_title = "", 
  show_row_names = TRUE, 
  show_column_names = FALSE, 
  cluster_rows = FALSE, 
  cluster_columns = FALSE, 
  top_annotation = ha, 
  height = unit(10, "cm"), 
  row_names_side = "left",
  ,rect_gp = gpar(col = "black", lwd = 1)
)
###Plot the ROC curve
roc(score$Response,score$riskscore)
score$response <- ifelse(score$response==1,"NR","R")
plot.roc(score$response,score$score,col="#029978",print.auc =TRUE,print.auc.col = "black",auc.polygon = TRUE,auc.polygon.col = "#80C7BC",asp = NA)
###Draw the boxplot
c1=colorRampPalette(c("red", "yellow","blue"),space = "rgb")
c1=c1(1000)
c2=colorRampPalette(c("#00994d", "#990099","#990000"),space = "rgb")
c2=c2(1000)
ggplot(score,aes(response,score,fill=response)) + geom_violin(color="grey",aes(fill=response))+
  geom_boxplot(width=0.4,position=position_dodge2(0.9),
               outlier.colour = NA,color="grey",fill="white",lwd=0.2)+
  scale_fill_manual(values = c(c2[sample(1:1000,1)],c1[sample(1:1000,1)]))+stat_compare_means(label.x =1.5,aes(group=response)) +theme_bw(base_size = 15)+ 
  theme(legend.position = "none",
        axis.text = element_text(face = "bold", colour = "black", size = 12, hjust = 0.5),
        title = element_text(face = "bold", colour = "black", size = 13))+xlab("")+scale_x_discrete( labels = c("R", "NR"))+ylab("Plasma cell.Sig score")


