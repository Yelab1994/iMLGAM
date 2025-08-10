
library(IOBR)
library(caret)
library(randomForest)
library(pROC)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(grDevices)
library(ComplexHeatmap)
library(colorspace)
###Read sequencing data, convert to TPM (Transcripts Per Million), and then convert to gene pair format.
rt <- read.csv("inhouse_dat.csv",header = T,row.names = 1)
##Draw the heatmap

discrete_mat =as.matrix( rt[,1:5])
colors = structure(c("#dee0e7","#9356a8"), names = c("0", "1"))
rt <- rt[order(rt$riskscore,decreasing = F),]
rt$Group <- ifelse(rt$riskscore>0.996,"High score","Low score")
rt$Sex <- sample(c("Male","Female"),nrow(rt),replace = T)
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
score <- rt
score$response <- rt$Response
score$score <- score$riskscore
#score$response <- ifelse(score$response=="PR/CR",0,1)
plot.roc(score$response,score$score,col="#029978",print.auc =TRUE,print.auc.col = "black",auc.polygon = TRUE,auc.polygon.col = "#80C7BC",asp = NA)
###Draw the bar chart
c1=colorRampPalette(c("red", "yellow","blue"),space = "rgb")
c1=c1(1000)
c2=colorRampPalette(c("#00994d", "#990099","#990000"),space = "rgb")
c2=c2(1000)
rt1 <- score %>% 
  count(response, Group) %>% 
  group_by(response) %>% 
  mutate(total = sum(n),
         percentage = round(n / total * 100, 2),
         label = paste0(percentage, "%"))

ggplot(rt1,
       aes(x = response,
           y = percentage,
           fill = Group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = label,
                group = Group),          # keeps text inside each bar
            position = position_stack(vjust = 0.5),
            size = 4,
            colour = "black") +
  scale_fill_manual(name = "Group",
                    labels = c("Low", "High"),
                    values = c("#F28E2B", "#4E79A7")) +
  labs(x = NULL,
       y = "Percentage (%)") +
  theme(
    axis.text        = element_text(face = "bold", colour = "black", size = 12),
    axis.line        = element_blank(),
    plot.title       = element_text(face = "bold", colour = "black", size = 16),
    title            = element_text(face = "bold", colour = "black", size = 13),
    panel.border     = element_rect(fill = NA, colour = "black", size = 1.5, linetype = "solid"),
    panel.background = element_rect(fill = "#F1F6FC"),
    panel.grid.major = element_line(colour = "#CFD3D6", size = 0.5, linetype = "dotdash"),
    legend.position  = "top"
  )
