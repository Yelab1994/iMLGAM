library(data.table)
library(tibble)
library(qs)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
###Read CRISPR data
load('crispr.Rdata')
df <- crispr
df$mean <- apply(df,1,function(z){mean(z,na.rm=T)})
df <- df[order(df$mean),]
df <- df[c(1:8, (nrow(df)-7):nrow(df)),]
range(df,na.rm=T)
column_ha = HeatmapAnnotation(cohort = colnames(df))
##Draw the heatmap
Heatmap(as.matrix(df), name = "z scores",
        column_title = NULL,
        row_title = NULL,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = colorRamp2(c(-2, 0, 2), c("#377EB8", "white", "#E41A1C")), # 新配色
        show_row_names = TRUE, 
        show_column_names = FALSE,
        rect_gp = gpar(col = "black", lwd = 2),
        width = ncol(df) * unit(5, "mm"), 
        height = nrow(df) * unit(5, "mm"),
        na_col = 'white',
        column_names_side = 'top',
        row_split = c(rep('a', 8), rep('b', 8)),
        column_split = c(rep('a', 17), 'b'),
        top_annotation = column_ha
)
###ranking of genes###
rank <- apply(crispr,1,function(z){ mean(z,na.rm=T)})
rank <- data.frame(meanZ=rank,row.names = rownames(crispr))
rank$genes <- rownames(crispr)
rank <- rank[order(rank$meanZ),]
rank$order <- order(rank$meanZ)
