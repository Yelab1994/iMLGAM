## ── 1. Load all required packages ---------------------------------------------
library(data.table)      # fast data manipulation
library(tibble)          # tidy data frames
library(qs)              # quick serialization / de-serialization
library(ComplexHeatmap)  # flexible heatmap generation
library(circlize)        # circular utilities (used by ComplexHeatmap)
library(ggpubr)          # ggplot2 extensions for publication
library(ggplot2)         # the grammar of graphics
library(data.table)      # repeated
library(tibble)          # repeated
library(qs)              # repeated
library(survival)        # survival analysis

## ── 2. Load CRISPR data and create a small demo heatmap -----------------------
load("crispr.Rdata")                      # loads object 'crispr' (matrix or DF)
df <- crispr                              # rename for convenience

# compute row-wise mean z-scores (ignoring NA)
df$mean <- apply(df, 1, function(z) mean(z, na.rm = TRUE))

# order rows by mean z-score
df <- df[order(df$mean), ]

# keep the 8 lowest- and 8 highest-mean rows (16 genes total)
df <- df[c(1:8, (nrow(df) - 7):nrow(df)), ]

# quick check of value range
range(df, na.rm = TRUE)

# column annotation: simple cohort label for each sample
column_ha = HeatmapAnnotation(cohort = colnames(df))

# draw the heatmap
Heatmap(
  as.matrix(df),
  name = "z scores",
  column_title = NULL,
  row_title    = NULL,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = colorRamp2(c(-2, 0, 2), c("#377EB8", "white", "#E41A1C")),  # blue-white-red
  show_row_names = TRUE,
  show_column_names = FALSE,
  rect_gp = gpar(col = "black", lwd = 2),        # thick black border for each cell
  width  = ncol(df) * unit(5, "mm"),             # 5 mm per column
  height = nrow(df) * unit(5, "mm"),             # 5 mm per row
  na_col = 'white',                              # color for missing values
  column_names_side = 'top',
  row_split    = c(rep('a', 8), rep('b', 8)),    # split 8 lowest vs 8 highest
  column_split = c(rep('a', 17), 'b'),           # split columns arbitrarily
  top_annotation = column_ha
)

## ── 3. Rank all genes by their mean z-score across samples --------------------
rank <- apply(crispr, 1, function(z) mean(z, na.rm = TRUE))  # vector of means
rank <- data.frame(meanZ = rank, row.names = rownames(crispr))  # turn into DF
rank$genes <- rownames(crispr)      # extra gene-name column
rank <- rank[order(rank$meanZ), ]   # order ascending by mean z-score
rank$order <- order(rank$meanZ)     # numeric rank

## ── 4. Load ICI clinical data and compute univariate Cox models ---------------
rt <- readRDS("ICI_exp_sva.rds")               # expression + clinical data
group <- readRDS("ICI_riskscore.rds")      # risk scores from earlier model
group <- rbind(group$trainset_7,
               group$trainset_3,
               group$testset)              # merge all risk-score sets

jj <- intersect(rownames(rt), rownames(group))  # samples in common
rt    <- rt[jj, ]
group <- group[jj, ]

# create binary risk group
group$Risk <- ifelse(group$riskscore > 0.996, "High", "Low")

rt <- cbind(Risk = group$Risk, rt)  # append risk label
rt <- rt[rt$Risk == "High", ]       # restrict to High-risk samples

## 4a. Perform univariate Cox regression for every gene ------------------------
outTab <- data.frame()            # results container
sigGenes <- c("futime", "fustat") # survival time & event columns

for (i in colnames(rt[, 9:ncol(rt)])) {        # iterate over genes (cols 9+)
  dat <- rt[, c("futime", "fustat", "Type", i)]
  colnames(dat)[4] <- "gene"                    # rename gene column
  cox <- coxph(Surv(futime, fustat) ~ ., data = dat)  # Cox model
  coxSummary <- summary(cox)
  
  # extract hazard ratio and 95% CI
  outTab <- rbind(outTab,
                  cbind(id     = i,
                        HR     = coxSummary$conf.int["gene", "exp(coef)"],
                        HR.95L = coxSummary$conf.int["gene", "lower .95"],
                        HR.95H = coxSummary$conf.int["gene", "upper .95"],
                        pvalue = coxSummary$coefficients["gene", "Pr(>|z|)"]))
}

outTab$pvalue <- as.numeric(outTab$pvalue)
outTab <- outTab[outTab$pvalue < 0.001, ]  # keep genes with p < 0.001

## ── 5. Merge Cox results with the ranking table ------------------------------
rank$id <- rownames(rank)                  # ensure gene ID column exists
library(dplyr)
rank <- inner_join(rank, outTab, by = "id")  # merge by gene symbol

## ── 6. Plot top 20 significant genes ordered by mean z-score -----------------
rank <- rank[1:20, ]                       # keep top 20 genes
rank$id <- factor(rank$id, levels = rank$id)  # lock factor order for ggplot

ggplot(rank, aes(x = id, y = order, group = 1)) +
  geom_line(color = "#1C4E80", size = 1) +   # line connecting points
  geom_point(color = "#E41A1C", size = 2.5) +# red points
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(face = "bold.italic", colour = "#441718",
                             size = 12, hjust = 0.5),
    axis.line = element_blank(),
    plot.title = element_text(face = "bold.italic", colour = "#441718",
                              size = 16, hjust = 0.5),
    title = element_text(face = "bold.italic", colour = "#441718",
                         size = 13),
    panel.border = element_rect(fill = NA, color = "black",
                                size = 1.5, linetype = "solid"),
    panel.background = element_rect(fill = "#F1F6FC"),
    panel.grid.major = element_line(color = "#CFD3D6",
                                    size = .5, linetype = "dotdash")
  ) +
  xlab("") +
  ylab("Order") +
  coord_flip()   # rotate so genes are on y-axis