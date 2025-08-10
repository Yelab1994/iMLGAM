## ── Load required packages ------------------------------------------
library(data.table)
library(tibble)
library(pROC)      # ROC / AUC utilities
library(ggplot2)   # plotting
library(IOBR)      # ssGSEA / deconvolution helper functions
library(easier)    # immune-response signatures
library(patchwork) # arranging multiple ggplot objects

## ── Import SVA-corrected expression matrix -------------------------
rt_sva <- readRDS("ICI_exp_sva.rds")
# File shared via cloud storage: ICI_exp_sva.rds
# Link: https://pan.baidu.com/s/1QGOW3XA4JBsEOXWMEk-sNg?pwd=3gsm
# Extraction code: 3gsm

## ── Import cohort list (train/test splits) -------------------------
cohort <- qs::qread("ICI_data.qs")

## ── Subset SVA matrix to common genes for each split ---------------
cohort$trainset_7 <- rt_sva[rownames(cohort$trainset_7), ]
cohort$trainset_3 <- rt_sva[rownames(cohort$trainset_3), ]
cohort$testset    <- rt_sva[rownames(cohort$testset),   ]

## ── CEP (Cytotoxic Effector Profile) -------------------------------
genes <- c("IDO1","CD3D","IL2RG","CXCL10","NKG7","CXCL9","CIITA","HLA-E",
           "HLA-DRA","CD3E","CXCR6","STAT1","LAG3","IFNG","GZMK","TAGAP","CD2")

## Check for any missing genes (should return character(0))
setdiff(genes, colnames(cohort$trainset_7))

## Compute CEP score (row-wise sum) for each split
for (i in names(cohort)) {
  dat           <- cohort[[i]]
  cep_score     <- rowSums(dat[, genes, drop = FALSE])   # sum of gene expression
  cohort[[i]]   <- cbind(CEP = cep_score[rownames(dat)], dat) # prepend to matrix
}

## Evaluate CEP performance via ROC
roc(cohort$trainset_7$Res, cohort$trainset_7$CEP)
roc(cohort$trainset_3$Res, cohort$trainset_3$CEP)
roc(cohort$testset$Res,    cohort$testset$CEP)

## ── PD-L1 (CD274) ---------------------------------------------------
roc(cohort$trainset_7$Res, cohort$trainset_7$CD274)
roc(cohort$trainset_3$Res, cohort$trainset_3$CD274)
roc(cohort$testset$Res,    cohort$testset$CD274)

## ── CAF (Cancer-Associated Fibroblast) -----------------------------
library(sparrow)   # eigenWeightedMean()

gene_vector <- c("MMP11","COL11A1","C1QTNF3","CTHRC1","COL12A1","COL10A1",
                 "COL5A2","THBS2","AEBP1","LRRC15","ITGA11")

for (i in names(cohort)) {
  dat         <- cohort[[i]]
  caf_score   <- eigenWeightedMean(dat[, gene_vector])$weights  # PCA-based score
  cohort[[i]] <- cbind(CAF = as.numeric(caf_score), dat)
}

roc(cohort$trainset_7$Res, cohort$trainset_7$CAF)
roc(cohort$trainset_3$Res, cohort$trainset_3$CAF)
roc(cohort$testset$Res,    cohort$testset$CAF)

## ── NLRP3 inflammasome signature (ssGSEA) --------------------------
all_data <- rbind(cohort$trainset_7,
                  cohort$trainset_3,
                  cohort$testset)

gene_vector <- c("ARRDC1-AS1","CARD8","GSDMD","ATAT1","CD36","CPTP","DHX33",
                 "EIF2AK2","GBP5","NLRC3","PYDC2","SIRT2","TLR4","TLR6","USP50",
                 "APP","CASP1","HSP90AB1","MEFV","NFKB1","NFKB2","NLRP3","P2RX7",
                 "PANX1","PSTPIP1","PYCARD","RELA","SUGT1","TXN","TXNIP")

dat <- t(all_data[, -c(1:10)])   # genes in rows, samples in columns

NLPR3 <- calculate_sig_score_ssgsea(
  eset         = dat,
  signature    = list(gene_vector),
  mini_gene_count = 3
)
NLPR3 <- NLPR3$V1
all_data <- cbind(NLPR3 = NLPR3, all_data)

## ── IPS (Immunophenoscore) -----------------------------------------
IPS <- deconvo_ips(dat, plot = FALSE)
all_data <- cbind(IPS = IPS$IPS_IPS, all_data)
roc(all_data$Res, all_data$IPS)

## ── Easier immune-response hallmarks -------------------------------
hallmarks_of_immune_response <- c("CYT","Roh_IS","chemokines","Davoli_IS",
                                  "IFNy","Ayers_expIS","Tcell_inflamed","RIR","TLS")

immune_response_scores <- compute_scores_immune_response(
  RNA_tpm         = dat,
  selected_scores = hallmarks_of_immune_response
)

jj <- intersect(rownames(all_data), rownames(immune_response_scores))
all_data               <- all_data[jj, ]
immune_response_scores <- immune_response_scores[jj, ]

all_data <- cbind(immune_response_scores, all_data)
all_data <- all_data[, -c(5, 9)]   # drop duplicated / unwanted columns

## ── TIDE score ------------------------------------------------------
tide <- read.csv("tide.csv", header = TRUE, row.names = 1, check.names = FALSE)
rownames(tide) <- gsub("[.]", "-", rownames(tide))   # harmonize sample IDs

jj <- intersect(rownames(all_data), rownames(tide))
all_data <- all_data[jj, ]
tide     <- tide[jj, ]
all_data <- cbind(TIDE = tide$TIDE, all_data)
roc(all_data$Res, all_data$TIDE)

## ── Merge back risk-score calculated elsewhere ---------------------
cohort <- readRDS("ICI_riskscore.rds")

cohort$trainset_7 <- cbind(riskscore = cohort$trainset_7$riskscore,
                           all_data[rownames(cohort$trainset_7), ])
cohort$trainset_3 <- cbind(riskscore = cohort$trainset_3$riskscore,
                           all_data[rownames(cohort$trainset_3), ])
cohort$testset    <- cbind(riskscore = cohort$testset$riskscore,
                           all_data[rownames(cohort$testset), ])

## ── Keep only first 13 columns + CD274 -----------------------------
for (i in names(cohort)) {
  cohort[[i]] <- cohort[[i]][, c(1:13, which(colnames(cohort$trainset_7) == "CD274"))]
}

## ── Compute AUC for every signature in each split ------------------
auc_list <- list()

for (i in names(cohort)) {
  auc_dat <- data.frame()
  dat     <- cohort[[i]]
  for (y in colnames(dat)[-length(colnames(dat))]) {   # skip last column (CD274)
    dat <- cbind(Res = rt_sva[rownames(dat), "Res"], dat)
    auc_test <- roc(dat[, "Res"], dat[, y])
    auc_dat  <- rbind(auc_dat,
                      data.frame(id = y, AUC = as.numeric(auc_test$auc)))
  }
  auc_list[[i]] <- auc_dat
}

## ── Generate bar-plot for each split -------------------------------
pic_list <- list()

for (i in names(auc_list)) {
  dat <- auc_list[[i]]
  dat$AUC <- as.numeric(dat$AUC)
  dat <- dat[order(dat$AUC, decreasing = TRUE), ]
  dat$id <- factor(dat$id, levels = dat$id)
  
  pic_list[[i]] <- ggplot(dat, aes(x = id, y = AUC, fill = id)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = "AUC") +
    scale_fill_manual(values = c("#E69F00","#56B4E9","#009E73","#F0E442",
                                 "#0072B2","#D55E00","#CC79A7","#999999",
                                 "#D73027","#1A9850","#31688E","#A50026",
                                 "#542788")) +
    theme(
      axis.text         = element_text(face = "bold.italic",
                                       colour = "#441718", size = 10),
      axis.title.y      = element_text(face = "bold.italic",
                                       colour = "#441718", size = 16),
      axis.line         = element_blank(),
      plot.title        = element_text(face = "bold.italic",
                                       colour = "#441718", size = 16, hjust = 0.5),
      legend.text       = element_text(face = "bold.italic"),
      panel.border      = element_rect(fill = NA, colour = "black",
                                       size = 1.5, linetype = "solid"),
      panel.background  = element_rect(fill = "white"),
      panel.grid.major  = element_line(colour = "#CFD3D6",
                                       size = .5, linetype = "dotdash"),
      legend.title      = element_text(face = "bold.italic", size = 13),
      legend.position   = "none"
    )
}

## ── Combine plots vertically and save ------------------------------
combined_plot <- pic_list$trainset_7 +
  pic_list$trainset_3 +
  pic_list$testset +
  plot_layout(ncol = 1)

ggsave("AUC比较.pdf", combined_plot,
       units = "mm", height = 180 * 1.3, width = 105 * 1.3)