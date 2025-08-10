################################################################################
## 0.  Working directory & packages
################################################################################
library(limma)                                      # Linear models for microarray/RNA-seq
library(ggplot2)                                    # Grammar-of-graphics plotting
library(ggrepel)                                    # Non-overlapping text labels
library(gsea)                                       # Gene-set enrichment utilities
library(clusterProfiler)                            # Functional enrichment wrapper
library(org.Hs.eg.db)                               # Human gene annotation DB
library(grDevices)                                  # Base graphic devices/colours
library(colorspace)                                 # Colour utilities
library(qs)                                         # Fast .qs file reader (for qread)

################################################################################
## 1.  Load sample metadata (risk groups) & TPM expression matrix
################################################################################
group <- readRDS("TCGA_group.rds")      # Data frame: sample → Risk (High / Low)
rt    <- qread("TCGA_tpm.qs")           # # File shared via cloud storage: TCGA_tpm.qs
# Link: https://pan.baidu.com/s/1_j_LztrdC5ZErq8_8ImJZw?pwd=np63
# Extraction code: np63

################################################################################
## 2.  Filter tumour samples (01–09) and log2-transform
################################################################################
mat <- rt[, as.numeric(substr(colnames(rt), 14, 15)) < 10]   # Keep tumour barcodes
colnames(mat) <- substr(colnames(mat), 1, 12)                # Trim to TCGA-XX-XXXX-01
mat <- log2(mat + 1)                                         # log2(TPM + 1)

################################################################################
## 3.  Match expression matrix with risk labels
################################################################################
mat  <- t(mat)                           # transpose → samples × genes
jj   <- intersect(rownames(mat), rownames(group))   # common samples
mat  <- mat[jj, ]                        # subset expression
group <- group[jj, ]                     # subset metadata

################################################################################
## 4.  limma differential-expression: Low-risk vs High-risk
################################################################################
list <- factor(group$Risk)               # Create factor: High / Low
design <- model.matrix(~0 + list)        # Design matrix without intercept
colnames(design) <- c("High", "Low")     # Rename coefficients

contrast <- makeContrasts(Low - High, levels = design)  # Contrast Low – High
fit      <- lmFit(t(mat), design)        # Fit linear model
fit      <- contrasts.fit(fit, contrast) # Apply contrast
fit      <- eBayes(fit)                  # Empirical Bayes moderation
tempOutput <- topTable(fit, n = Inf, adjust = "fdr")    # All genes + FDR
tempOutput$SYMBOL <- rownames(tempOutput)               # Preserve gene symbol

################################################################################
## 5.  Convert gene symbols → Entrez IDs for GSEA
################################################################################
df <- bitr(geneID = tempOutput$SYMBOL,
           fromType = "SYMBOL",
           toType   = "ENTREZID",
           OrgDb    = org.Hs.eg.db)      # Annotation mapping

dat <- tempOutput[, c("logFC", "adj.P.Val")]  # Keep effect size & FDR
dat <- dplyr::inner_join(dat, df, by = "SYMBOL")  # Merge with Entrez IDs
df2 <- dat[order(dat$logFC, decreasing = TRUE), ] # Sort by descending logFC

## Build named numeric vector required by clusterProfiler
genelist <- df2$logFC
names(genelist) <- df2$ENTREZID

################################################################################
## 6.  Gene-set enrichment analysis (GSEA) against KEGG pathways
################################################################################
KEGG <- gseKEGG(genelist, organism = "hsa")  # Perform preranked GSEA
KEGG_dat <- as.data.frame(KEGG)              # Convert result to data.frame

## IDs of four pathways of interest (immune-related)
id <- c("hsa04612", "hsa05169", "hsa05332", "hsa04650")

## Build two random colour gradients for aesthetics
c1 <- colorRampPalette(c("red", "yellow", "blue"), space = "rgb")(1000)
c2 <- colorRampPalette(c("#00994d", "#990099", "#990000"), space = "rgb")(1000)

################################################################################
## 7.  Generate GSEA running-score plots for the four chosen pathways
################################################################################
gse1 <- gseaNb(object   = KEGG,
               geneSetID = id[1],
               newGsea   = TRUE,
               rmHt      = TRUE,
               addPval   = TRUE,
               pvalX     = 0.9,
               pvalY     = 0.6,
               lineSize  = 1) +
  scale_color_gradient(low = c2[sample(1:1000, 1)],
                       high = c1[sample(1:1000, 1)]) +
  ylab("")

gse2 <- gseaNb(object   = KEGG,
               geneSetID = id[2],
               newGsea   = TRUE,
               rmHt      = TRUE,
               addPval   = TRUE,
               pvalX     = 0.9,
               pvalY     = 0.6,
               lineSize  = 1) +
  scale_color_gradient(low = c2[sample(1:1000, 1)],
                       high = c1[sample(1:1000, 1)]) +
  ylab("")

gse3 <- gseaNb(object   = KEGG,
               geneSetID = id[3],
               newGsea   = TRUE,
               rmHt      = TRUE,
               addPval   = TRUE,
               pvalX     = 0.9,
               pvalY     = 0.6,
               lineSize  = 1) +
  scale_color_gradient(low = c2[sample(1:1000, 1)],
                       high = c1[sample(1:1000, 1)]) +
  ylab("")

gse4 <- gseaNb(object   = KEGG,
               geneSetID = id[4],
               newGsea   = TRUE,
               rmHt      = TRUE,
               addPval   = TRUE,
               pvalX     = 0.9,
               pvalY     = 0.6,
               lineSize  = 1) +
  scale_color_gradient(low = c2[sample(1:1000, 1)],
                       high = c1[sample(1:1000, 1)]) +
  ylab("")