## 0. Set working directory, load packages
library(limma)                                          # Differential-expression engine
library(ggplot2)                                        # Plotting
library(ggrepel)                                        # Non-overlapping labels
library(corrplot)                                       # Correlation matrix plots

## 1. Import data
rt   <- read.csv("ssGSEA_29.csv", header = TRUE, row.names = 1, check.names = FALSE)  # 29-pathway ssGSEA scores
risk <- readRDS("TCGA_group.rds")                       # Pre-computed risk groups (High / Low)

## 2. Split samples into normal and tumour based on TCGA barcode digits 14-15
nor <- rt[, as.numeric(substr(colnames(rt), 14, 15)) >= 10]  # Normal samples (NT)
tum <- rt[, as.numeric(substr(colnames(rt), 14, 15)) <  10]  # Tumour samples (TP)

## 3. Harmonise tumour sample names to 12-character TCGA ID for matching with risk object
colnames(tum) <- substr(colnames(tum), 1, 12)

## 4. Further split tumour samples by risk group
tum_high <- tum[, rownames(risk)[risk$Risk == "High"]]  # High-risk tumours
tum_low  <- tum[, rownames(risk)[risk$Risk == "Low"]]   # Low-risk tumours

###############################################################################
############### Part 1: High-risk tumours vs Normal (Volcano 1) ###############
###############################################################################
rt <- cbind(nor, tum_high)                              # Combine Normal + High-risk tumour matrix

# Build design matrix (no intercept)
group  <- factor(c(rep("Norml_group", ncol(nor)),
                   rep("High_risk_group", ncol(tum_high))))
design <- model.matrix(~0 + group)
colnames(design) <- c("High_risk_group", "Norml_group")

# limma differential-expression pipeline
fit      <- lmFit(rt, design)                           # Fit linear model
contrast <- makeContrasts(High_risk_group - Norml_group, levels = design)
fit      <- contrasts.fit(fit, contrast)                # Apply contrast
fit      <- eBayes(fit)                                 # Empirical Bayes moderation
tempOutput <- topTable(fit, n = Inf, adjust = "fdr")    # Retrieve full results with FDR
tempOutput$id <- rownames(tempOutput)                   # Add ID column

# Prepare data frame for plotting
df <- tempOutput
df$Group <- factor(
  ifelse(df$adj.P.Val < 0.05 & abs(df$logFC) >= 0.005,
         ifelse(df$logFC > 0.005, "Up", "Down"),
         "NotSignifi")
)

# Construct volcano plot 1
p1 <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val),
                     colour = Group, size = -log10(adj.P.Val))) +
  geom_rect(aes(xmin = -Inf, xmax = -0.005, ymin = -Inf, ymax = Inf),
            fill = "#ADD8E6", alpha = 0.2, color = NA) +
  geom_rect(aes(xmin = 0.005, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "#FFC0CB", alpha = 0.2, color = NA) +
  geom_point(shape = 20, stroke = 1) +
  scale_color_manual(values = c("steelblue", "#999999", "#CD5C5C")) +
  labs(x = "log2 (FoldChange)", y = "-log10 (adj.Pvalue)", fill = "", size = "") +
  geom_text_repel(
    data = subset(df, adj.P.Val < 0.05 & abs(logFC) > 0.005),
    aes(label = id), size = 3.5, color = "black", segment.color = "black", show.legend = FALSE
  ) +
  geom_vline(xintercept = -0.005, lty = 2, col = "#1E90FF", lwd = 1) +
  geom_vline(xintercept =  0.005, lty = 2, col = "#FF4500", lwd = 1) +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "#FF7F50", lwd = 1) +
  guides(fill = guide_legend(override.aes = list(size = 8))) +
  xlim(c(-0.08, 0.08)) +
  theme(
    axis.text        = element_text(face = "bold", colour = "#333333", size = 12, hjust = 0.5),
    axis.line        = element_blank(),
    plot.title       = element_text(face = "bold", colour = "#333333", size = 16),
    title            = element_text(face = "bold", colour = "#333333", size = 13),
    panel.border     = element_rect(fill = NA, color = "#333333", size = 1.2, linetype = "solid"),
    panel.background = element_rect(fill = "#F5F5F5"),
    panel.grid.major = element_line(color = "#D3D3D3", size = 0.5, linetype = "dotdash"),
    legend.position  = "right"
  )

###############################################################################
############### Part 2: Low-risk tumours vs Normal (Volcano 2) ################
###############################################################################
rt <- cbind(nor, tum_low)                               # Combine Normal + Low-risk tumour matrix

# Re-build design matrix (no intercept)
group  <- factor(c(rep("Norml_group", ncol(nor)),
                   rep("Low_risk_group", ncol(tum_low))))
design <- model.matrix(~0 + group)
colnames(design) <- c("Low_risk_group", "Norml_group")

# limma differential-expression pipeline (same as above)
fit      <- lmFit(rt, design)
contrast <- makeContrasts(Low_risk_group - Norml_group, levels = design)
fit      <- contrasts.fit(fit, contrast)
fit      <- eBayes(fit)
tempOutput <- topTable(fit, n = Inf, adjust = "fdr")
tempOutput$id <- rownames(tempOutput)

# Prepare data frame for plotting
df <- tempOutput
df$Group <- factor(
  ifelse(df$adj.P.Val < 0.05 & abs(df$logFC) >= 0.005,
         ifelse(df$logFC > 0.005, "Up", "Down"),
         "NotSignifi")
)

# Construct volcano plot 2 (identical formatting to p1)
p2 <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val),
                     colour = Group, size = -log10(adj.P.Val))) +
  geom_rect(aes(xmin = -Inf, xmax = -0.005, ymin = -Inf, ymax = Inf),
            fill = "#ADD8E6", alpha = 0.2, color = NA) +
  geom_rect(aes(xmin = 0.005, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "#FFC0CB", alpha = 0.2, color = NA) +
  geom_point(shape = 20, stroke = 1) +
  scale_color_manual(values = c("steelblue", "#999999", "#CD5C5C")) +
  labs(x = "log2 (FoldChange)", y = "-log10 (adj.Pvalue)", fill = "", size = "") +
  geom_text_repel(
    data = subset(df, adj.P.Val < 0.05 & abs(logFC) > 0.005),
    aes(label = id), size = 3.5, color = "black", segment.color = "black", show.legend = FALSE
  ) +
  geom_vline(xintercept = -0.005, lty = 2, col = "#1E90FF", lwd = 1) +
  geom_vline(xintercept =  0.005, lty = 2, col = "#FF4500", lwd = 1) +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "#FF7F50", lwd = 1) +
  guides(fill = guide_legend(override.aes = list(size = 8))) +
  xlim(c(-0.08, 0.08)) +
  theme(
    axis.text        = element_text(face = "bold", colour = "#333333", size = 12, hjust = 0.5),
    axis.line        = element_blank(),
    plot.title       = element_text(face = "bold", colour = "#333333", size = 16),
    title            = element_text(face = "bold", colour = "#333333", size = 13),
    panel.border     = element_rect(fill = NA, color = "#333333", size = 1.2, linetype = "solid"),
    panel.background = element_rect(fill = "#F5F5F5"),
    panel.grid.major = element_line(color = "#D3D3D3", size = 0.5, linetype = "dotdash"),
    legend.position  = "right"
  )

## Display the two volcano plots
p1
p2

###############################################################################
############### Part 3: Correlation matrices of High- vs Low-risk #############
###############################################################################
## Compute Pearson correlation matrices between pathways for each risk group
cor_group1 <- cor(t(tum_high))   # High-risk tumours (rows = pathways, cols = samples)
cor_group2 <- cor(t(tum_low))    # Low-risk tumours

## Define colour palettes
col1 <- colorRampPalette(c("purple", "white", "orange"))(200)  # Palette 1
col2 <- colorRampPalette(c("green",  "white", "red"))   (200)  # Palette 2

## Plot lower-triangle correlation for High-risk group
corrplot(cor_group1,
         method = "circle",
         type   = "lower",
         col    = col1,
         tl.col = "black",
         title  = "",
         tl.pos = "l")          # Label position = left

## Overlay upper-triangle correlation for Low-risk group
corrplot(cor_group2,
         method  = "circle",
         type    = "upper",
         col     = col2,
         add     = TRUE,        # Add to existing plot
         title   = "",
         tl.pos  = "t")         # Label position = top
