# README.md

## Differential Expression Analysis and Volcano Plot

This repository contains R scripts for performing differential gene expression analysis using the `limma` package and visualizing results with a volcano plot.

### Files:
- `differential_expression_analysis.R` - Performs differential expression analysis using `limma`.
- `volcano_plot.R` - Generates a volcano plot of differentially expressed genes.
- `requirements.R` - Installs necessary R packages.

### Usage:
1. Install required packages by running `source("requirements.R")` in R.
2. Execute `differential_expression_analysis.R` to generate differential expression results.
3. Run `volcano_plot.R` to create a volcano plot.

---

# requirements.R

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install(c('limma', 'EnhancedVolcano'))
install.packages(c('openxlsx', 'ggplot2', 'ggrepel'))

---

# differential_expression_analysis.R

library(Biobase)
library(limma)
library(openxlsx)

setwd("~/Documents/datasets/")
outDir <- "Results_DEGs"
dir.create(outDir, showWarnings = FALSE)

# Import metadata and count table
data <- read.xlsx("counts_brain.xlsx")
colnames(data) <- data[1, ]
data <- data[-1, ]
rownames(data) <- data[, 1]
data <- data[, -1]
metadata <- read.xlsx("AMY42_metadata.xlsx")

# Convert to numeric
my_data <- data.frame(sapply(data, as.numeric))
rownames(my_data) <- rownames(data)

design_df <- data.frame(Condition = factor(c(rep("MD", 13), rep("Control", 29))))
rownames(design_df) <- colnames(my_data)
design <- model.matrix(~ Condition, data = design_df)
colnames(design) <- c("Control", "MD")

fit <- lmFit(my_data, design)
contrast_matrix <- makeContrasts(Control - MD, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust = "BH", number = Inf)
write.xlsx(results, file.path(outDir, "differential_expression_results.xlsx"), rowNames = TRUE)

---

# volcano_plot.R

library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)

# Load results
deg_results <- read.xlsx("Results_DEGs/differential_expression_results.xlsx")

# Define colors
keyvals <- ifelse(deg_results$log2FoldChange > 1 & deg_results$P.Value < 0.05, '#23908F',
                   ifelse(deg_results$log2FoldChange < -1 & deg_results$P.Value < 0.05, '#FCE630', 'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#23908F'] <- 'up-regulated'
names(keyvals)[keyvals == 'gray'] <- 'nonDEGs'
names(keyvals)[keyvals == '#FCE630'] <- 'down-regulated'

# Plot volcano
p <- EnhancedVolcano(deg_results, lab = "", x = "log2FoldChange", y = "P.Value",
                     pCutoff = 0.05, FCcutoff = 1, xlim = c(-6, 6), ylim = c(0, 10),
                     xlab = "log2 Fold Change", ylab = "-log10 p-value", colCustom = keyvals)

# Save plots
ggsave("volcano_plot.tiff", plot = p, device = "tiff", width = 7, height = 5)
ggsave("volcano_plot.svg", plot = p, device = "svg", width = 6, height = 5)



