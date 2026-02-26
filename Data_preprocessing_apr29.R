# WGCNA preprocessing and module discovery workflow
#
# Usage:
#   Rscript Data_preprocessing_apr29.R <input_csv> [output_dir]
#
# input_csv: Gene expression matrix in CSV format.
#            Rows = genes, columns = samples.
#            First column is gene ID (used as row names).
# output_dir: Directory where outputs will be written (default: ./output)

suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
  library(dplyr)
  library(reshape2)
  library(tidyr)
  library(ggplot2)
})

allowWGCNAThreads()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Missing required argument: <input_csv>")
}

input_csv <- args[[1]]
output_dir <- ifelse(length(args) >= 2, args[[2]], "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading expression matrix: ", input_csv)
expr_raw <- read.csv(input_csv, row.names = 1, check.names = FALSE)

# Remove rows with missing values and keep only numeric columns.
expr_raw <- expr_raw[, vapply(expr_raw, is.numeric, logical(1)), drop = FALSE]
datExpr <- na.omit(expr_raw)

if (nrow(datExpr) < 100 || ncol(datExpr) < 4) {
  stop("Input data is too small for robust WGCNA analysis.")
}

# WGCNA expects samples as rows and genes as columns.
datExpr_t <- t(datExpr)

powers <- c(1:10, seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(datExpr_t, powerVector = powers, verbose = 5)

soft_power <- sft$powerEstimate
if (is.na(soft_power)) {
  warning("pickSoftThreshold did not return a power estimate; defaulting to 16")
  soft_power <- 16
}
message("Using soft-thresholding power: ", soft_power)

png(file.path(output_dir, "soft_threshold_diagnostics.png"), width = 1200, height = 500)
par(mfrow = c(1, 2))
cex1 <- 0.9

plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = "Scale independence"
)
text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red"
)
abline(h = 0.9, col = "red")

plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean connectivity"
)
text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1,
  col = "red"
)
dev.off()

# Avoid namespace conflicts with cor() used inside WGCNA.
temp_cor <- cor
cor <- WGCNA::cor
on.exit({
  cor <<- temp_cor
}, add = TRUE)

netwk <- blockwiseModules(
  datExpr_t,
  power = soft_power,
  networkType = "signed",
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minModuleSize = 30,
  maxBlockSize = 4000,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  saveTOMs = TRUE,
  saveTOMFileBase = file.path(output_dir, "TOM"),
  numericLabels = TRUE,
  verbose = 3
)

merged_colors <- labels2colors(netwk$colors)
module_df <- data.frame(
  gene_id = names(netwk$colors),
  module = merged_colors,
  stringsAsFactors = FALSE
)

write.csv(
  module_df,
  file = file.path(output_dir, "gene_modules.csv"),
  row.names = FALSE
)

# Eigengene matrix and heatmap
MEs0 <- moduleEigengenes(datExpr_t, merged_colors)$eigengenes
MEs0 <- orderMEs(MEs0)
MEs0$sample <- row.names(MEs0)

mME_base <- reshape2::melt(
  MEs0,
  id.vars = "sample",
  variable.name = "module",
  value.name = "eigengene"
)
mME_base$module <- gsub("ME", "", mME_base$module)

p <- ggplot(mME_base, aes(x = sample, y = module, fill = eigengene)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "Module Eigengene Heatmap",
    x = "Sample",
    y = "Modules",
    fill = "Eigengene"
  )

ggsave(
  filename = file.path(output_dir, "module_eigengene_heatmap.png"),
  plot = p,
  width = 12,
  height = 6,
  dpi = 300
)

message("Analysis complete. Outputs written to: ", normalizePath(output_dir))
