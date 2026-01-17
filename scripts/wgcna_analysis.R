#!/usr/bin/env Rscript
# WGCNA Analysis Script
# This script performs Weighted Gene Co-expression Network Analysis
# using the R WGCNA package

# Load required libraries
suppressPackageStartupMessages({
  library(WGCNA)
})

# Allow multi-threading
allowWGCNAThreads()

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript wgcna_analysis.R <input_csv> <module_output_csv> <eigengene_output_csv> [n_modules]")
}

input_file <- args[1]
module_output_file <- args[2]
eigengene_output_file <- args[3]
n_modules <- if (length(args) >= 4) as.integer(args[4]) else 4

cat("Performing WGCNA-like Co-expression Network Analysis\n\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# Read expression data
# Expected format: rows are samples, columns are genes
expr_data <- read.csv(input_file, row.names = 1, check.names = FALSE)

cat(sprintf("Loaded expression data: %d samples x %d genes\n", nrow(expr_data), ncol(expr_data)))

# Check for missing values
if (any(is.na(expr_data))) {
  cat("Warning: Found missing values, removing them\n")
  expr_data <- na.omit(expr_data)
}

# Check for genes with zero variance
good_genes <- apply(expr_data, 2, var) > 0
if (sum(!good_genes) > 0) {
  cat(sprintf("Removing %d genes with zero variance\n", sum(!good_genes)))
  expr_data <- expr_data[, good_genes]
}

# Calculate correlation matrix (genes x genes)
expr_corr <- cor(expr_data, method = "pearson", use = "pairwise.complete.obs")

# Convert correlation to distance
expr_dist <- 1 - abs(expr_corr)

cat(sprintf("Gene correlation matrix shape: (%d, %d)\n", nrow(expr_corr), ncol(expr_corr)))
cat(sprintf("Gene distance matrix shape: (%d, %d)\n", nrow(expr_dist), ncol(expr_dist)))

# Perform hierarchical clustering on genes
gene_tree <- hclust(as.dist(expr_dist), method = "average")

# Cut tree to get modules
module_labels <- cutree(gene_tree, k = n_modules)

# Adjust module labels to be 0-indexed (to match Python implementation)
module_labels <- module_labels - 1

# Create module assignment data frame
module_assignment <- data.frame(
  Gene = colnames(expr_data),
  Module = module_labels,
  stringsAsFactors = FALSE
)

cat(sprintf("\nIdentified %d co-expression modules\n", n_modules))
cat("Module sizes:\n")
module_counts <- table(module_assignment$Module)
for (i in sort(unique(module_assignment$Module))) {
  cat(sprintf("%d    %d\n", i, module_counts[as.character(i)]))
}

# Print genes in each module
for (mod in sort(unique(module_labels))) {
  genes_in_module <- module_assignment$Gene[module_assignment$Module == mod]
  cat(sprintf("Module %d: %s\n", mod, paste(genes_in_module, collapse=", ")))
}

# Calculate module eigengenes
# WGCNA moduleEigengenes function expects samples as rows, genes as columns
module_colors <- paste0("ME", module_labels)
names(module_colors) <- colnames(expr_data)

# Calculate eigengenes using WGCNA's moduleEigengenes function
MEs <- moduleEigengenes(expr_data, colors = module_colors)$eigengenes

# Rename columns to match expected format (ME0, ME1, etc.)
colnames(MEs) <- gsub("MEME", "ME", colnames(MEs))

# Reorder columns by module number
module_nums <- as.integer(gsub("ME", "", colnames(MEs)))
MEs <- MEs[, order(module_nums)]

cat("\n✓ Module eigengenes calculated\n")
cat(sprintf("  Shape: (%d, %d)\n", nrow(MEs), ncol(MEs)))

# Write results to CSV files
write.csv(module_assignment, module_output_file, row.names = FALSE, quote = FALSE)
cat(sprintf("\n✓ Module assignments saved to: %s\n", module_output_file))

# Add sample IDs as first column for eigengenes
eigengenes_output <- cbind(SampleID = rownames(expr_data), MEs)
write.csv(eigengenes_output, eigengene_output_file, row.names = FALSE, quote = FALSE)
cat(sprintf("✓ Module eigengenes saved to: %s\n", eigengene_output_file))

cat("\nWGCNA analysis complete!\n")
