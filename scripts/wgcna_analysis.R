#!/usr/bin/env Rscript
# WGCNA on an expression matrix.
#
# Usage:
#   Rscript wgcna_analysis.R <expression.csv> <modules.csv> <eigengenes.csv> \
#                            <soft_threshold.csv> <power|NA> <min_module_size> \
#                            <merge_cut_height> <network_type>
#
# Inputs:
#   expression.csv   samples x genes, first column = sample IDs
#   power            integer or "NA" (auto-pick via pickSoftThreshold)
#   network_type     "unsigned" or "signed"
#
# Outputs:
#   modules.csv      Gene,Module
#   eigengenes.csv   SampleID,ME0,ME1,...
#   soft_threshold.csv  pickSoftThreshold fitIndices

suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript wgcna_analysis.R <input.csv> <modules.csv> <eigengenes.csv> <sft.csv> [power] [min_module_size] [merge_cut_height] [network_type]")
}

input_file       <- args[1]
modules_file     <- args[2]
eigengenes_file  <- args[3]
sft_file         <- args[4]
power_arg        <- if (length(args) >= 5) args[5] else "NA"
min_module_size  <- if (length(args) >= 6) as.integer(args[6]) else 10L
merge_cut_height <- if (length(args) >= 7) as.numeric(args[7]) else 0.25
network_type     <- if (length(args) >= 8) args[8] else "unsigned"

cat(sprintf("[wgcna_analysis.R] reading %s\n", input_file))
datExpr <- read.csv(input_file, row.names = 1, check.names = FALSE)
cat(sprintf("[wgcna_analysis.R] %d samples x %d genes\n",
            nrow(datExpr), ncol(datExpr)))

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  cat(sprintf("[wgcna_analysis.R] dropping %d genes / %d samples failing QC\n",
              sum(!gsg$goodGenes), sum(!gsg$goodSamples)))
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                         networkType = network_type, verbose = 0)
write.csv(sft$fitIndices, sft_file, row.names = FALSE, quote = FALSE)

if (toupper(power_arg) == "NA") {
  power <- sft$powerEstimate
  if (is.na(power)) {
    power <- if (network_type == "signed") 12L else 6L
    cat(sprintf("[wgcna_analysis.R] pickSoftThreshold returned NA, default beta=%d\n", power))
  } else {
    cat(sprintf("[wgcna_analysis.R] auto-selected beta=%d\n", power))
  }
} else {
  power <- as.integer(power_arg)
  cat(sprintf("[wgcna_analysis.R] using user-supplied beta=%d\n", power))
}

tom_type <- if (network_type == "signed") "signed" else "unsigned"
net <- blockwiseModules(
  datExpr,
  power            = power,
  networkType      = network_type,
  TOMType          = tom_type,
  minModuleSize    = min_module_size,
  mergeCutHeight   = merge_cut_height,
  numericLabels    = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs         = FALSE,
  verbose          = 0
)

modules <- net$colors
MEs     <- net$MEs

gene_names <- if (!is.null(names(modules))) names(modules) else colnames(datExpr)
module_df <- data.frame(
  Gene   = gene_names,
  Module = as.integer(modules)
)
write.csv(module_df, modules_file, row.names = FALSE, quote = FALSE)

eigengenes_out <- data.frame(SampleID = rownames(datExpr), MEs, check.names = FALSE)
write.csv(eigengenes_out, eigengenes_file, row.names = FALSE, quote = FALSE)

n_mod <- length(unique(modules))
cat(sprintf("[wgcna_analysis.R] PARAMS power=%d,modules=%d,network=%s,min_module_size=%d,merge_cut_height=%.3f\n",
            power, n_mod, network_type, min_module_size, merge_cut_height))
cat(sprintf("[wgcna_analysis.R] wrote %s, %s, %s\n",
            modules_file, eigengenes_file, sft_file))
