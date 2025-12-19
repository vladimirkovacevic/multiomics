#!/usr/bin/env Rscript
#
# DIABLO Multi-Omics Integration Analysis
# Uses mixOmics package for supervised multi-block integration
#
# This script performs DIABLO analysis on expression and genotype data
# to identify multi-omics biomarker signatures that discriminate between
# ASD cases and controls.

# Suppress warnings for cleaner output
options(warn=-1)

# Load required libraries
suppressPackageStartupMessages({
  library(mixOmics)
  library(ggplot2)
})

cat("\n")
cat("================================================================================\n")
cat("  DIABLO Multi-Omics Integration Analysis\n")
cat("================================================================================\n")
cat("\n")

# ============================================================================
# 1. PARSE COMMAND LINE ARGUMENTS
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  cat("ERROR: Insufficient arguments\n\n")
  cat("Usage: Rscript diablo_analysis.R <expression_csv> <genotypes_csv> <covariates_csv> <output_prefix> [n_components] [keepX_expr] [keepX_geno]\n\n")
  cat("Arguments:\n")
  cat("  expression_csv  - Expression data (samples x genes)\n")
  cat("  genotypes_csv   - Genotype data (samples x SNPs)\n")
  cat("  covariates_csv  - Covariates with 'ASD' column (0=Control, 1=ASD)\n")
  cat("  output_prefix   - Prefix for output files\n")
  cat("  n_components    - Number of components (default: 2)\n")
  cat("  keepX_expr      - Genes to select per component (default: 5)\n")
  cat("  keepX_geno      - SNPs to select per component (default: 5)\n")
  cat("\n")
  quit(status = 1)
}

# Required arguments
expression_file <- args[1]
genotypes_file <- args[2]
covariates_file <- args[3]
output_prefix <- args[4]

# Optional arguments with defaults
n_components <- ifelse(length(args) >= 5, as.integer(args[5]), 2)
keepX_expr <- ifelse(length(args) >= 6, as.integer(args[6]), 5)
keepX_geno <- ifelse(length(args) >= 7, as.integer(args[7]), 5)

cat(sprintf("Configuration:\n"))
cat(sprintf("  Expression file:  %s\n", expression_file))
cat(sprintf("  Genotypes file:   %s\n", genotypes_file))
cat(sprintf("  Covariates file:  %s\n", covariates_file))
cat(sprintf("  Output prefix:    %s\n", output_prefix))
cat(sprintf("  Components:       %d\n", n_components))
cat(sprintf("  KeepX Expression: %d\n", keepX_expr))
cat(sprintf("  KeepX Genotypes:  %d\n", keepX_geno))
cat("\n")

# ============================================================================
# 2. LOAD AND VALIDATE DATA
# ============================================================================

cat("Loading data...\n")

# Check file existence
if (!file.exists(expression_file)) {
  stop(sprintf("Expression file not found: %s", expression_file))
}
if (!file.exists(genotypes_file)) {
  stop(sprintf("Genotypes file not found: %s", genotypes_file))
}
if (!file.exists(covariates_file)) {
  stop(sprintf("Covariates file not found: %s", covariates_file))
}

# Load data
expression <- read.csv(expression_file, row.names = 1, check.names = FALSE)
genotypes <- read.csv(genotypes_file, row.names = 1, check.names = FALSE)
covariates <- read.csv(covariates_file, row.names = 1, check.names = FALSE)

cat(sprintf("  Expression: %d samples × %d genes\n", nrow(expression), ncol(expression)))
cat(sprintf("  Genotypes:  %d samples × %d SNPs\n", nrow(genotypes), ncol(genotypes)))
cat(sprintf("  Covariates: %d samples × %d variables\n", nrow(covariates), ncol(covariates)))

# Verify sample alignment
if (!identical(rownames(expression), rownames(genotypes)) ||
    !identical(rownames(expression), rownames(covariates))) {
  stop("ERROR: Sample names do not match across datasets")
}

# Extract outcome variable
if (!"ASD" %in% colnames(covariates)) {
  stop("ERROR: 'ASD' column not found in covariates file")
}

Y <- factor(covariates$ASD, levels = c(0, 1), labels = c("Control", "ASD"))
cat(sprintf("  Outcome:    %d Controls, %d ASD cases\n", sum(Y == "Control"), sum(Y == "ASD")))
cat("\n")

# ============================================================================
# 3. PREPARE DATA FOR DIABLO
# ============================================================================

cat("Preparing data for DIABLO analysis...\n")

# Normalize expression data if needed (assume raw counts)
# Log2 transform if values are large (likely raw counts)
if (max(expression, na.rm = TRUE) > 100) {
  cat("  Applying log2(x + 1) transformation to expression data\n")
  expression <- log2(expression + 1)
}

# Center and scale both datasets
expression_scaled <- scale(expression, center = TRUE, scale = TRUE)
genotypes_scaled <- scale(genotypes, center = TRUE, scale = TRUE)

# Create multi-block data list
X <- list(
  Expression = expression_scaled,
  Genotypes = genotypes_scaled
)

cat("  Data prepared and scaled\n")
cat("\n")

# ============================================================================
# 4. DEFINE DESIGN MATRIX
# ============================================================================

cat("Setting up design matrix...\n")

# Design matrix: define correlation structure between blocks
# 0.1 = weak correlation (allows complementary information)
design <- matrix(0.1, ncol = 2, nrow = 2,
                dimnames = list(c("Expression", "Genotypes"),
                               c("Expression", "Genotypes")))
diag(design) <- 0  # Diagonal must be 0

cat("  Design matrix (block correlation):\n")
print(design)
cat("\n")

# ============================================================================
# 5. RUN INITIAL DIABLO MODEL (Full model)
# ============================================================================

cat("Running initial DIABLO model (no feature selection)...\n")

# Run full model to assess overall structure
diablo_full <- block.plsda(
  X = X,
  Y = Y,
  ncomp = n_components,
  design = design
)

cat(sprintf("  ✓ Initial model fitted with %d components\n", n_components))
cat("\n")

# ============================================================================
# 6. RUN SPARSE DIABLO MODEL (With feature selection)
# ============================================================================

cat("Running sparse DIABLO model with feature selection...\n")

# Define keepX for each block (number of features to select)
keepX_list <- list(
  Expression = rep(keepX_expr, n_components),
  Genotypes = rep(keepX_geno, n_components)
)

cat("  Feature selection per component:\n")
cat(sprintf("    Expression: %s\n", paste(keepX_list$Expression, collapse = ", ")))
cat(sprintf("    Genotypes:  %s\n", paste(keepX_list$Genotypes, collapse = ", ")))
cat("\n")

# Run sparse DIABLO
diablo_final <- block.splsda(
  X = X,
  Y = Y,
  ncomp = n_components,
  keepX = keepX_list,
  design = design
)

cat("  ✓ Sparse DIABLO model fitted\n")
cat("\n")

# ============================================================================
# 7. CROSS-VALIDATION FOR PERFORMANCE ASSESSMENT
# ============================================================================

cat("Performing cross-validation (this may take a moment)...\n")

# Perform 5-fold cross-validation with 10 repeats
set.seed(42)  # For reproducibility
perf_diablo <- perf(
  diablo_final,
  validation = "Mfold",
  folds = 5,
  nrepeat = 10,
  dist = "centroids.dist"
)

cat("  ✓ Cross-validation completed\n")
cat("\n")

# ============================================================================
# 8. EXTRACT RESULTS
# ============================================================================

cat("Extracting results...\n")

# -------------------------
# 8.1 Component Scores (Sample Projections)
# -------------------------
components_expr <- as.data.frame(diablo_final$variates$Expression)
components_geno <- as.data.frame(diablo_final$variates$Genotypes)

# Average components across blocks for consensus representation
components_avg <- (components_expr + components_geno) / 2
colnames(components_avg) <- paste0("Comp", 1:n_components)
components_avg$SampleID <- rownames(components_avg)
components_avg <- components_avg[, c("SampleID", paste0("Comp", 1:n_components))]

# -------------------------
# 8.2 Feature Loadings
# -------------------------
loadings_expr <- as.data.frame(diablo_final$loadings$Expression)
loadings_expr$Feature <- rownames(loadings_expr)
loadings_expr <- loadings_expr[, c("Feature", paste0("comp", 1:n_components))]
colnames(loadings_expr) <- c("Feature", paste0("Comp", 1:n_components))

loadings_geno <- as.data.frame(diablo_final$loadings$Genotypes)
loadings_geno$Feature <- rownames(loadings_geno)
loadings_geno <- loadings_geno[, c("Feature", paste0("comp", 1:n_components))]
colnames(loadings_geno) <- c("Feature", paste0("Comp", 1:n_components))

# -------------------------
# 8.3 Selected Features
# -------------------------
selected_features <- list()
for (comp in 1:n_components) {
  # Expression
  expr_vars <- selectVar(diablo_final, block = "Expression", comp = comp)
  if (!is.null(expr_vars$Expression)) {
    selected_features[[paste0("Expression_Comp", comp)]] <-
      data.frame(
        Component = comp,
        Block = "Expression",
        Feature = as.character(expr_vars$Expression$name),
        Loading = as.numeric(expr_vars$Expression$value$value.var),
        stringsAsFactors = FALSE
      )
  }

  # Genotypes
  geno_vars <- selectVar(diablo_final, block = "Genotypes", comp = comp)
  if (!is.null(geno_vars$Genotypes)) {
    selected_features[[paste0("Genotypes_Comp", comp)]] <-
      data.frame(
        Component = comp,
        Block = "Genotypes",
        Feature = as.character(geno_vars$Genotypes$name),
        Loading = as.numeric(geno_vars$Genotypes$value$value.var),
        stringsAsFactors = FALSE
      )
  }
}

selected_features_df <- do.call(rbind, selected_features)
rownames(selected_features_df) <- NULL

# -------------------------
# 8.4 Performance Metrics
# -------------------------
# Extract error rates from DIABLO perf object
# Use MajorityVote for multi-block classification performance
error_rates <- perf_diablo$MajorityVote.error.rate$centroids.dist

# Create performance dataframe
performance <- data.frame(
  Component = 1:n_components,
  Overall_Error = error_rates["Overall.ER", ],
  Control_Error = error_rates["Control", ],
  ASD_Error = error_rates["ASD", ],
  BER = error_rates["Overall.BER", ]
)

# -------------------------
# 8.5 Block Correlations
# -------------------------
# Compute correlation between blocks' components
correlations <- matrix(NA, nrow = n_components, ncol = 1)
for (i in 1:n_components) {
  correlations[i, 1] <- cor(components_expr[, i], components_geno[, i])
}

correlations_df <- data.frame(
  Component = 1:n_components,
  Expression_Genotypes_Cor = correlations[, 1]
)

cat("  ✓ Results extracted\n")
cat("\n")

# ============================================================================
# 9. SAVE RESULTS TO CSV FILES
# ============================================================================

cat("Saving results to CSV files...\n")

# Save all outputs with proper naming
write.csv(components_avg,
          file = paste0(output_prefix, "_components.csv"),
          row.names = FALSE)

write.csv(loadings_expr,
          file = paste0(output_prefix, "_loadings_expression.csv"),
          row.names = FALSE)

write.csv(loadings_geno,
          file = paste0(output_prefix, "_loadings_genotypes.csv"),
          row.names = FALSE)

write.csv(selected_features_df,
          file = paste0(output_prefix, "_selected_features.csv"),
          row.names = FALSE)

write.csv(performance,
          file = paste0(output_prefix, "_performance.csv"),
          row.names = FALSE)

write.csv(correlations_df,
          file = paste0(output_prefix, "_correlations.csv"),
          row.names = FALSE)

cat(sprintf("  ✓ %s_components.csv\n", output_prefix))
cat(sprintf("  ✓ %s_loadings_expression.csv\n", output_prefix))
cat(sprintf("  ✓ %s_loadings_genotypes.csv\n", output_prefix))
cat(sprintf("  ✓ %s_selected_features.csv\n", output_prefix))
cat(sprintf("  ✓ %s_performance.csv\n", output_prefix))
cat(sprintf("  ✓ %s_correlations.csv\n", output_prefix))
cat("\n")

# ============================================================================
# 10. PRINT SUMMARY STATISTICS
# ============================================================================

cat("================================================================================\n")
cat("  ANALYSIS SUMMARY\n")
cat("================================================================================\n")
cat("\n")

cat("Selected Features:\n")
cat(sprintf("  Expression: %d genes\n", sum(selected_features_df$Block == "Expression")))
cat(sprintf("  Genotypes:  %d SNPs\n", sum(selected_features_df$Block == "Genotypes")))
cat("\n")

cat("Classification Performance (Cross-Validation):\n")
for (i in 1:n_components) {
  cat(sprintf("  Component %d:\n", i))
  cat(sprintf("    Overall Error Rate: %.2f%%\n", performance$Overall_Error[i] * 100))
  cat(sprintf("    Balanced Error Rate: %.2f%%\n", performance$BER[i] * 100))
  cat(sprintf("    Control Error: %.2f%%\n", performance$Control_Error[i] * 100))
  cat(sprintf("    ASD Error: %.2f%%\n", performance$ASD_Error[i] * 100))
  cat(sprintf("    Block Correlation: %.3f\n", correlations_df$Expression_Genotypes_Cor[i]))
  cat("\n")
}

cat("Top Selected Features (Component 1):\n")
# Filter and sort expression features (base R)
top_expr <- selected_features_df[selected_features_df$Block == "Expression" & selected_features_df$Component == 1, ]
top_expr <- top_expr[order(-abs(top_expr$Loading)), ]
top_expr <- head(top_expr, 5)

# Filter and sort genotype features (base R)
top_geno <- selected_features_df[selected_features_df$Block == "Genotypes" & selected_features_df$Component == 1, ]
top_geno <- top_geno[order(-abs(top_geno$Loading)), ]
top_geno <- head(top_geno, 5)

cat("  Expression:\n")
for (i in 1:min(5, nrow(top_expr))) {
  cat(sprintf("    %s (loading: %.3f)\n", top_expr$Feature[i], top_expr$Loading[i]))
}

cat("  Genotypes:\n")
for (i in 1:min(5, nrow(top_geno))) {
  cat(sprintf("    %s (loading: %.3f)\n", top_geno$Feature[i], top_geno$Loading[i]))
}

cat("\n")
cat("================================================================================\n")
cat("  DIABLO ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("================================================================================\n")
cat("\n")

# Exit successfully
quit(status = 0)
