# Run DESeq2 Analysis using Docker (getwilds/deseq2)
import subprocess
import pandas as pd
import numpy as np
import os

print("Preparing data for DESeq2 Analysis\n")
print("="*80)

# Load the data files
try:
    expression_indexed = pd.read_csv('data/ASD_dataset/ASD_expression1.csv', index_col=0)
    # Load covariates
    covariates_indexed = pd.read_csv('data/ASD_dataset/ASD_covariates.csv', index_col=0)

    print(f"✓ Loaded expression data: {expression_indexed.shape}")
    print(f"✓ Loaded covariates data: {covariates_indexed.shape}")
except FileNotFoundError as e:
    print(f"ERROR: Could not load data files: {e}")
    print("Please ensure the ASD_dataset files exist.")
    exit(1)

# Note: DESeq2 works best with raw integer counts
# Create sample metadata
sample_metadata = pd.DataFrame({
    'Sample': covariates_indexed.index,
    'Condition': ['ASD' if x == 1 else 'Control' for x in covariates_indexed['ASD']],
    'Sex': covariates_indexed['Sex'],
    'Age': covariates_indexed['Age']
})

# Ensure the data directory exists
os.makedirs('data/ASD_dataset', exist_ok=True)

# Save to CSV for R to read
sample_metadata.to_csv('data/ASD_dataset/temp_sample_metadata_for_deseq2.csv', index=False)

print("✓ Data prepared for DESeq2")
print(f"  Count matrix shape: {expression_indexed.shape}")
print(f"  Sample metadata shape: {sample_metadata.shape}")

# R script for DESeq2
r_script = """
library(DESeq2)

# Read data
count_data <- read.csv('/data/ASD_expression1.csv', row.names=1)
sample_info <- read.csv('/data/temp_sample_metadata_for_deseq2.csv', row.names=1)

# Data is already in raw count format (integers)
# Transpose so genes are rows and samples are columns
count_data <- t(count_data)

# Ensure sample order matches
count_data <- count_data[, rownames(sample_info)]

# Verify we have valid count data
cat(sprintf("Count data range: %.0f to %.0f\\n", min(count_data), max(count_data)))
cat(sprintf("Mean count: %.1f\\n", mean(count_data)))
cat(sprintf("Any NA values: %s\\n", any(is.na(count_data))))

# Create DESeq2 dataset
# Using Condition as the main factor
sample_info$Condition <- factor(sample_info$Condition, levels = c("Control", "ASD"))

dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = sample_info,
    design = ~ Condition
)

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("Condition", "ASD", "Control"))

# Convert to data frame
res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)

# Reorder columns
res_df <- res_df[, c("Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

# Add significance column
# Using padj < 0.05 threshold (FDR-corrected p-value)
# No fold change threshold - let the statistics speak for themselves
res_df$Significant <- (!is.na(res_df$padj)) & (res_df$padj < 0.05)

# Sort by p-value
res_df <- res_df[order(res_df$pvalue), ]

# Save results
write.csv(res_df, '/data/Phase2_DESeq2_results.csv', row.names = FALSE)

# Print summary
cat("\\nDESeq2 Analysis Complete\\n")
cat("========================================\\n")
cat(sprintf("Total genes tested: %d\\n", nrow(res_df)))
cat(sprintf("Significant DEGs (padj < 0.05): %d\\n",
            sum(res_df$Significant, na.rm = TRUE)))
cat(sprintf("Upregulated in ASD: %d\\n",
            sum(res_df$Significant & res_df$log2FoldChange > 0, na.rm = TRUE)))
cat(sprintf("Downregulated in ASD: %d\\n",
            sum(res_df$Significant & res_df$log2FoldChange < 0, na.rm = TRUE)))

cat("\\nTop 10 Differentially Expressed Genes:\\n")
print(head(res_df[, c("Gene", "log2FoldChange", "pvalue", "padj", "Significant")], 10))
"""

# Save R script
with open('data/ASD_dataset/run_deseq2.R', 'w') as f:
    f.write(r_script)

print("✓ R script created")
print("\nRunning DESeq2 using Docker container...")
print("(This may take a few minutes on first run while pulling the image)\n")

# Get absolute path to data directory
data_dir = os.path.abspath('data/ASD_dataset')

# Run DESeq2 using Docker
try:
    # Docker command to run DESeq2
    # Mount the data directory to /data in the container
    docker_cmd = [
        'docker', 'run', '--rm',
        '-v', f'{data_dir}:/data',
        'getwilds/deseq2',
        'Rscript', '/data/run_deseq2.R'
    ]

    print(f"Running command: {' '.join(docker_cmd)}\n")

    # Run with real-time output
    result = subprocess.run(
        docker_cmd,
        capture_output=True,
        text=True,
        timeout=3600  # 60 minute timeout (includes Docker image pull time)
    )

    print(result.stdout)

    if result.returncode != 0:
        print("Docker/R Script Errors:")
        print(result.stderr)
        raise Exception("DESeq2 analysis failed")

    # Load results
    deseq2_results = pd.read_csv('data/ASD_dataset/Phase2_DESeq2_results.csv')
    
    if os.path.exists('data/ASD_dataset/temp_sample_metadata_for_deseq2.csv'):
        os.remove('data/ASD_dataset/temp_sample_metadata_for_deseq2.csv')

    print("\n" + "="*80)
    print("DESeq2 Results Summary:")
    print("="*80)
    print(deseq2_results[['Gene', 'log2FoldChange', 'pvalue', 'padj', 'Significant']].head(10).to_string(index=False))

    print("\n✓ DESeq2 analysis complete!")
    print("✓ Results saved to 'data/ASD_dataset/DESeq2_results.csv'")

    # Display some statistics
    print("\n" + "="*80)
    print("Final Statistics:")
    print("="*80)
    print(f"Total genes: {len(deseq2_results)}")
    print(f"Significant DEGs: {deseq2_results['Significant'].sum()}")
    print(f"Upregulated in ASD: {((deseq2_results['Significant']) & (deseq2_results['log2FoldChange'] > 0)).sum()}")
    print(f"Downregulated in ASD: {((deseq2_results['Significant']) & (deseq2_results['log2FoldChange'] < 0)).sum()}")

except FileNotFoundError:
    print("ERROR: Docker not found. Please ensure Docker is installed and running.")
    print("To install Docker: https://docs.docker.com/get-docker/")
    exit(1)
except subprocess.TimeoutExpired:
    print("ERROR: DESeq2 analysis timed out (>10 minutes)")
    exit(1)
except Exception as e:
    print(f"ERROR: {str(e)}")
    exit(1)
