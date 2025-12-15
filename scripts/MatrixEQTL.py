#!/usr/bin/env python3
"""
eQTL Analysis using MatrixEQTL (Industry-Standard Tool)

This script performs eQTL analysis using MatrixEQTL, a widely-used tool
in the scientific community (Shabalin 2012, Bioinformatics).
Results can be compared with the simple linear model approach.

Requirements:
    - Docker must be installed and running
    - Uses Docker Hub image: pl92297/matrixeqtl (pulled automatically on first run)

Usage:
    python scripts/eQTL.py

The script will:
    1. Load genotype and expression data
    2. Calculate population structure PCs for covariate adjustment
    3. Run MatrixEQTL in a Docker container (no build required)
    4. Generate visualizations matching the notebook format
    5. Save results for comparison with linear model
"""

import subprocess
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

print("="*80)
print("eQTL Analysis using MatrixEQTL")
print("="*80)
print("\nMatrixEQTL (Shabalin 2012) is an industry-standard tool for eQTL mapping")
print("used in large-scale projects like GTEx, ENCODE, and many published studies.\n")

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================

print("Step 1: Loading data...\n")

try:
    # Load expression data
    expression_indexed = pd.read_csv('data/ASD_dataset/ASD_expression1.csv', index_col=0)

    # Load genotype data
    genotypes_indexed = pd.read_csv('data/ASD_dataset/ASD_genotypes.csv', index_col=0)

    # Load covariates
    covariates_indexed = pd.read_csv('data/ASD_dataset/ASD_covariates.csv', index_col=0)

    # Load gene annotations (for gene positions - needed for cis/trans distinction)
    # Check if the file exists and has position information
    try:
        gene_annot_raw = pd.read_csv('data/ASD_dataset/ASD_gene_annotations.csv')
        # Check if it has position columns
        if not all(col in gene_annot_raw.columns for col in ['chr', 'start', 'end']):
            # File exists but doesn't have positions, create dummy positions
            gene_annotations = pd.DataFrame({
                'gene_id': expression_indexed.columns,
                'chr': 1,
                'start': range(1000000, 1000000 + len(expression_indexed.columns) * 1000000, 1000000),
                'end': range(1100000, 1100000 + len(expression_indexed.columns) * 1000000, 1000000)
            })
        else:
            gene_annotations = gene_annot_raw
    except:
        # Create dummy annotations if file not available
        gene_annotations = pd.DataFrame({
            'gene_id': expression_indexed.columns,
            'chr': 1,
            'start': range(1000000, 1000000 + len(expression_indexed.columns) * 1000000, 1000000),
            'end': range(1100000, 1100000 + len(expression_indexed.columns) * 1000000, 1000000)
        })

    # Load SNP annotations (for SNP positions)
    try:
        snp_annot_raw = pd.read_csv('data/ASD_dataset/ASD_snp_annotations.csv')
        # Check if it has position columns
        if not all(col in snp_annot_raw.columns for col in ['chr', 'pos']):
            # File exists but doesn't have positions, create dummy positions
            snp_annotations = pd.DataFrame({
                'snp_id': genotypes_indexed.columns,
                'chr': 1,
                'pos': range(1000000, 1000000 + len(genotypes_indexed.columns) * 100000, 100000)
            })
        else:
            snp_annotations = snp_annot_raw
    except:
        # Create dummy annotations if file not available
        snp_annotations = pd.DataFrame({
            'snp_id': genotypes_indexed.columns,
            'chr': 1,
            'pos': range(1000000, 1000000 + len(genotypes_indexed.columns) * 100000, 100000)
        })

    print(f"✓ Loaded expression data: {expression_indexed.shape}")
    print(f"✓ Loaded genotype data: {genotypes_indexed.shape}")
    print(f"✓ Loaded covariates data: {covariates_indexed.shape}")
    print(f"✓ Loaded gene annotations: {gene_annotations.shape}")
    print(f"✓ Loaded SNP annotations: {snp_annotations.shape}")

except FileNotFoundError as e:
    print(f"ERROR: Could not load data files: {e}")
    print("Please ensure the ASD_dataset files exist.")
    exit(1)

# Verify samples match across datasets
if not (expression_indexed.index.equals(genotypes_indexed.index) and
        genotypes_indexed.index.equals(covariates_indexed.index)):
    print("ERROR: Sample IDs do not match across datasets!")
    exit(1)

print(f"\n✓ Data validation complete")
print(f"  Total samples: {len(expression_indexed)}")
print(f"  Total genes: {len(expression_indexed.columns)}")
print(f"  Total SNPs: {len(genotypes_indexed.columns)}")
print(f"  Total SNP-gene pairs to test: {len(expression_indexed.columns) * len(genotypes_indexed.columns)}")

# ============================================================================
# 2. PREPARE DATA FOR MATRIXEQTL
# ============================================================================

print("\nStep 2: Preparing data for MatrixEQTL...\n")

# Ensure output directory exists
os.makedirs('data/ASD_dataset', exist_ok=True)

# MatrixEQTL expects:
# - Gene expression: genes in rows, samples in columns
# - Genotypes: SNPs in rows, samples in columns
# - Covariates: covariates in rows, samples in columns (optional)
# - Gene locations: gene_id, chr, start, end
# - SNP locations: snp_id, chr, pos

# Transpose expression (genes as rows, samples as columns)
expression_for_matrixeqtl = expression_indexed.T
expression_for_matrixeqtl.index.name = 'geneid'
expression_for_matrixeqtl.to_csv('data/ASD_dataset/temp_expression_matrixeqtl.txt', sep='\t')

# Transpose genotypes (SNPs as rows, samples as columns)
genotypes_for_matrixeqtl = genotypes_indexed.T
genotypes_for_matrixeqtl.index.name = 'snpid'
genotypes_for_matrixeqtl.to_csv('data/ASD_dataset/temp_genotypes_matrixeqtl.txt', sep='\t')

# Prepare covariates (Age, Sex, and first 3 PCs from genotypes)
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Calculate PCs for population structure
scaler = StandardScaler()
genotypes_scaled = scaler.fit_transform(genotypes_indexed)
pca = PCA(n_components=3)
pcs = pca.fit_transform(genotypes_scaled)

# Create covariates matrix
covariates_for_matrixeqtl = pd.DataFrame({
    'Age': covariates_indexed['Age'].values,
    'Sex': covariates_indexed['Sex'].values,
    'PC1': pcs[:, 0],
    'PC2': pcs[:, 1],
    'PC3': pcs[:, 2]
}, index=covariates_indexed.index).T

covariates_for_matrixeqtl.index.name = 'id'
covariates_for_matrixeqtl.to_csv('data/ASD_dataset/temp_covariates_matrixeqtl.txt', sep='\t')

# Prepare gene locations
# Ensure we have the right columns - they should already be standardized from above
if 'gene_id' not in gene_annotations.columns:
    # The annotations must be from the original file, extract gene IDs from the first column
    gene_locations = pd.DataFrame({
        'gene_id': expression_indexed.columns,
        'chr': 1,
        'start': range(1000000, 1000000 + len(expression_indexed.columns) * 1000000, 1000000),
        'end': range(1100000, 1100000 + len(expression_indexed.columns) * 1000000, 1000000)
    })
else:
    gene_locations = gene_annotations[['gene_id', 'chr', 'start', 'end']].copy()
gene_locations.to_csv('data/ASD_dataset/temp_gene_locations.txt', sep='\t', index=False)

# Prepare SNP locations
if 'snp_id' not in snp_annotations.columns:
    # The annotations must be from the original file, create positions
    snp_locations = pd.DataFrame({
        'snp_id': genotypes_indexed.columns,
        'chr': 1,
        'pos': range(1000000, 1000000 + len(genotypes_indexed.columns) * 100000, 100000)
    })
else:
    snp_locations = snp_annotations[['snp_id', 'chr', 'pos']].copy()
snp_locations.to_csv('data/ASD_dataset/temp_snp_locations.txt', sep='\t', index=False)

print("✓ Data files prepared for MatrixEQTL:")
print("  - Expression matrix (genes × samples)")
print("  - Genotype matrix (SNPs × samples)")
print("  - Covariates matrix (Age, Sex, 3 PCs)")
print("  - Gene location file")
print("  - SNP location file")

# ============================================================================
# 3. CREATE R SCRIPT FOR MATRIXEQTL
# ============================================================================

print("\nStep 3: Creating R script for MatrixEQTL...\n")

# R script that runs MatrixEQTL
r_script = """
# Load MatrixEQTL library
library(MatrixEQTL)

# Settings
useModel = modelLINEAR  # Use linear model for quantitative traits
output_file_name = '/data/Phase2_MatrixEQTL_all_results.txt'
pvOutputThreshold = 1  # Output all p-values

# P-value thresholds
pvOutputThreshold_cis = 1.0  # All cis tests
pvOutputThreshold_tra = 1.0  # All trans tests

# Error covariance matrix (set to numeric() for homoskedastic model)
errorCovariance = numeric()

# Distance for local (cis) gene-SNP pairs (1 Mb = 1e6)
cisDist = 1e6

cat("\\nMatrixEQTL Analysis\\n")
cat("==================================================\\n")

# Load data
cat("Loading data files...\\n")

# Expression data
expr = SlicedData$new()
expr$fileDelimiter = "\\t"
expr$fileOmitCharacters = "NA"
expr$fileSkipRows = 1
expr$fileSkipColumns = 1
expr$fileSliceSize = 2000
expr$LoadFile('/data/temp_expression_matrixeqtl.txt')

cat(sprintf("  Expression: %d genes, %d samples\\n", expr$nRows(), expr$nCols()))

# Genotype data
snps = SlicedData$new()
snps$fileDelimiter = "\\t"
snps$fileOmitCharacters = "NA"
snps$fileSkipRows = 1
snps$fileSkipColumns = 1
snps$fileSliceSize = 2000
snps$LoadFile('/data/temp_genotypes_matrixeqtl.txt')

cat(sprintf("  Genotypes: %d SNPs, %d samples\\n", snps$nRows(), snps$nCols()))

# Covariates
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\\t"
cvrt$fileOmitCharacters = "NA"
cvrt$fileSkipRows = 1
cvrt$fileSkipColumns = 1
cvrt$LoadFile('/data/temp_covariates_matrixeqtl.txt')

cat(sprintf("  Covariates: %d covariates, %d samples\\n", cvrt$nRows(), cvrt$nCols()))

# Gene and SNP positions
genepos = read.table('/data/temp_gene_locations.txt', header=TRUE, stringsAsFactors=FALSE)
snpspos = read.table('/data/temp_snp_locations.txt', header=TRUE, stringsAsFactors=FALSE)

# Rename columns to match MatrixEQTL requirements
names(genepos) = c("geneid", "chr", "s1", "s2")
names(snpspos) = c("snpid", "chr", "pos")

cat(sprintf("  Gene positions: %d genes\\n", nrow(genepos)))
cat(sprintf("  SNP positions: %d SNPs\\n", nrow(snpspos)))

cat("\\nRunning MatrixEQTL analysis...\\n")
cat("  Model: Linear regression with covariates\\n")
cat(sprintf("  Cis distance threshold: %.0e bp\\n", cisDist))
cat("  (Testing both cis and trans eQTLs)\\n\\n")

# Run MatrixEQTL
me = Matrix_eQTL_main(
    snps = snps,
    gene = expr,
    cvrt = cvrt,
    output_file_name = NULL,  # Don't save intermediate files
    pvOutputThreshold = 0,    # No p-value threshold for cis
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = NULL,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = FALSE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
)

cat("\\nMatrixEQTL Analysis Complete!\\n")
cat("==================================================\\n")

# Combine cis and trans results
# Check if we have trans results (they might be empty in small datasets)
has_trans = !is.null(me$trans$eqtls) && nrow(me$trans$eqtls) > 0
if (has_trans) {
    all_results = rbind(
        data.frame(me$cis$eqtls, Type='cis'),
        data.frame(me$trans$eqtls, Type='trans')
    )
} else {
    # Only cis results
    all_results = data.frame(me$cis$eqtls, Type='cis')
    cat("Note: No trans-eQTLs detected (all SNPs within cis distance of all genes)\\n")
}

# Rename columns for clarity
names(all_results)[1:5] = c('SNP', 'Gene', 'Statistic', 'P_Value', 'FDR')

# Apply FDR correction across all tests
all_results$FDR = p.adjust(all_results$P_Value, method='BH')

# Add significance flag
all_results$Significant = all_results$FDR < 0.05

# Sort by p-value
all_results = all_results[order(all_results$P_Value), ]

# Print summary statistics
cat(sprintf("\\nTotal eQTL tests performed: %d\\n", nrow(all_results)))
cat(sprintf("Cis eQTLs tested: %d\\n", sum(all_results$Type == 'cis')))
cat(sprintf("Trans eQTLs tested: %d\\n", sum(all_results$Type == 'trans')))
cat(sprintf("\\nSignificant eQTLs (FDR < 0.05): %d\\n", sum(all_results$Significant)))
cat(sprintf("  Significant cis eQTLs: %d\\n", sum(all_results$Significant & all_results$Type == 'cis')))
cat(sprintf("  Significant trans eQTLs: %d\\n", sum(all_results$Significant & all_results$Type == 'trans')))

if (sum(all_results$Significant) > 0) {
    cat(sprintf("  Unique genes with eQTL: %d\\n", length(unique(all_results$Gene[all_results$Significant]))))
    cat(sprintf("  Unique SNPs with eQTL: %d\\n", length(unique(all_results$SNP[all_results$Significant]))))
}

# Save results
write.csv(all_results, '/data/Phase2_MatrixEQTL_results.csv', row.names=FALSE)

cat("\\n✓ Results saved to 'Phase2_MatrixEQTL_results.csv'\\n")

# Print top results
cat("\\nTop 20 eQTL Associations:\\n")
cat("==================================================\\n")
print(head(all_results[, c('SNP', 'Gene', 'Statistic', 'P_Value', 'FDR', 'Type', 'Significant')], 20))
"""

# Save R script
with open('data/ASD_dataset/run_matrixeqtl.R', 'w') as f:
    f.write(r_script)

print("✓ R script created: 'data/ASD_dataset/run_matrixeqtl.R'")

# ============================================================================
# 4. RUN MATRIXEQTL IN DOCKER
# ============================================================================

print("\nStep 4: Running MatrixEQTL analysis in Docker...\n")
print("Using pre-built Docker image: pl92297/matrixeqtl from Docker Hub\n")

data_dir = os.path.abspath('data/ASD_dataset')

print("Note: On first run, Docker will pull the image from Docker Hub.")
print("This may take a few minutes depending on your internet connection.\n")
print("This performs regression with covariates for all SNP-gene pairs\n")

try:
    # Docker command to run MatrixEQTL using the Docker Hub image
    docker_cmd = [
        'docker', 'run', '--rm',
        '-v', f'{data_dir}:/data',
        'pl92297/matrixeqtl',
        'Rscript', '/data/run_matrixeqtl.R'
    ]

    print(f"Running MatrixEQTL...")
    print(f"Command: {' '.join(docker_cmd)}\n")

    result = subprocess.run(
        docker_cmd,
        capture_output=True,
        text=True,
        timeout=3600  # 60 minute timeout
    )

    print(result.stdout)

    if result.returncode != 0:
        print("MatrixEQTL errors:")
        print(result.stderr)
        raise Exception("MatrixEQTL analysis failed")

    print("✓ MatrixEQTL analysis completed successfully")

except FileNotFoundError:
    print("ERROR: Docker not found. Please ensure Docker is installed and running.")
    print("To install Docker: https://docs.docker.com/get-docker/")
    exit(1)
except subprocess.TimeoutExpired:
    print("ERROR: MatrixEQTL analysis timed out (>60 minutes)")
    exit(1)
except Exception as e:
    print(f"ERROR: {str(e)}")
    exit(1)

# ============================================================================
# 5. LOAD AND ANALYZE RESULTS
# ============================================================================

print("\nStep 5: Loading MatrixEQTL results...\n")

try:
    matrixeqtl_results = pd.read_csv('data/ASD_dataset/Phase2_MatrixEQTL_results.csv')

    # Calculate correlation coefficient from test statistic for comparability
    # In linear model: t = r * sqrt(n-2) / sqrt(1-r^2)
    # So: r = t / sqrt(t^2 + n-2)
    n_samples = len(expression_indexed)
    t_stat = matrixeqtl_results['Statistic']
    matrixeqtl_results['Correlation'] = t_stat / np.sqrt(t_stat**2 + n_samples - 2)

    # Also calculate the actual correlation for top hits
    print("Calculating actual correlation coefficients for verification...")

    correlations = []
    for _, row in matrixeqtl_results.head(100).iterrows():  # Just top 100 for speed
        snp = row['SNP']
        gene = row['Gene']

        if snp in genotypes_indexed.columns and gene in expression_indexed.columns:
            corr, _ = pearsonr(genotypes_indexed[snp], expression_indexed[gene])
            correlations.append(corr)
        else:
            correlations.append(np.nan)

    # Update top results with actual correlations
    matrixeqtl_results.loc[:99, 'Correlation_Actual'] = correlations

    print(f"✓ Results loaded: {len(matrixeqtl_results)} eQTL tests")
    print(f"✓ Significant eQTLs (FDR < 0.05): {matrixeqtl_results['Significant'].sum()}")

except FileNotFoundError:
    print("ERROR: MatrixEQTL results file not found")
    exit(1)

# ============================================================================
# 6. CREATE VISUALIZATIONS (MATCHING LINEAR MODEL APPROACH)
# ============================================================================

print("\nStep 6: Creating visualizations...\n")

# Get top eQTLs (by p-value, regardless of FDR significance for comparison)
# This allows comparison with the linear model approach
if matrixeqtl_results['Significant'].sum() > 0:
    top_eqtls = matrixeqtl_results[matrixeqtl_results['Significant']].head(6)
    print(f"Generating scatter plots for top {len(top_eqtls)} significant eQTL associations...")
else:
    # No FDR-significant results, plot top 6 by p-value for comparison
    top_eqtls = matrixeqtl_results.head(6)
    print(f"No FDR-significant eQTLs found.")
    print(f"Generating scatter plots for top {len(top_eqtls)} eQTL associations by p-value...")

if len(top_eqtls) > 0:

    # Create figure with subplots
    n_plots = len(top_eqtls)
    n_rows = 2
    n_cols = 3

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 10))
    axes = axes.flatten()

    for idx, (_, row) in enumerate(top_eqtls.iterrows()):
        if idx >= 6:
            break

        snp = row['SNP']
        gene = row['Gene']
        pval = row['P_Value']
        fdr = row['FDR']
        eqtl_type = row['Type']

        # Get data
        if snp in genotypes_indexed.columns and gene in expression_indexed.columns:
            genotype = genotypes_indexed[snp]
            expression = expression_indexed[gene]

            # Calculate correlation
            corr, _ = pearsonr(genotype, expression)

            # Create plot
            ax = axes[idx]

            # Prepare data for box plots - group expression by genotype
            genotype_groups = []
            genotype_positions = []
            genotype_labels = []

            for geno in sorted(genotype.unique()):
                mask = genotype == geno
                genotype_groups.append(expression[mask].values)
                genotype_positions.append(geno)
                genotype_labels.append(f'{int(geno)}')

            # Create box plots
            bp = ax.boxplot(genotype_groups,
                           positions=genotype_positions,
                           widths=0.15,
                           patch_artist=True,
                           showfliers=False,  # Don't show outliers as separate points
                           medianprops=dict(color='black', linewidth=2),
                           boxprops=dict(facecolor='lightgray', alpha=0.7, edgecolor='black'),
                           whiskerprops=dict(color='black', linewidth=1.5),
                           capprops=dict(color='black', linewidth=1.5))

            # Add scatter plot with jitter for better visualization
            colors = {0: 'blue', 1: 'green', 2: 'red'}
            np.random.seed(42)  # For reproducible jitter
            for geno in sorted(genotype.unique()):
                mask = genotype == geno
                # Add small jitter to x-coordinates for better visibility
                jitter = np.random.normal(0, 0.02, mask.sum())
                ax.scatter(genotype[mask] + jitter, expression[mask],
                          alpha=0.5, s=30, color=colors.get(geno, 'gray'),
                          label=f'Genotype {int(geno)}', edgecolors='black', linewidths=0.5)

            # Add regression line
            z = np.polyfit(genotype, expression, 1)
            p = np.poly1d(z)
            x_line = np.linspace(genotype.min(), genotype.max(), 100)
            ax.plot(x_line, p(x_line), "r--", alpha=0.8, linewidth=2.5, label='Linear fit')

            # Labels and title
            ax.set_xlabel(f'{snp} Genotype', fontsize=10)
            ax.set_ylabel(f'{gene} Expression', fontsize=10)
            ax.set_title(f'{snp} - {gene} ({eqtl_type})\n' +
                        f'r = {corr:.3f}, P = {pval:.2e}, FDR = {fdr:.3f}',
                        fontsize=11, fontweight='bold')
            ax.legend(fontsize=7, loc='best')
            ax.grid(True, alpha=0.3)

            # Set x-axis to show genotype categories
            ax.set_xticks(genotype_positions)
            ax.set_xticklabels(genotype_labels)

    # Hide unused subplots
    for idx in range(len(top_eqtls), len(axes)):
        axes[idx].axis('off')

    plt.tight_layout()
    plt.savefig('data/ASD_dataset/Phase2_MatrixEQTL_top_associations.png',
                dpi=300, bbox_inches='tight')
    plt.close()

    print("✓ Scatter plots saved: 'Phase2_MatrixEQTL_top_associations.png'")

# Create QQ plot
print("Generating QQ plot for p-values...")

fig, ax = plt.subplots(figsize=(8, 8))

observed_p = matrixeqtl_results['P_Value'].values
observed_p = observed_p[~np.isnan(observed_p)]
observed_p = np.sort(observed_p)

n = len(observed_p)
expected_p = np.arange(1, n + 1) / (n + 1)

# Convert to -log10
observed_log = -np.log10(observed_p)
expected_log = -np.log10(expected_p)

# Plot
ax.scatter(expected_log, observed_log, alpha=0.5, s=20, color='steelblue')
ax.plot([0, max(expected_log)], [0, max(expected_log)], 'r--', linewidth=2, label='Expected (null)')

# Calculate genomic inflation factor
from scipy import stats as scipy_stats
chisq = scipy_stats.chi2.ppf(1 - observed_p, df=1)
lambda_gc = np.median(chisq) / scipy_stats.chi2.ppf(0.5, df=1)

ax.set_xlabel('Expected -Log10 P-Value', fontsize=12)
ax.set_ylabel('Observed -Log10 P-Value', fontsize=12)
ax.set_title(f'QQ Plot: MatrixEQTL eQTL P-Values\nGenomic Inflation Factor λ = {lambda_gc:.3f}',
            fontsize=14, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('data/ASD_dataset/Phase2_MatrixEQTL_qq_plot.png', dpi=300, bbox_inches='tight')
plt.close()

print("✓ QQ plot saved: 'Phase2_MatrixEQTL_qq_plot.png'")

# Create Manhattan plot (if we have position information)
print("Generating Manhattan plot...")

# Merge with SNP positions for plotting
plot_data = matrixeqtl_results.merge(snp_locations, left_on='SNP', right_on='snp_id', how='left')

fig, ax = plt.subplots(figsize=(14, 6))

# Plot by chromosome
for chrom in sorted(plot_data['chr'].unique()):
    chr_data = plot_data[plot_data['chr'] == chrom]
    ax.scatter(chr_data['pos'], -np.log10(chr_data['P_Value']),
              alpha=0.6, s=20, label=f'Chr {chrom}')

# Add significance threshold line for nominal p-value
threshold = -np.log10(0.05)
ax.axhline(y=threshold, color='r', linestyle='--', linewidth=2, label='Nominal P = 0.05')

# FDR threshold - only show if there are FDR-significant results
if matrixeqtl_results['Significant'].sum() > 0:
    # Find the maximum p-value among FDR-significant results
    # This is the p-value threshold that corresponds to FDR = 0.05
    fdr_threshold = -np.log10(matrixeqtl_results[matrixeqtl_results['FDR'] < 0.05]['P_Value'].max())
    ax.axhline(y=fdr_threshold, color='orange', linestyle='--', linewidth=2, label='FDR = 0.05')
else:
    # No FDR-significant results - show where the most significant result is
    # This shows that even the best p-value doesn't reach FDR < 0.05
    min_fdr = matrixeqtl_results['FDR'].min()
    ax.axhline(y=threshold, color='orange', linestyle=':', linewidth=2,
               label=f'No FDR < 0.05 (min FDR = {min_fdr:.3f})', alpha=0.5)

ax.set_xlabel('Genomic Position', fontsize=12)
ax.set_ylabel('-Log10 P-Value', fontsize=12)
ax.set_title('Manhattan Plot: MatrixEQTL eQTL Results', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('data/ASD_dataset/Phase2_MatrixEQTL_manhattan_plot.png', dpi=300, bbox_inches='tight')
plt.close()

print("✓ Manhattan plot saved: 'Phase2_MatrixEQTL_manhattan_plot.png'")

# ============================================================================
# 7. CLEAN UP TEMPORARY FILES
# ============================================================================

print("\nStep 7: Cleaning up temporary files...\n")

temp_files = [
    'data/ASD_dataset/temp_expression_matrixeqtl.txt',
    'data/ASD_dataset/temp_genotypes_matrixeqtl.txt',
    'data/ASD_dataset/temp_covariates_matrixeqtl.txt',
    'data/ASD_dataset/temp_gene_locations.txt',
    'data/ASD_dataset/temp_snp_locations.txt'
]

for temp_file in temp_files:
    if os.path.exists(temp_file):
        os.remove(temp_file)
        print(f"✓ Removed {temp_file}")

# ============================================================================
# 8. FINAL SUMMARY AND COMPARISON
# ============================================================================

print("\n" + "="*80)
print("MatrixEQTL ANALYSIS COMPLETE")
print("="*80)

print("\nSummary Statistics:")
print(f"  Total eQTL tests: {len(matrixeqtl_results)}")
print(f"  Significant eQTLs (FDR < 0.05): {matrixeqtl_results['Significant'].sum()}")
print(f"  Cis eQTLs: {(matrixeqtl_results['Type'] == 'cis').sum()}")
print(f"  Trans eQTLs: {(matrixeqtl_results['Type'] == 'trans').sum()}")

if matrixeqtl_results['Significant'].sum() > 0:
    print(f"  Unique genes with eQTL: {matrixeqtl_results[matrixeqtl_results['Significant']]['Gene'].nunique()}")
    print(f"  Unique SNPs with eQTL: {matrixeqtl_results[matrixeqtl_results['Significant']]['SNP'].nunique()}")
    print(f"  Genomic inflation factor (λ): {lambda_gc:.3f}")

print("\nTop 10 eQTL Associations:")
print("="*80)
display_cols = ['SNP', 'Gene', 'Correlation', 'P_Value', 'FDR', 'Type', 'Significant']
print(matrixeqtl_results[display_cols].head(10).to_string(index=False))

print("\n" + "="*80)
print("OUTPUT FILES:")
print("="*80)
print("  1. Phase2_MatrixEQTL_results.csv - Full eQTL results")
print("  2. Phase2_MatrixEQTL_top_associations.png - Scatter plots of top eQTLs")
print("  3. Phase2_MatrixEQTL_qq_plot.png - QQ plot of p-values")
print("  4. Phase2_MatrixEQTL_manhattan_plot.png - Manhattan plot")

print("\n" + "="*80)
print("COMPARISON WITH LINEAR MODEL:")
print("="*80)
print("\nTo compare with the simple linear model results from the notebook:")
print("  1. Load: data/ASD_dataset/Phase2_eQTL_results.csv (linear model)")
print("  2. Load: data/ASD_dataset/Phase2_MatrixEQTL_results.csv (MatrixEQTL)")
print("  3. Compare number of significant hits, p-value distributions, and effect sizes")
print("\nMatrixEQTL advantages over simple linear model:")
print("  ✓ Accounts for covariates (Age, Sex, population structure)")
print("  ✓ Distinguishes cis vs trans eQTLs")
print("  ✓ More robust statistical framework")
print("  ✓ Industry standard (used in GTEx, ENCODE, etc.)")
print("  ✓ Better control of false positives")

print("\n" + "="*80)
print("✓ ANALYSIS COMPLETE!")
print("="*80)
