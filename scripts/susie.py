#!/usr/bin/env python3
"""
SuSiE Fine-Mapping Analysis

This script performs SuSiE (Sum of Single Effects) fine-mapping on eQTL
associations to identify credible sets of likely causal variants.
SuSiE complements MatrixEQTL by narrowing broad association signals
to small sets of variants most likely to be causal.

Requirements:
    - Docker must be installed and running
    - Uses Docker Hub image: eqtlcatalogue/susie-finemapping

Usage:
    python scripts/susie.py

The script will:
    1. Load genotype and expression data
    2. Regress out covariates to get residualized expression
    3. Run SuSiE fine-mapping in a Docker container
    4. Generate visualizations (PIP plots, credible sets, comparisons)
    5. Save results for comparison with MatrixEQTL
"""

import subprocess
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression

print("=" * 80)
print("SuSiE Fine-Mapping Analysis")
print("=" * 80)
print("\nSuSiE (Wang et al. 2020, JRSS-B) identifies credible sets of causal variants")
print("from eQTL associations, resolving LD to pinpoint likely causal SNPs.\n")

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================

print("Step 1: Loading data...\n")

try:
    # Load raw expression counts
    expression_raw = pd.read_csv('data/ASD_dataset/ASD_expression.csv', index_col=0)
    genotypes = pd.read_csv('data/ASD_dataset/ASD_genotypes.csv', index_col=0)
    covariates = pd.read_csv('data/ASD_dataset/ASD_covariates.csv', index_col=0)

    # Inverse Normal Transformation (INT) per gene — industry standard (GTEx, eQTLGen)
    from scipy.stats import rankdata
    from scipy.stats import norm as sp_norm

    def inverse_normal_transform(x):
        """Rank-based inverse normal transformation (Blom's method)."""
        n = len(x)
        ranks = rankdata(x)
        quantiles = (ranks - 3/8) / (n + 1/4)
        return sp_norm.ppf(quantiles)

    expression = expression_raw.apply(inverse_normal_transform, axis=0)
    print(f"  Expression normalized via INT (mean~0, std~1 per gene)")
    print(f"  Expression data: {expression.shape}")
    print(f"  Genotype data: {genotypes.shape}")
    print(f"  Covariates data: {covariates.shape}")

except FileNotFoundError as e:
    print(f"ERROR: Could not load data files: {e}")
    exit(1)

if not expression.index.equals(genotypes.index):
    print("ERROR: Sample IDs do not match across datasets!")
    exit(1)

# ============================================================================
# 2. REGRESS OUT COVARIATES
# ============================================================================

print("\nStep 2: Regressing out covariates from expression...\n")

# Calculate genotype PCs for population structure
scaler = StandardScaler()
genotypes_scaled = scaler.fit_transform(genotypes)
pca = PCA(n_components=3)
pcs = pca.fit_transform(genotypes_scaled)

# Build covariate matrix: Age, Sex, PC1-3
cov_matrix = np.column_stack([
    covariates['Age'].values,
    covariates['Sex'].values,
    pcs
])

print(f"  Covariates: Age, Sex, PC1, PC2, PC3")

# Residualize expression
residualized = pd.DataFrame(index=expression.index, columns=expression.columns, dtype=float)
for gene in expression.columns:
    y = expression[gene].values.astype(float)
    reg = LinearRegression().fit(cov_matrix, y)
    residualized[gene] = y - reg.predict(cov_matrix)

print(f"  Residualized expression for {len(expression.columns)} genes")

# ============================================================================
# 3. SAVE DATA FOR R
# ============================================================================

print("\nStep 3: Preparing data for SuSiE R script...\n")

data_dir = os.path.abspath('data/ASD_dataset')

# Save genotype matrix (samples x SNPs) - R will read as matrix X
genotypes.to_csv(os.path.join(data_dir, 'temp_susie_genotypes.csv'))

# Save residualized expression
residualized.to_csv(os.path.join(data_dir, 'temp_susie_expression.csv'))

print(f"  Saved genotype matrix: {genotypes.shape}")
print(f"  Saved residualized expression: {residualized.shape}")

# ============================================================================
# 4. CREATE R SCRIPT FOR SUSIE
# ============================================================================

print("\nStep 4: Creating R script for SuSiE...\n")

r_script = r"""
library(susieR)

cat("SuSiE Fine-Mapping Analysis\n")
cat("==================================================\n\n")

# Load data
genotypes <- read.csv('/data/temp_susie_genotypes.csv', row.names=1)
expression <- read.csv('/data/temp_susie_expression.csv', row.names=1)

X <- as.matrix(genotypes)
storage.mode(X) <- "double"
genes <- colnames(expression)
snps <- colnames(genotypes)
n_samples <- nrow(X)

cat(sprintf("Loaded: %d samples, %d SNPs, %d genes\n\n", n_samples, ncol(X), length(genes)))

# Results storage
all_results <- data.frame()
cs_results <- data.frame()
gene_summary <- data.frame()

for (gene in genes) {
    y <- as.numeric(expression[[gene]])

    cat(sprintf("Fine-mapping %s...\n", gene))

    # Run SuSiE
    fit <- tryCatch({
        # Use scaled_prior_variance to set a fixed prior (0.2 * var(y))
        # rather than estimating it, which shrinks to 0 when signals are weak.
        # This is appropriate for educational/exploratory fine-mapping.
        susie(X, y,
              L = 5,
              scaled_prior_variance = 0.2,
              estimate_residual_variance = TRUE,
              estimate_prior_variance = FALSE,
              min_abs_corr = 0.0,
              coverage = 0.5,
              verbose = FALSE)
    }, error = function(e) {
        cat(sprintf("  WARNING: SuSiE failed for %s: %s\n", gene, e$message))
        return(NULL)
    })

    if (is.null(fit)) {
        # Record gene with no results
        gene_summary <- rbind(gene_summary, data.frame(
            Gene = gene,
            N_Credible_Sets = 0,
            Top_SNP = NA,
            Top_PIP = 0,
            Converged = FALSE
        ))
        # Still record PIPs (all zero)
        for (j in seq_along(snps)) {
            all_results <- rbind(all_results, data.frame(
                Gene = gene,
                SNP = snps[j],
                PIP = 0,
                Effect = 0,
                CS_ID = NA,
                stringsAsFactors = FALSE
            ))
        }
        next
    }

    # Extract PIPs
    pips <- fit$pip

    # Extract credible sets
    cs <- fit$sets
    n_cs <- 0
    if (!is.null(cs) && !is.null(cs$cs) && length(cs$cs) > 0) {
        n_cs <- length(cs$cs)
    }

    # Build per-SNP results
    for (j in seq_along(snps)) {
        cs_id <- NA
        if (n_cs > 0) {
            for (k in seq_along(cs$cs)) {
                if (j %in% cs$cs[[k]]) {
                    cs_id <- k
                    break
                }
            }
        }

        # Effect size from posterior mean
        effect <- coef(fit)[j + 1]  # +1 because coef includes intercept

        all_results <- rbind(all_results, data.frame(
            Gene = gene,
            SNP = snps[j],
            PIP = pips[j],
            Effect = effect,
            CS_ID = ifelse(is.na(cs_id), NA, cs_id),
            stringsAsFactors = FALSE
        ))
    }

    # Build credible set details
    if (n_cs > 0) {
        for (k in seq_along(cs$cs)) {
            cs_snp_indices <- cs$cs[[k]]
            cs_coverage <- cs$coverage[k]
            cs_min_corr <- ifelse(!is.null(cs$purity) && nrow(cs$purity) >= k,
                                   cs$purity[k, "min.abs.corr"], NA)

            for (idx in cs_snp_indices) {
                cs_results <- rbind(cs_results, data.frame(
                    Gene = gene,
                    CS_ID = k,
                    SNP = snps[idx],
                    PIP = pips[idx],
                    Effect = coef(fit)[idx + 1],
                    CS_Size = length(cs_snp_indices),
                    CS_Coverage = cs_coverage,
                    CS_Min_Corr = cs_min_corr,
                    stringsAsFactors = FALSE
                ))
            }
        }
    }

    # Gene summary
    top_idx <- which.max(pips)
    gene_summary <- rbind(gene_summary, data.frame(
        Gene = gene,
        N_Credible_Sets = n_cs,
        Top_SNP = snps[top_idx],
        Top_PIP = pips[top_idx],
        Converged = fit$converged
    ))

    if (n_cs > 0) {
        cat(sprintf("  -> %d credible set(s) found, top PIP = %.3f (%s)\n",
                     n_cs, pips[top_idx], snps[top_idx]))
    } else {
        cat(sprintf("  -> No credible sets, top PIP = %.3f (%s)\n",
                     pips[top_idx], snps[top_idx]))
    }
}

# Save results
write.csv(all_results, '/data/Phase2_SuSiE_results.csv', row.names=FALSE)
write.csv(cs_results, '/data/Phase2_SuSiE_credible_sets.csv', row.names=FALSE)
write.csv(gene_summary, '/data/Phase2_SuSiE_gene_summary.csv', row.names=FALSE)

cat("\n==================================================\n")
cat("SuSiE Analysis Complete!\n")
cat(sprintf("Total genes analyzed: %d\n", length(genes)))
cat(sprintf("Genes with credible sets: %d\n", sum(gene_summary$N_Credible_Sets > 0)))
cat(sprintf("Total credible sets found: %d\n", sum(gene_summary$N_Credible_Sets)))
cat(sprintf("Total SNPs in credible sets: %d\n", nrow(cs_results)))
cat("Results saved to Phase2_SuSiE_*.csv\n")
"""

with open(os.path.join(data_dir, 'run_susie.R'), 'w') as f:
    f.write(r_script)

print("  R script created: 'data/ASD_dataset/run_susie.R'")

# ============================================================================
# 5. RUN SUSIE IN DOCKER
# ============================================================================

print("\nStep 5: Running SuSiE analysis in Docker...\n")
print("Using Docker image: eqtlcatalogue/susie-finemapping\n")

# The R script uses /data/ paths inside Docker; for local R, rewrite paths
r_script_local = r_script.replace('/data/', f'{data_dir}/')

with open(os.path.join(data_dir, 'run_susie_local.R'), 'w') as f:
    f.write(r_script_local)

docker_cmd = [
    'docker', 'run', '--rm',
    '-v', f'{data_dir}:/data',
    'eqtlcatalogue/susie-finemapping',
    'Rscript', '/data/run_susie.R'
]

succeeded = False

# Try Docker first
try:
    print(f"Attempting Docker: {' '.join(docker_cmd)}\n")
    result = subprocess.run(docker_cmd, capture_output=True, text=True, timeout=3600)
    print(result.stdout)
    if result.returncode == 0:
        succeeded = True
    else:
        print("Docker run failed:")
        print(result.stderr)
except (FileNotFoundError, subprocess.TimeoutExpired, Exception) as e:
    print(f"Docker not available: {e}")

# Fallback: local Rscript
if not succeeded:
    print("\nFalling back to local Rscript...\n")
    try:
        result = subprocess.run(
            ['Rscript', os.path.join(data_dir, 'run_susie_local.R')],
            capture_output=True,
            text=True,
            timeout=3600
        )
        print(result.stdout)
        if result.returncode != 0:
            print("Local Rscript errors:")
            print(result.stderr)
            raise Exception("SuSiE analysis failed")
        succeeded = True
    except FileNotFoundError:
        print("ERROR: Neither Docker nor local Rscript available.")
        exit(1)
    except subprocess.TimeoutExpired:
        print("ERROR: SuSiE analysis timed out")
        exit(1)
    except Exception as e:
        print(f"ERROR: {str(e)}")
        exit(1)

# Clean up local R script
local_r_path = os.path.join(data_dir, 'run_susie_local.R')
if os.path.exists(local_r_path):
    os.remove(local_r_path)

print("\n  SuSiE analysis completed successfully")

# ============================================================================
# 6. LOAD AND ANALYZE RESULTS
# ============================================================================

print("\nStep 6: Loading SuSiE results...\n")

try:
    susie_results = pd.read_csv(os.path.join(data_dir, 'Phase2_SuSiE_results.csv'))
    cs_details = pd.read_csv(os.path.join(data_dir, 'Phase2_SuSiE_credible_sets.csv'))
    gene_summary = pd.read_csv(os.path.join(data_dir, 'Phase2_SuSiE_gene_summary.csv'))

    print(f"  Total SNP-gene results: {len(susie_results)}")
    print(f"  SNPs in credible sets: {len(cs_details)}")
    print(f"  Genes with credible sets: {(gene_summary['N_Credible_Sets'] > 0).sum()}/{len(gene_summary)}")

except FileNotFoundError as e:
    print(f"ERROR: Results file not found: {e}")
    exit(1)

# Load MatrixEQTL results for comparison
try:
    matrixeqtl_results = pd.read_csv(os.path.join(data_dir, 'Phase2_MatrixEQTL_results.csv'))
    has_matrixeqtl = True
    print(f"  MatrixEQTL results loaded: {len(matrixeqtl_results)} tests")
except FileNotFoundError:
    has_matrixeqtl = False
    print("  MatrixEQTL results not found (comparison plots will be skipped)")

# ============================================================================
# 7. VISUALIZATIONS
# ============================================================================

print("\nStep 7: Creating visualizations...\n")

# --- Plot 1: PIP Manhattan Plot ---
print("  Generating PIP Manhattan plot...")

genes_with_cs = gene_summary[gene_summary['N_Credible_Sets'] > 0]['Gene'].tolist()
# Show top 6 genes by max PIP (whether or not they have credible sets)
top_genes = (susie_results.groupby('Gene')['PIP'].max()
             .sort_values(ascending=False).head(6).index.tolist())

n_plots = len(top_genes)
if n_plots > 0:
    n_cols = min(3, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
    if n_plots == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    for idx, gene in enumerate(top_genes):
        ax = axes[idx]
        gene_data = susie_results[susie_results['Gene'] == gene].copy()
        gene_data = gene_data.sort_values('SNP')

        # Extract SNP number for x-axis positioning
        gene_data['SNP_num'] = gene_data['SNP'].str.extract(r'(\d+)').astype(int)
        gene_data = gene_data.sort_values('SNP_num')

        # Color by credible set membership
        colors = []
        cs_palette = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00']
        for _, row in gene_data.iterrows():
            if pd.isna(row['CS_ID']):
                colors.append('gray')
            else:
                cs_idx = int(row['CS_ID']) - 1
                colors.append(cs_palette[cs_idx % len(cs_palette)])

        ax.scatter(gene_data['SNP_num'], gene_data['PIP'],
                   c=colors, s=60, edgecolors='black', linewidths=0.5, zorder=3)
        ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.5, label='PIP = 0.5')
        ax.axhline(y=0.1, color='orange', linestyle=':', alpha=0.5, label='PIP = 0.1')

        ax.set_xlabel('SNP Index', fontsize=10)
        ax.set_ylabel('Posterior Inclusion Probability', fontsize=10)
        has_cs = gene in genes_with_cs
        n_cs = gene_summary[gene_summary['Gene'] == gene]['N_Credible_Sets'].values[0]
        ax.set_title(f'{gene}: {n_cs} credible set(s)', fontsize=12, fontweight='bold')
        ax.set_ylim(-0.02, 1.02)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    for idx in range(n_plots, len(axes)):
        axes[idx].axis('off')

    plt.suptitle('SuSiE Fine-Mapping: Posterior Inclusion Probabilities',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(data_dir, 'Phase2_SuSiE_pip_manhattan.png'),
                dpi=300, bbox_inches='tight')
    plt.close()
    print("    Saved: Phase2_SuSiE_pip_manhattan.png")

# --- Plot 2: Credible Set Summary ---
print("  Generating credible set summary plot...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Bar chart of credible sets per gene
ax = axes[0]
summary_sorted = gene_summary.sort_values('N_Credible_Sets', ascending=True)
bars = ax.barh(summary_sorted['Gene'], summary_sorted['N_Credible_Sets'],
               color=['steelblue' if n > 0 else 'lightgray'
                      for n in summary_sorted['N_Credible_Sets']],
               edgecolor='black')

for i, (_, row) in enumerate(summary_sorted.iterrows()):
    if row['N_Credible_Sets'] > 0:
        ax.text(row['N_Credible_Sets'] + 0.1, i,
                f"top: {row['Top_SNP']} (PIP={row['Top_PIP']:.2f})",
                va='center', fontsize=8)

ax.set_xlabel('Number of Credible Sets', fontsize=11)
ax.set_title('Credible Sets per Gene', fontsize=12, fontweight='bold')
ax.grid(True, alpha=0.3, axis='x')

# Top PIPs across all genes
ax = axes[1]
top_pips = susie_results.nlargest(15, 'PIP')
colors = ['#e41a1c' if not pd.isna(cs) else 'gray' for cs in top_pips['CS_ID']]
ax.barh([f"{r['SNP']}-{r['Gene']}" for _, r in top_pips.iterrows()],
        top_pips['PIP'], color=colors, edgecolor='black')
ax.axvline(x=0.5, color='red', linestyle='--', alpha=0.5, label='PIP = 0.5')
ax.set_xlabel('Posterior Inclusion Probability', fontsize=11)
ax.set_title('Top 15 SNP-Gene Pairs by PIP', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig(os.path.join(data_dir, 'Phase2_SuSiE_credible_sets.png'),
            dpi=300, bbox_inches='tight')
plt.close()
print("    Saved: Phase2_SuSiE_credible_sets.png")

# --- Plot 3: SuSiE vs MatrixEQTL Comparison ---
if has_matrixeqtl:
    print("  Generating SuSiE vs MatrixEQTL comparison plot...")

    # Merge results
    merged = susie_results.merge(matrixeqtl_results[['SNP', 'Gene', 'P_Value', 'FDR']],
                                 on=['SNP', 'Gene'], how='inner')

    if len(merged) > 0:
        merged['neg_log10_p'] = -np.log10(merged['P_Value'].clip(lower=1e-300))

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Scatter: -log10(p) vs PIP
        ax = axes[0]
        in_cs = ~merged['CS_ID'].isna()
        ax.scatter(merged.loc[~in_cs, 'neg_log10_p'], merged.loc[~in_cs, 'PIP'],
                   alpha=0.4, s=30, color='gray', label='Not in credible set')
        if in_cs.any():
            ax.scatter(merged.loc[in_cs, 'neg_log10_p'], merged.loc[in_cs, 'PIP'],
                       alpha=0.8, s=80, color='red', edgecolors='black',
                       linewidths=0.5, label='In credible set', zorder=5)

        ax.set_xlabel('MatrixEQTL -log10(P-value)', fontsize=11)
        ax.set_ylabel('SuSiE PIP', fontsize=11)
        ax.set_title('MatrixEQTL P-value vs SuSiE PIP', fontsize=12, fontweight='bold')
        ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.3)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

        # Comparison summary
        ax = axes[1]
        categories = ['All pairs', 'MatrixEQTL\nnominal\n(P<0.05)',
                       'SuSiE\nPIP>0.1', 'In credible\nset']
        counts = [
            len(merged),
            (merged['P_Value'] < 0.05).sum(),
            (merged['PIP'] > 0.1).sum(),
            in_cs.sum()
        ]
        bars = ax.bar(categories, counts, color=['steelblue', '#ff7f0e', '#2ca02c', '#d62728'],
                       edgecolor='black')
        for bar, count in zip(bars, counts):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                    str(count), ha='center', fontweight='bold')
        ax.set_ylabel('Number of SNP-Gene Pairs', fontsize=11)
        ax.set_title('Fine-Mapping Refinement', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')

        plt.tight_layout()
        plt.savefig(os.path.join(data_dir, 'Phase2_SuSiE_vs_MatrixEQTL.png'),
                    dpi=300, bbox_inches='tight')
        plt.close()
        print("    Saved: Phase2_SuSiE_vs_MatrixEQTL.png")

# --- Plot 4: LD and Credible Sets Heatmap ---
print("  Generating LD + credible sets heatmap...")

# Pick the gene with the most interesting SuSiE results
best_gene = gene_summary.sort_values(['N_Credible_Sets', 'Top_PIP'],
                                      ascending=[False, False]).iloc[0]['Gene']

fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# LD heatmap (r^2 between SNPs)
ax = axes[0]
geno_matrix = genotypes.values.astype(float)
ld_matrix = np.corrcoef(geno_matrix.T) ** 2

sns.heatmap(ld_matrix, ax=ax, cmap='YlOrRd', vmin=0, vmax=1,
            xticklabels=5, yticklabels=5,
            cbar_kws={'label': 'r²'})
ax.set_title('Linkage Disequilibrium (r²) Between SNPs', fontsize=12, fontweight='bold')
ax.set_xlabel('SNP Index')
ax.set_ylabel('SNP Index')

# PIP + credible set overlay for best gene
ax = axes[1]
gene_data = susie_results[susie_results['Gene'] == best_gene].copy()
gene_data['SNP_num'] = gene_data['SNP'].str.extract(r'(\d+)').astype(int)
gene_data = gene_data.sort_values('SNP_num')

# Build a matrix view: PIP as bar heights, color by CS
cs_palette = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00']
bar_colors = []
for _, row in gene_data.iterrows():
    if pd.isna(row['CS_ID']):
        bar_colors.append('lightgray')
    else:
        bar_colors.append(cs_palette[int(row['CS_ID']) - 1 % len(cs_palette)])

bars = ax.bar(gene_data['SNP_num'], gene_data['PIP'],
              color=bar_colors, edgecolor='black', linewidth=0.5)
ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.5, label='PIP = 0.5')
ax.set_xlabel('SNP Index', fontsize=11)
ax.set_ylabel('Posterior Inclusion Probability', fontsize=11)
n_cs = gene_summary[gene_summary['Gene'] == best_gene]['N_Credible_Sets'].values[0]
ax.set_title(f'SuSiE PIPs for {best_gene} ({n_cs} credible set(s))',
             fontsize=12, fontweight='bold')
ax.set_ylim(0, 1.05)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(os.path.join(data_dir, 'Phase2_SuSiE_ld_credible_sets.png'),
            dpi=300, bbox_inches='tight')
plt.close()
print("    Saved: Phase2_SuSiE_ld_credible_sets.png")

# ============================================================================
# 8. CLEAN UP TEMPORARY FILES
# ============================================================================

print("\nStep 8: Cleaning up temporary files...\n")

temp_files = [
    os.path.join(data_dir, 'temp_susie_genotypes.csv'),
    os.path.join(data_dir, 'temp_susie_expression.csv'),
    os.path.join(data_dir, 'run_susie.R'),
]

for temp_file in temp_files:
    if os.path.exists(temp_file):
        os.remove(temp_file)
        print(f"  Removed {os.path.basename(temp_file)}")

# ============================================================================
# 9. FINAL SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("SuSiE FINE-MAPPING ANALYSIS COMPLETE")
print("=" * 80)

print("\nSummary Statistics:")
print(f"  Genes analyzed: {len(gene_summary)}")
print(f"  Genes with credible sets: {(gene_summary['N_Credible_Sets'] > 0).sum()}")
print(f"  Total credible sets: {gene_summary['N_Credible_Sets'].sum()}")
print(f"  SNPs in credible sets: {len(cs_details)}")
if len(cs_details) > 0:
    print(f"  Average credible set size: {cs_details.groupby(['Gene', 'CS_ID']).size().mean():.1f}")

print("\nPer-Gene Summary:")
print("=" * 80)
print(gene_summary.to_string(index=False))

if len(cs_details) > 0:
    print("\nCredible Set Details:")
    print("=" * 80)
    print(cs_details.to_string(index=False))

if has_matrixeqtl:
    print("\n" + "=" * 80)
    print("COMPARISON WITH MatrixEQTL:")
    print("=" * 80)

    merged = susie_results.merge(matrixeqtl_results[['SNP', 'Gene', 'P_Value', 'FDR']],
                                 on=['SNP', 'Gene'], how='inner')
    if len(merged) > 0:
        nominal_sig = merged[merged['P_Value'] < 0.05]
        in_cs = merged[~merged['CS_ID'].isna()]

        print(f"\n  MatrixEQTL nominal significant (P<0.05): {len(nominal_sig)}")
        print(f"  SuSiE SNPs with PIP > 0.1: {(merged['PIP'] > 0.1).sum()}")
        print(f"  SuSiE SNPs in credible sets: {len(in_cs)}")

        if len(nominal_sig) > 0 and len(in_cs) > 0:
            overlap = nominal_sig.merge(in_cs[['SNP', 'Gene']], on=['SNP', 'Gene'])
            print(f"  Overlap (nominal sig AND in credible set): {len(overlap)}")

    print("\n  Key insight: SuSiE narrows broad association signals to small")
    print("  credible sets of variants most likely to be causal.")
    print("  MatrixEQTL finds many correlated SNPs; SuSiE resolves LD.")

print("\n" + "=" * 80)
print("OUTPUT FILES:")
print("=" * 80)
print("  1. Phase2_SuSiE_results.csv - Full results (all SNPs, all genes, PIPs)")
print("  2. Phase2_SuSiE_credible_sets.csv - Credible set membership")
print("  3. Phase2_SuSiE_gene_summary.csv - Per-gene summary")
print("  4. Phase2_SuSiE_pip_manhattan.png - PIP Manhattan plot")
print("  5. Phase2_SuSiE_credible_sets.png - Credible set summary")
if has_matrixeqtl:
    print("  6. Phase2_SuSiE_vs_MatrixEQTL.png - Comparison scatter plot")
print("  7. Phase2_SuSiE_ld_credible_sets.png - LD + credible sets heatmap")

print("\n" + "=" * 80)
print("ANALYSIS COMPLETE!")
print("=" * 80)
