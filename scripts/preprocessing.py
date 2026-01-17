#!/usr/bin/env python3
"""
Multiomics Data Preprocessing Pipeline for eQTL Analysis

This script processes multi-sample VCF, gene expression, and covariate data
to prepare inputs for eQTL mapping tools like MatrixEQTL.

Phases:
1. Sample Validation - Match samples across all input files
2. DNA Processing - Filter variants and compute genotype PCA
3. RNA Processing - Normalize expression and compute PEER factors
4. Final Validation - Verify output consistency

python scripts/preprocessing.py \
      --vcf data/raw_dummy_dataset/variants.vcf.gz \
      --expression data/raw_dummy_dataset/expression_counts.csv \
      --covariates data/raw_dummy_dataset/covariates.csv \
      --gtf data/raw_dummy_dataset/genes.gtf \
      --output-dir data/processed
      
Author: Vladimir Kovacevic
"""

import argparse
import gzip
import json
import logging
import os
import subprocess
import sys
import tempfile
from collections import OrderedDict
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import rankdata

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - [%(name)s] - %(message)s'
)
logger = logging.getLogger('preprocessing')


# =============================================================================
# Configuration and Data Classes
# =============================================================================

@dataclass
class PipelineConfig:
    """Configuration parameters for the preprocessing pipeline."""
    # Input files
    vcf_file: str
    expression_file: str
    covariates_file: str
    gtf_file: str

    # Output directory
    output_dir: str = "data/processed"

    # DNA processing parameters
    maf_threshold: float = 0.05
    call_rate_threshold: float = 0.95
    hwe_pvalue: float = 1e-6
    n_pca_components: Optional[int] = None  # Auto-select if None

    # RNA processing parameters
    cpm_threshold: float = 1.0
    min_samples_expressed: float = 0.2  # 20% of samples
    n_peer_factors: Optional[int] = None  # Auto-select if None
    protein_coding_only: bool = False

    # Docker image for PLINK2
    plink2_docker: str = "miguelpmachado/plink_2.0"

    def to_dict(self) -> dict:
        """Convert config to dictionary for logging."""
        return {
            'vcf_file': self.vcf_file,
            'expression_file': self.expression_file,
            'covariates_file': self.covariates_file,
            'gtf_file': self.gtf_file,
            'output_dir': self.output_dir,
            'maf_threshold': self.maf_threshold,
            'call_rate_threshold': self.call_rate_threshold,
            'hwe_pvalue': self.hwe_pvalue,
            'n_pca_components': self.n_pca_components,
            'cpm_threshold': self.cpm_threshold,
            'min_samples_expressed': self.min_samples_expressed,
            'n_peer_factors': self.n_peer_factors,
            'protein_coding_only': self.protein_coding_only,
            'plink2_docker': self.plink2_docker,
        }


@dataclass
class QCReport:
    """Quality control report tracking filtering steps."""
    initial_samples: int = 0
    initial_variants: int = 0
    initial_genes: int = 0

    # Sample filtering
    samples_in_vcf: int = 0
    samples_in_expression: int = 0
    samples_in_covariates: int = 0
    samples_common: int = 0

    # Variant filtering
    variants_biallelic: int = 0
    variants_autosome: int = 0
    variants_maf_pass: int = 0
    variants_call_rate_pass: int = 0
    variants_hwe_pass: int = 0
    variants_final: int = 0

    # Gene filtering
    genes_matched_gtf: int = 0
    genes_protein_coding: int = 0
    genes_expressed: int = 0
    genes_final: int = 0

    # PCA info
    pca_components_used: int = 0
    pca_variance_explained: List[float] = field(default_factory=list)

    # PEER info
    peer_factors_used: int = 0

    # Expression outliers
    expression_outliers: List[str] = field(default_factory=list)

    # Covariate info
    covariate_vif_warnings: List[str] = field(default_factory=list)

    # Parameters used
    parameters: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        """Convert report to dictionary."""
        return {
            'samples': {
                'initial_vcf': self.samples_in_vcf,
                'initial_expression': self.samples_in_expression,
                'initial_covariates': self.samples_in_covariates,
                'common': self.samples_common,
            },
            'variants': {
                'initial': self.initial_variants,
                'biallelic': self.variants_biallelic,
                'autosome': self.variants_autosome,
                'maf_pass': self.variants_maf_pass,
                'call_rate_pass': self.variants_call_rate_pass,
                'hwe_pass': self.variants_hwe_pass,
                'final': self.variants_final,
            },
            'genes': {
                'initial': self.initial_genes,
                'matched_gtf': self.genes_matched_gtf,
                'protein_coding': self.genes_protein_coding,
                'expressed': self.genes_expressed,
                'final': self.genes_final,
            },
            'pca': {
                'components_used': self.pca_components_used,
                'variance_explained': self.pca_variance_explained,
            },
            'peer': {
                'factors_used': self.peer_factors_used,
            },
            'expression_outliers': self.expression_outliers,
            'covariate_vif_warnings': self.covariate_vif_warnings,
            'parameters': self.parameters,
        }


# =============================================================================
# Phase 1: Sample Validation
# =============================================================================

class SampleValidator:
    """Validates and harmonizes samples across input files."""

    def __init__(self, config: PipelineConfig, qc_report: QCReport):
        self.config = config
        self.qc = qc_report
        self.logger = logging.getLogger('preprocessing.sample_validation')

    def extract_vcf_samples(self) -> Set[str]:
        """Extract sample IDs from VCF file."""
        self.logger.info("Extracting sample IDs from VCF...")

        vcf_path = self.config.vcf_file
        open_func = gzip.open if vcf_path.endswith('.gz') else open
        mode = 'rt' if vcf_path.endswith('.gz') else 'r'

        with open_func(vcf_path, mode) as f:
            for line in f:
                if line.startswith('#CHROM'):
                    parts = line.strip().split('\t')
                    samples = set(parts[9:])
                    self.logger.info(f"Found {len(samples)} samples in VCF")
                    return samples

        raise ValueError("Could not find header line in VCF file")

    def extract_expression_samples(self) -> Set[str]:
        """Extract sample IDs from expression matrix (samples as rows)."""
        self.logger.info("Extracting sample IDs from expression matrix...")

        # Read only the first column (index) to get sample IDs efficiently
        df = pd.read_csv(self.config.expression_file, usecols=[0])
        samples = set(df.iloc[:, 0].astype(str))
        self.logger.info(f"Found {len(samples)} samples in expression matrix")
        return samples

    def extract_covariate_samples(self) -> Set[str]:
        """Extract sample IDs from covariates file."""
        self.logger.info("Extracting sample IDs from covariates file...")

        df = pd.read_csv(self.config.covariates_file)
        # Assume first column is sample ID
        samples = set(df.iloc[:, 0].astype(str))
        self.logger.info(f"Found {len(samples)} samples in covariates file")
        return samples

    def detect_genome_build(self) -> Tuple[str, str]:
        """Detect genome build from VCF and GTF files."""
        self.logger.info("Detecting genome build...")

        # Check VCF header for reference
        vcf_build = None
        vcf_path = self.config.vcf_file
        open_func = gzip.open if vcf_path.endswith('.gz') else open
        mode = 'rt' if vcf_path.endswith('.gz') else 'r'

        with open_func(vcf_path, mode) as f:
            for line in f:
                if not line.startswith('#'):
                    break
                if 'GRCh38' in line or 'hg38' in line:
                    vcf_build = 'GRCh38'
                    break
                elif 'GRCh37' in line or 'hg19' in line:
                    vcf_build = 'GRCh37'
                    break

        # Check GTF header
        gtf_build = None
        with open(self.config.gtf_file, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    break
                if 'GRCh38' in line or 'hg38' in line:
                    gtf_build = 'GRCh38'
                    break
                elif 'GRCh37' in line or 'hg19' in line:
                    gtf_build = 'GRCh37'
                    break

        vcf_build = vcf_build or 'unknown'
        gtf_build = gtf_build or 'unknown'

        self.logger.info(f"Detected builds - VCF: {vcf_build}, GTF: {gtf_build}")

        if vcf_build != gtf_build and vcf_build != 'unknown' and gtf_build != 'unknown':
            self.logger.warning(f"Genome build mismatch! VCF: {vcf_build}, GTF: {gtf_build}")

        return vcf_build, gtf_build

    def detect_chromosome_naming(self) -> Tuple[str, str]:
        """Detect chromosome naming convention (chr1 vs 1)."""
        self.logger.info("Detecting chromosome naming convention...")

        # Check VCF
        vcf_chr_prefix = None
        vcf_path = self.config.vcf_file
        open_func = gzip.open if vcf_path.endswith('.gz') else open
        mode = 'rt' if vcf_path.endswith('.gz') else 'r'

        with open_func(vcf_path, mode) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                chrom = line.split('\t')[0]
                vcf_chr_prefix = 'chr' if chrom.startswith('chr') else 'num'
                break

        # Check GTF
        gtf_chr_prefix = None
        with open(self.config.gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                chrom = line.split('\t')[0]
                gtf_chr_prefix = 'chr' if chrom.startswith('chr') else 'num'
                break

        self.logger.info(f"Chromosome naming - VCF: {vcf_chr_prefix}, GTF: {gtf_chr_prefix}")
        return vcf_chr_prefix, gtf_chr_prefix

    def run(self) -> Tuple[List[str], str]:
        """
        Run sample validation phase.

        Returns:
            Tuple of (common_samples, chromosome_prefix_to_use)
        """
        self.logger.info("=" * 60)
        self.logger.info("PHASE 1: Sample Validation")
        self.logger.info("=" * 60)

        # Extract samples from all files
        vcf_samples = self.extract_vcf_samples()
        expr_samples = self.extract_expression_samples()
        cov_samples = self.extract_covariate_samples()

        self.qc.samples_in_vcf = len(vcf_samples)
        self.qc.samples_in_expression = len(expr_samples)
        self.qc.samples_in_covariates = len(cov_samples)

        # Find common samples
        common_samples = vcf_samples & expr_samples & cov_samples
        self.qc.samples_common = len(common_samples)

        self.logger.info(f"Common samples across all files: {len(common_samples)}")

        # Check for significant sample loss
        for name, samples in [('VCF', vcf_samples), ('Expression', expr_samples), ('Covariates', cov_samples)]:
            missing = samples - common_samples
            pct_missing = len(missing) / len(samples) * 100 if samples else 0
            if pct_missing > 10:
                self.logger.warning(f"{name}: {len(missing)} samples ({pct_missing:.1f}%) not in common set")
            if missing and len(missing) <= 10:
                self.logger.info(f"{name} samples not in common: {missing}")

        if len(common_samples) == 0:
            raise ValueError("No common samples found across all input files!")

        # Detect genome build
        self.detect_genome_build()

        # Detect chromosome naming
        vcf_chr, gtf_chr = self.detect_chromosome_naming()
        chr_prefix = 'chr' if vcf_chr == 'chr' else ''

        # Sort samples for consistent ordering
        common_samples = sorted(common_samples)

        self.logger.info(f"Phase 1 complete. {len(common_samples)} samples will be processed.")

        return common_samples, chr_prefix


# =============================================================================
# Phase 2: DNA Processing
# =============================================================================

class DNAProcessor:
    """Processes VCF data for eQTL analysis."""

    def __init__(self, config: PipelineConfig, qc_report: QCReport, common_samples: List[str]):
        self.config = config
        self.qc = qc_report
        self.samples = common_samples
        self.logger = logging.getLogger('preprocessing.dna_processing')

    def parse_vcf(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Parse VCF file and apply initial filters.

        Returns:
            Tuple of (genotype_matrix, variant_info)
        """
        self.logger.info("Parsing VCF file...")

        vcf_path = self.config.vcf_file
        open_func = gzip.open if vcf_path.endswith('.gz') else open
        mode = 'rt' if vcf_path.endswith('.gz') else 'r'

        variants = []
        genotypes = []

        with open_func(vcf_path, mode) as f:
            for line in f:
                if line.startswith('##'):
                    continue
                if line.startswith('#CHROM'):
                    header = line.strip().split('\t')
                    vcf_samples = header[9:]
                    # Map sample indices
                    sample_indices = [vcf_samples.index(s) for s in self.samples if s in vcf_samples]
                    continue

                parts = line.strip().split('\t')
                chrom = parts[0]
                pos = int(parts[1])
                var_id = parts[2]
                ref = parts[3]
                alt = parts[4]
                qual = parts[5]
                filt = parts[6]
                info = parts[7]
                fmt = parts[8]
                sample_data = parts[9:]

                # Skip multi-allelic
                if ',' in alt:
                    continue

                # Extract genotypes for common samples
                gt_index = fmt.split(':').index('GT')
                gts = []
                for idx in sample_indices:
                    gt_field = sample_data[idx].split(':')[gt_index]
                    # Convert to dosage (0, 1, 2, or NaN)
                    if gt_field in ['./.', '.|.', '.']:
                        gts.append(np.nan)
                    else:
                        alleles = gt_field.replace('|', '/').split('/')
                        dosage = sum(int(a) for a in alleles if a != '.')
                        gts.append(dosage)

                genotypes.append(gts)

                # Standardize variant ID
                std_id = f"{chrom}:{pos}:{ref}:{alt}"
                variants.append({
                    'variant_id': std_id,
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'original_id': var_id
                })

        # Create DataFrames
        var_df = pd.DataFrame(variants)
        gt_df = pd.DataFrame(genotypes, columns=self.samples)
        gt_df.index = var_df['variant_id']

        self.qc.initial_variants = len(var_df)
        self.qc.variants_biallelic = len(var_df)  # Already filtered

        self.logger.info(f"Parsed {len(var_df)} biallelic variants")

        return gt_df, var_df

    def filter_autosomes(self, gt_df: pd.DataFrame, var_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Filter to autosomes only (chr1-22)."""
        self.logger.info("Filtering to autosomes...")

        autosome_chroms = set([f"chr{i}" for i in range(1, 23)] + [str(i) for i in range(1, 23)])
        mask = var_df['chrom'].isin(autosome_chroms)

        var_df = var_df[mask].copy()
        gt_df = gt_df.loc[var_df['variant_id']].copy()

        self.qc.variants_autosome = len(var_df)
        self.logger.info(f"Variants on autosomes: {len(var_df)}")

        return gt_df, var_df

    def filter_maf(self, gt_df: pd.DataFrame, var_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Filter by minor allele frequency."""
        self.logger.info(f"Filtering by MAF >= {self.config.maf_threshold}...")

        # Calculate MAF
        alt_freq = gt_df.mean(axis=1) / 2
        maf = np.minimum(alt_freq, 1 - alt_freq)

        mask = maf >= self.config.maf_threshold
        var_df = var_df[mask.values].copy()
        gt_df = gt_df.loc[mask].copy()

        self.qc.variants_maf_pass = len(var_df)
        self.logger.info(f"Variants passing MAF filter: {len(var_df)}")

        return gt_df, var_df

    def filter_call_rate(self, gt_df: pd.DataFrame, var_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Filter by call rate."""
        self.logger.info(f"Filtering by call rate >= {self.config.call_rate_threshold}...")

        call_rate = gt_df.notna().mean(axis=1)
        mask = call_rate >= self.config.call_rate_threshold

        var_df = var_df[mask.values].copy()
        gt_df = gt_df.loc[mask].copy()

        self.qc.variants_call_rate_pass = len(var_df)
        self.logger.info(f"Variants passing call rate filter: {len(var_df)}")

        return gt_df, var_df

    def filter_hwe(self, gt_df: pd.DataFrame, var_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Filter by Hardy-Weinberg equilibrium."""
        self.logger.info(f"Filtering by HWE p-value >= {self.config.hwe_pvalue}...")

        hwe_pass = []

        for var_id in gt_df.index:
            gts = gt_df.loc[var_id].dropna().values

            # Count genotypes
            n0 = np.sum(gts == 0)
            n1 = np.sum(gts == 1)
            n2 = np.sum(gts == 2)
            n = n0 + n1 + n2

            if n == 0:
                hwe_pass.append(False)
                continue

            # Calculate allele frequencies
            p = (2 * n0 + n1) / (2 * n)
            q = 1 - p

            # Expected counts under HWE
            exp0 = p * p * n
            exp1 = 2 * p * q * n
            exp2 = q * q * n

            # Chi-square test
            obs = np.array([n0, n1, n2])
            exp = np.array([exp0, exp1, exp2])

            # Avoid division by zero
            exp = np.maximum(exp, 0.001)

            chi2 = np.sum((obs - exp) ** 2 / exp)
            p_value = 1 - stats.chi2.cdf(chi2, df=1)

            hwe_pass.append(p_value >= self.config.hwe_pvalue)

        mask = pd.Series(hwe_pass, index=gt_df.index)
        var_df = var_df[mask.values].copy()
        gt_df = gt_df.loc[mask].copy()

        self.qc.variants_hwe_pass = len(var_df)
        self.logger.info(f"Variants passing HWE filter: {len(var_df)}")

        return gt_df, var_df

    def compute_sample_qc(self, gt_df: pd.DataFrame) -> List[str]:
        """
        Compute sample-level QC metrics (heterozygosity).

        Returns list of samples flagged as potential outliers.
        """
        self.logger.info("Computing sample QC metrics...")

        flagged_samples = []

        # Calculate heterozygosity per sample
        het_rates = (gt_df == 1).sum(axis=0) / gt_df.notna().sum(axis=0)

        mean_het = het_rates.mean()
        std_het = het_rates.std()

        for sample in self.samples:
            if sample in het_rates:
                het = het_rates[sample]
                if abs(het - mean_het) > 3 * std_het:
                    flagged_samples.append(sample)
                    self.logger.warning(f"Sample {sample} has unusual heterozygosity: {het:.4f}")

        self.logger.info(f"Sample QC complete. {len(flagged_samples)} samples flagged.")
        return flagged_samples

    def compute_pca(self, gt_df: pd.DataFrame) -> pd.DataFrame:
        """
        Compute PCA on genotype data.

        Returns DataFrame with PCA components for each sample.
        """
        self.logger.info("Computing genotype PCA...")

        # Impute missing values with mean
        gt_imputed = gt_df.fillna(gt_df.mean())

        # Standardize
        gt_scaled = (gt_imputed - gt_imputed.mean()) / gt_imputed.std()
        gt_scaled = gt_scaled.fillna(0)

        # Compute PCA using SVD
        X = gt_scaled.values.T  # samples x variants

        # Center
        X_centered = X - X.mean(axis=0)

        # SVD
        U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)

        # Variance explained
        var_explained = (S ** 2) / (S ** 2).sum()

        # Determine number of components
        if self.config.n_pca_components is None:
            # Auto-select: keep components explaining > 1% variance, max 20
            n_components = min(20, sum(var_explained > 0.01))
            n_components = max(3, n_components)  # At least 3
        else:
            n_components = self.config.n_pca_components

        self.logger.info(f"Using {n_components} PCA components")

        # Store in QC report
        self.qc.pca_components_used = n_components
        self.qc.pca_variance_explained = var_explained[:n_components].tolist()

        # Create PCA DataFrame
        pca_df = pd.DataFrame(
            U[:, :n_components],
            index=self.samples,
            columns=[f'PC{i+1}' for i in range(n_components)]
        )

        self.logger.info(f"Top 5 PCs explain {sum(var_explained[:5])*100:.1f}% of variance")

        return pca_df

    def run(self) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Run DNA processing phase.

        Returns:
            Tuple of (genotype_dosage_matrix, variant_locations, pca_components)
        """
        self.logger.info("=" * 60)
        self.logger.info("PHASE 2: DNA Processing")
        self.logger.info("=" * 60)

        # Parse VCF
        gt_df, var_df = self.parse_vcf()

        # Apply filters
        gt_df, var_df = self.filter_autosomes(gt_df, var_df)
        gt_df, var_df = self.filter_maf(gt_df, var_df)
        gt_df, var_df = self.filter_call_rate(gt_df, var_df)
        gt_df, var_df = self.filter_hwe(gt_df, var_df)

        self.qc.variants_final = len(var_df)

        # Sample QC
        flagged = self.compute_sample_qc(gt_df)

        # Compute PCA
        pca_df = self.compute_pca(gt_df)

        # Create SNP location file
        snp_loc = var_df[['variant_id', 'chrom', 'pos']].copy()
        snp_loc.columns = ['snp_id', 'chr', 'pos']

        # Transpose genotype matrix: samples as rows, variants as columns
        gt_df_transposed = gt_df.T
        gt_df_transposed.index.name = 'sample_id'

        self.logger.info(f"Phase 2 complete. {len(gt_df_transposed.columns)} variants, {len(pca_df)} PCA components")

        return gt_df_transposed, snp_loc, pca_df


# =============================================================================
# Phase 3: RNA Processing
# =============================================================================

class RNAProcessor:
    """Processes RNA expression data for eQTL analysis."""

    def __init__(self, config: PipelineConfig, qc_report: QCReport, common_samples: List[str]):
        self.config = config
        self.qc = qc_report
        self.samples = common_samples
        self.logger = logging.getLogger('preprocessing.rna_processing')

    def load_expression(self) -> pd.DataFrame:
        """Load expression matrix (samples as rows, genes as columns)."""
        self.logger.info("Loading expression matrix...")

        df = pd.read_csv(self.config.expression_file, index_col=0)
        # Filter to common samples
        df = df.loc[self.samples]

        self.qc.initial_genes = len(df.columns)
        self.logger.info(f"Loaded expression for {len(df.columns)} genes, {len(df)} samples")

        return df

    def load_gtf(self) -> pd.DataFrame:
        """Load gene annotations from GTF."""
        self.logger.info("Loading GTF annotations...")

        genes = []
        with open(self.config.gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                if parts[2] != 'gene':
                    continue

                chrom = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attrs = parts[8]

                # Parse attributes
                gene_id = None
                gene_name = None
                gene_biotype = None

                for attr in attrs.split(';'):
                    attr = attr.strip()
                    if attr.startswith('gene_id'):
                        gene_id = attr.split('"')[1]
                    elif attr.startswith('gene_name'):
                        gene_name = attr.split('"')[1]
                    elif attr.startswith('gene_biotype'):
                        gene_biotype = attr.split('"')[1]

                if gene_id:
                    genes.append({
                        'gene_id': gene_id,
                        'gene_name': gene_name,
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'biotype': gene_biotype
                    })

        gtf_df = pd.DataFrame(genes)
        self.logger.info(f"Loaded {len(gtf_df)} genes from GTF")

        return gtf_df

    def match_gene_ids(self, expr_df: pd.DataFrame, gtf_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Match gene IDs between expression and GTF, handling version suffixes."""
        self.logger.info("Matching gene IDs...")

        # Expression has genes as columns
        expr_genes = set(expr_df.columns)
        gtf_genes = set(gtf_df['gene_id'])

        # Direct matches
        direct_matches = expr_genes & gtf_genes

        # Try removing version suffix from expression
        expr_no_version = {g.split('.')[0]: g for g in expr_genes}
        gtf_no_version = {g.split('.')[0]: g for g in gtf_genes}

        # Match by base ID
        matched_pairs = {}
        for base_id, expr_id in expr_no_version.items():
            if base_id in gtf_no_version:
                matched_pairs[expr_id] = gtf_no_version[base_id]
            elif expr_id in gtf_genes:
                matched_pairs[expr_id] = expr_id

        # Filter expression to matched genes (columns)
        matched_expr_genes = list(matched_pairs.keys())
        expr_df = expr_df[matched_expr_genes].copy()

        # Rename expression columns to match GTF if needed
        expr_df.columns = [matched_pairs[g] for g in expr_df.columns]

        # Filter GTF
        gtf_df = gtf_df[gtf_df['gene_id'].isin(expr_df.columns)].copy()

        self.qc.genes_matched_gtf = len(expr_df.columns)
        self.logger.info(f"Matched {len(expr_df.columns)} genes between expression and GTF")

        return expr_df, gtf_df

    def filter_protein_coding(self, expr_df: pd.DataFrame, gtf_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Filter to protein-coding genes only (genes are columns)."""
        if not self.config.protein_coding_only:
            self.qc.genes_protein_coding = len(expr_df.columns)
            return expr_df, gtf_df

        self.logger.info("Filtering to protein-coding genes...")

        protein_coding = gtf_df[gtf_df['biotype'] == 'protein_coding']['gene_id']
        genes_to_keep = [g for g in expr_df.columns if g in protein_coding.values]

        expr_df = expr_df[genes_to_keep].copy()
        gtf_df = gtf_df[gtf_df['gene_id'].isin(expr_df.columns)].copy()

        self.qc.genes_protein_coding = len(expr_df.columns)
        self.logger.info(f"Protein-coding genes: {len(expr_df.columns)}")

        return expr_df, gtf_df

    def filter_low_expression(self, expr_df: pd.DataFrame) -> pd.DataFrame:
        """Filter lowly expressed genes based on CPM threshold (samples as rows, genes as columns)."""
        self.logger.info(f"Filtering low-expression genes (CPM > {self.config.cpm_threshold} in {self.config.min_samples_expressed*100:.0f}% samples)...")

        # Calculate CPM: library size per sample (sum across genes = axis 1)
        lib_sizes = expr_df.sum(axis=1)
        cpm = expr_df.div(lib_sizes, axis=0) * 1e6

        # Filter genes (columns): count samples (rows) where gene passes threshold
        min_samples = int(len(self.samples) * self.config.min_samples_expressed)
        genes_pass = (cpm > self.config.cpm_threshold).sum(axis=0) >= min_samples

        genes_to_keep = genes_pass[genes_pass].index.tolist()
        expr_df = expr_df[genes_to_keep].copy()

        self.qc.genes_expressed = len(expr_df.columns)
        self.logger.info(f"Genes passing expression filter: {len(expr_df.columns)}")

        return expr_df

    def detect_expression_outliers(self, expr_df: pd.DataFrame) -> List[str]:
        """Detect expression outliers using PCA (samples as rows, genes as columns)."""
        self.logger.info("Detecting expression outliers...")

        # Log transform for PCA
        log_expr = np.log2(expr_df + 1)

        # Standardize each gene (column)
        log_scaled = (log_expr - log_expr.mean(axis=0)) / log_expr.std(axis=0)
        log_scaled = log_scaled.fillna(0)

        # PCA: X is already samples x genes
        X = log_scaled.values
        X_centered = X - X.mean(axis=0)

        U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)

        # Check first 3 PCs for outliers
        outliers = []
        sample_list = list(expr_df.index)
        for i in range(min(3, U.shape[1])):
            pc = U[:, i]
            mean_pc = pc.mean()
            std_pc = pc.std()

            for j, sample in enumerate(sample_list):
                if abs(pc[j] - mean_pc) > 3 * std_pc:
                    if sample not in outliers:
                        outliers.append(sample)
                        self.logger.warning(f"Sample {sample} is outlier on expression PC{i+1}")

        self.qc.expression_outliers = outliers
        self.logger.info(f"Expression outliers detected: {len(outliers)}")

        return outliers

    def normalize_expression(self, expr_df: pd.DataFrame) -> pd.DataFrame:
        """Apply log2(CPM+1) and TMM normalization (samples as rows, genes as columns)."""
        self.logger.info("Normalizing expression...")

        # Calculate CPM: library size per sample (sum across genes = axis 1)
        lib_sizes = expr_df.sum(axis=1)
        cpm = expr_df.div(lib_sizes, axis=0) * 1e6

        # Log2(CPM + 1)
        log_cpm = np.log2(cpm + 1)

        # Quantile normalization across samples (for each gene, rank samples)
        # Simplified approach: for each gene, transform to standard normal
        normalized = log_cpm.copy()

        # Quantile normalize each gene across samples
        for gene in normalized.columns:
            ranked = log_cpm[gene].rank(method='average')
            n = len(ranked)
            quantiles = (ranked - 0.5) / n
            normalized[gene] = stats.norm.ppf(quantiles)
            # Handle inf values
            normalized[gene] = normalized[gene].replace([np.inf, -np.inf], np.nan)
            normalized[gene] = normalized[gene].fillna(0)

        self.logger.info("Log2(CPM+1) and quantile normalization complete")

        return normalized

    def compute_peer_factors(self, expr_df: pd.DataFrame) -> pd.DataFrame:
        """
        Compute PEER factors for hidden confounder correction (samples as rows, genes as columns).

        Uses PCA as a PEER approximation since peer package requires R.
        """
        self.logger.info("Computing PEER factors (using PCA approximation)...")

        # Determine number of factors
        n_samples = len(expr_df)
        if self.config.n_peer_factors is None:
            if n_samples < 150:
                n_factors = 15
            elif n_samples < 250:
                n_factors = 30
            elif n_samples < 350:
                n_factors = 45
            else:
                n_factors = 60
        else:
            n_factors = self.config.n_peer_factors

        # Limit to available samples/genes
        n_factors = min(n_factors, n_samples - 1, len(expr_df.columns))

        self.logger.info(f"Using {n_factors} PEER factors for {n_samples} samples")

        # Standardize expression (per gene/column)
        expr_scaled = (expr_df - expr_df.mean(axis=0)) / expr_df.std(axis=0)
        expr_scaled = expr_scaled.fillna(0)

        # PCA: X is already samples x genes
        X = expr_scaled.values
        X_centered = X - X.mean(axis=0)

        U, _, _ = np.linalg.svd(X_centered, full_matrices=False)

        # PEER factors
        peer_df = pd.DataFrame(
            U[:, :n_factors],
            index=expr_df.index,
            columns=[f'PEER{i+1}' for i in range(n_factors)]
        )

        self.qc.peer_factors_used = n_factors

        return peer_df

    def inverse_normal_transform(self, expr_df: pd.DataFrame) -> pd.DataFrame:
        """Apply inverse normal transformation to expression values (samples as rows, genes as columns)."""
        self.logger.info("Applying inverse normal transformation...")

        def int_transform(x):
            """Inverse normal transform a single gene across samples."""
            ranks = rankdata(x, method='average')
            n = len(x)
            # Blom transform
            quantiles = (ranks - 0.375) / (n + 0.25)
            return stats.norm.ppf(quantiles)

        # Apply to each column (gene)
        int_df = expr_df.apply(int_transform, axis=0)

        self.logger.info("Inverse normal transformation complete")

        return int_df

    def run(self) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Run RNA processing phase.

        Returns:
            Tuple of (normalized_expression, gene_locations, peer_factors)
        """
        self.logger.info("=" * 60)
        self.logger.info("PHASE 3: RNA Processing")
        self.logger.info("=" * 60)

        # Load data
        expr_df = self.load_expression()
        gtf_df = self.load_gtf()

        # Match gene IDs
        expr_df, gtf_df = self.match_gene_ids(expr_df, gtf_df)

        # Filter protein-coding
        expr_df, gtf_df = self.filter_protein_coding(expr_df, gtf_df)

        # Filter low expression
        expr_df = self.filter_low_expression(expr_df)

        # Update GTF to match filtered genes (genes are columns)
        gtf_df = gtf_df[gtf_df['gene_id'].isin(expr_df.columns)].copy()

        # Detect outliers
        self.detect_expression_outliers(expr_df)

        # Normalize
        norm_expr = self.normalize_expression(expr_df)

        # Compute PEER factors
        peer_df = self.compute_peer_factors(norm_expr)

        # Inverse normal transform
        int_expr = self.inverse_normal_transform(norm_expr)

        # Create gene location file
        gene_loc = gtf_df[['gene_id', 'chrom', 'start', 'end']].copy()
        gene_loc.columns = ['gene_id', 'chr', 'start', 'end']

        self.qc.genes_final = len(int_expr.columns)

        self.logger.info(f"Phase 3 complete. {len(int_expr.columns)} genes processed")

        return int_expr, gene_loc, peer_df


# =============================================================================
# Phase 4: Final Validation
# =============================================================================

class FinalValidator:
    """Performs final validation checks on outputs."""

    def __init__(self, config: PipelineConfig, qc_report: QCReport):
        self.config = config
        self.qc = qc_report
        self.logger = logging.getLogger('preprocessing.final_validation')

    def check_multicollinearity(self, covariates_df: pd.DataFrame) -> List[str]:
        """Check for multicollinearity using VIF."""
        self.logger.info("Checking covariate multicollinearity...")

        warnings = []

        # Only check numeric columns
        numeric_cols = covariates_df.select_dtypes(include=[np.number]).columns.tolist()

        if len(numeric_cols) < 2:
            return warnings

        X = covariates_df[numeric_cols].values

        # Calculate VIF for each variable
        for i, col in enumerate(numeric_cols):
            # Get all other columns
            other_cols = [j for j in range(len(numeric_cols)) if j != i]
            if not other_cols:
                continue

            X_other = X[:, other_cols]
            y = X[:, i]

            # Add constant
            X_other = np.column_stack([np.ones(len(y)), X_other])

            # Fit linear regression
            try:
                beta = np.linalg.lstsq(X_other, y, rcond=None)[0]
                y_pred = X_other @ beta
                ss_res = np.sum((y - y_pred) ** 2)
                ss_tot = np.sum((y - y.mean()) ** 2)

                if ss_tot > 0:
                    r_squared = 1 - ss_res / ss_tot
                    if r_squared < 1:
                        vif = 1 / (1 - r_squared)
                        if vif > 5:
                            warnings.append(f"{col}: VIF = {vif:.2f}")
                            self.logger.warning(f"High VIF for {col}: {vif:.2f}")
            except Exception:
                pass

        self.qc.covariate_vif_warnings = warnings
        return warnings

    def verify_sample_consistency(self, genotypes: pd.DataFrame, expression: pd.DataFrame,
                                   covariates: pd.DataFrame) -> bool:
        """Verify samples are consistent across all output files (samples as rows/index)."""
        self.logger.info("Verifying sample consistency...")

        gt_samples = list(genotypes.index)
        expr_samples = list(expression.index)
        cov_samples = list(covariates.index)

        # Check count
        if not (len(gt_samples) == len(expr_samples) == len(cov_samples)):
            self.logger.error(f"Sample count mismatch: GT={len(gt_samples)}, Expr={len(expr_samples)}, Cov={len(cov_samples)}")
            return False

        # Check order
        if gt_samples != expr_samples or gt_samples != cov_samples:
            self.logger.error("Sample order mismatch across files!")
            return False

        self.logger.info(f"Sample consistency verified: {len(gt_samples)} samples in all files")
        return True

    def run(self, genotypes: pd.DataFrame, expression: pd.DataFrame,
            covariates: pd.DataFrame) -> bool:
        """
        Run final validation phase.

        Returns True if all checks pass.
        """
        self.logger.info("=" * 60)
        self.logger.info("PHASE 4: Final Validation")
        self.logger.info("=" * 60)

        # Check multicollinearity
        self.check_multicollinearity(covariates)

        # Verify sample consistency
        consistent = self.verify_sample_consistency(genotypes, expression, covariates)

        self.logger.info(f"Phase 4 complete. Validation {'passed' if consistent else 'FAILED'}")

        return consistent


# =============================================================================
# Main Pipeline
# =============================================================================

class PreprocessingPipeline:
    """Main preprocessing pipeline orchestrator."""

    def __init__(self, config: PipelineConfig):
        self.config = config
        self.qc = QCReport()
        self.logger = logging.getLogger('preprocessing.pipeline')

    def save_outputs(self, genotypes: pd.DataFrame, snp_loc: pd.DataFrame,
                     expression: pd.DataFrame, gene_loc: pd.DataFrame,
                     covariates: pd.DataFrame):
        """Save all output files."""
        self.logger.info("Saving output files...")

        os.makedirs(self.config.output_dir, exist_ok=True)

        # Genotype dosage matrix
        gt_path = os.path.join(self.config.output_dir, 'genotype_dosage.csv')
        genotypes.to_csv(gt_path)
        self.logger.info(f"Saved genotype dosage matrix: {gt_path}")

        # SNP locations
        snp_path = os.path.join(self.config.output_dir, 'snp_locations.csv')
        snp_loc.to_csv(snp_path, index=False)
        self.logger.info(f"Saved SNP locations: {snp_path}")

        # Expression matrix
        expr_path = os.path.join(self.config.output_dir, 'expression_normalized.csv')
        expression.to_csv(expr_path)
        self.logger.info(f"Saved normalized expression: {expr_path}")

        # Gene locations
        gene_path = os.path.join(self.config.output_dir, 'gene_locations.csv')
        gene_loc.to_csv(gene_path, index=False)
        self.logger.info(f"Saved gene locations: {gene_path}")

        # Covariates
        cov_path = os.path.join(self.config.output_dir, 'covariates_combined.csv')
        covariates.to_csv(cov_path)
        self.logger.info(f"Saved combined covariates: {cov_path}")

        # QC report
        self.qc.parameters = self.config.to_dict()
        qc_path = os.path.join(self.config.output_dir, 'qc_report.json')
        with open(qc_path, 'w') as f:
            json.dump(self.qc.to_dict(), f, indent=2)
        self.logger.info(f"Saved QC report: {qc_path}")

    def run(self):
        """Run the complete preprocessing pipeline."""
        start_time = datetime.now()

        self.logger.info("=" * 60)
        self.logger.info("MULTIOMICS PREPROCESSING PIPELINE")
        self.logger.info(f"Started: {start_time}")
        self.logger.info("=" * 60)

        # Phase 1: Sample Validation
        validator = SampleValidator(self.config, self.qc)
        common_samples, chr_prefix = validator.run()

        # Phase 2: DNA Processing
        dna_processor = DNAProcessor(self.config, self.qc, common_samples)
        genotypes, snp_loc, pca_df = dna_processor.run()

        # Phase 3: RNA Processing
        rna_processor = RNAProcessor(self.config, self.qc, common_samples)
        expression, gene_loc, peer_df = rna_processor.run()

        # Combine covariates
        self.logger.info("Combining covariates...")
        user_cov = pd.read_csv(self.config.covariates_file)
        user_cov = user_cov.set_index(user_cov.columns[0])
        user_cov = user_cov.loc[common_samples]

        # Encode categorical variables
        for col in user_cov.select_dtypes(include=['object']).columns:
            user_cov[col] = pd.Categorical(user_cov[col]).codes

        combined_cov = pd.concat([user_cov, pca_df, peer_df], axis=1)

        # Phase 4: Final Validation
        final_validator = FinalValidator(self.config, self.qc)
        final_validator.run(genotypes, expression, combined_cov)

        # Save outputs
        self.save_outputs(genotypes, snp_loc, expression, gene_loc, combined_cov)

        end_time = datetime.now()
        duration = end_time - start_time

        self.logger.info("=" * 60)
        self.logger.info("PIPELINE COMPLETE")
        self.logger.info(f"Duration: {duration}")
        self.logger.info(f"Final samples: {len(common_samples)}")
        self.logger.info(f"Final variants: {self.qc.variants_final}")
        self.logger.info(f"Final genes: {self.qc.genes_final}")
        self.logger.info(f"Output directory: {self.config.output_dir}")
        self.logger.info("=" * 60)


# =============================================================================
# CLI Entry Point
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Multiomics Preprocessing Pipeline for eQTL Analysis',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required inputs
    parser.add_argument('--vcf', '-v', required=True,
                        help='Multi-sample VCF file (can be gzipped)')
    parser.add_argument('--expression', '-e', required=True,
                        help='Raw gene expression counts (CSV, genes x samples)')
    parser.add_argument('--covariates', '-c', required=True,
                        help='Covariates CSV file')
    parser.add_argument('--gtf', '-g', required=True,
                        help='Gene annotation GTF file')

    # Output
    parser.add_argument('--output-dir', '-o', default='data/processed',
                        help='Output directory')

    # DNA processing parameters
    parser.add_argument('--maf', type=float, default=0.05,
                        help='Minor allele frequency threshold')
    parser.add_argument('--call-rate', type=float, default=0.95,
                        help='Minimum call rate threshold')
    parser.add_argument('--hwe-pvalue', type=float, default=1e-6,
                        help='HWE p-value threshold')
    parser.add_argument('--n-pca', type=int, default=None,
                        help='Number of PCA components (auto-select if not specified)')

    # RNA processing parameters
    parser.add_argument('--cpm-threshold', type=float, default=1.0,
                        help='CPM threshold for low expression filter')
    parser.add_argument('--min-samples', type=float, default=0.2,
                        help='Minimum fraction of samples for expression filter')
    parser.add_argument('--n-peer', type=int, default=None,
                        help='Number of PEER factors (auto-select if not specified)')
    parser.add_argument('--protein-coding-only', action='store_true',
                        help='Filter to protein-coding genes only')

    args = parser.parse_args()

    # Create config
    config = PipelineConfig(
        vcf_file=args.vcf,
        expression_file=args.expression,
        covariates_file=args.covariates,
        gtf_file=args.gtf,
        output_dir=args.output_dir,
        maf_threshold=args.maf,
        call_rate_threshold=args.call_rate,
        hwe_pvalue=args.hwe_pvalue,
        n_pca_components=args.n_pca,
        cpm_threshold=args.cpm_threshold,
        min_samples_expressed=args.min_samples,
        n_peer_factors=args.n_peer,
        protein_coding_only=args.protein_coding_only,
    )

    # Run pipeline
    pipeline = PreprocessingPipeline(config)
    pipeline.run()


if __name__ == '__main__':
    main()
