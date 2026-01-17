#!/usr/bin/env python3
"""
Create a simulated dummy dataset for testing the multiomics preprocessing pipeline.

Generates:
- VCF file with 200 samples and 100 biallelic SNPs on chromosomes 1-22
- Gene expression matrix (raw counts) with 200 samples and 100 genes
- GTF annotation file for the 100 genes
- Covariates CSV with age and sex

Author: Multiomics Pipeline
"""

import argparse
import gzip
import logging
import os
import random
from datetime import datetime

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def create_sample_ids(n_samples: int) -> list:
    """Generate sample IDs."""
    return [f"SAMPLE_{i:03d}" for i in range(1, n_samples + 1)]


def create_vcf(output_path: str, sample_ids: list, n_variants: int = 100, seed: int = 42):
    """
    Create a simulated VCF file with biallelic SNPs.

    Variants are distributed across chromosomes 1-22 with realistic MAFs.
    """
    logger.info(f"Creating VCF file with {len(sample_ids)} samples and {n_variants} variants")
    np.random.seed(seed)
    random.seed(seed)

    # VCF header
    header_lines = [
        '##fileformat=VCFv4.2',
        f'##fileDate={datetime.now().strftime("%Y%m%d")}',
        '##source=create_dummy_dataset.py',
        '##reference=GRCh38',
        '##contig=<ID=chr1,length=248956422>',
        '##contig=<ID=chr2,length=242193529>',
        '##contig=<ID=chr3,length=198295559>',
        '##contig=<ID=chr4,length=190214555>',
        '##contig=<ID=chr5,length=181538259>',
        '##contig=<ID=chr6,length=170805979>',
        '##contig=<ID=chr7,length=159345973>',
        '##contig=<ID=chr8,length=145138636>',
        '##contig=<ID=chr9,length=138394717>',
        '##contig=<ID=chr10,length=133797422>',
        '##contig=<ID=chr11,length=135086622>',
        '##contig=<ID=chr12,length=133275309>',
        '##contig=<ID=chr13,length=114364328>',
        '##contig=<ID=chr14,length=107043718>',
        '##contig=<ID=chr15,length=101991189>',
        '##contig=<ID=chr16,length=90338345>',
        '##contig=<ID=chr17,length=83257441>',
        '##contig=<ID=chr18,length=80373285>',
        '##contig=<ID=chr19,length=58617616>',
        '##contig=<ID=chr20,length=64444167>',
        '##contig=<ID=chr21,length=46709983>',
        '##contig=<ID=chr22,length=50818468>',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    ]

    # Column header
    column_header = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(sample_ids)

    # Generate variants distributed across chromosomes
    variants_per_chrom = n_variants // 22
    extra_variants = n_variants % 22

    variant_records = []
    nucleotides = ['A', 'C', 'G', 'T']

    for chrom in range(1, 23):
        n_vars = variants_per_chrom + (1 if chrom <= extra_variants else 0)

        # Random positions within reasonable chromosome ranges
        chrom_length = 100_000_000  # Simplified
        positions = sorted(np.random.randint(10000, chrom_length, size=n_vars))

        for pos in positions:
            # Random ref and alt alleles (different from each other)
            ref = random.choice(nucleotides)
            alt = random.choice([n for n in nucleotides if n != ref])

            # MAF between 0.05 and 0.5
            maf = np.random.uniform(0.05, 0.5)

            # Generate genotypes based on HWE
            p = 1 - maf  # Frequency of reference allele
            q = maf      # Frequency of alternate allele

            # Genotype probabilities under HWE
            prob_00 = p * p
            prob_01 = 2 * p * q
            prob_11 = q * q

            genotypes = []
            for _ in sample_ids:
                r = np.random.random()
                if r < prob_00:
                    gt = '0/0'
                elif r < prob_00 + prob_01:
                    gt = '0/1'
                else:
                    gt = '1/1'
                genotypes.append(gt)

            # Calculate actual AF from genotypes
            alt_count = sum(gt.count('1') for gt in genotypes)
            actual_af = alt_count / (2 * len(sample_ids))

            variant_id = f"chr{chrom}:{pos}:{ref}:{alt}"
            record = f"chr{chrom}\t{pos}\t{variant_id}\t{ref}\t{alt}\t.\tPASS\tAF={actual_af:.4f}\tGT\t" + '\t'.join(genotypes)
            variant_records.append(record)

    # Write VCF
    with gzip.open(output_path, 'wt') as f:
        f.write('\n'.join(header_lines) + '\n')
        f.write(column_header + '\n')
        f.write('\n'.join(variant_records) + '\n')

    logger.info(f"VCF file written to {output_path}")
    return len(variant_records)


def create_expression_matrix(output_path: str, sample_ids: list, gene_ids: list, seed: int = 42):
    """
    Create a simulated raw count expression matrix.

    Uses negative binomial distribution to simulate RNA-seq counts.
    Format: samples as rows, genes as columns.
    """
    logger.info(f"Creating expression matrix with {len(sample_ids)} samples and {len(gene_ids)} genes")
    np.random.seed(seed)

    n_samples = len(sample_ids)
    n_genes = len(gene_ids)

    # Simulate raw counts using negative binomial
    # Mean expression varies by gene (some highly expressed, some lowly)
    gene_means = np.random.lognormal(mean=5, sigma=2, size=n_genes)

    # Dispersion parameter
    dispersion = 0.1

    # counts array: samples x genes
    counts = np.zeros((n_samples, n_genes), dtype=int)

    for i, mean_expr in enumerate(gene_means):
        # Add sample-specific variation (library size effects)
        sample_factors = np.random.lognormal(mean=0, sigma=0.3, size=n_samples)
        adjusted_means = mean_expr * sample_factors

        # Negative binomial parameters
        # variance = mean + mean^2 / dispersion
        # p = dispersion / (dispersion + mean)
        # n = dispersion
        for j, adj_mean in enumerate(adjusted_means):
            if adj_mean > 0:
                p = dispersion / (dispersion + adj_mean)
                counts[j, i] = np.random.negative_binomial(dispersion, p)

    # Create DataFrame: samples as rows, genes as columns
    df = pd.DataFrame(counts, index=sample_ids, columns=gene_ids)
    df.index.name = 'sample_id'

    # Save to CSV
    df.to_csv(output_path)
    logger.info(f"Expression matrix written to {output_path} (samples x genes)")

    return df


def create_gtf(output_path: str, gene_ids: list, seed: int = 42):
    """
    Create a simulated GTF annotation file.

    Distributes genes across chromosomes 1-22 with realistic coordinates.
    """
    logger.info(f"Creating GTF file with {len(gene_ids)} genes")
    np.random.seed(seed)
    random.seed(seed)

    n_genes = len(gene_ids)
    genes_per_chrom = n_genes // 22
    extra_genes = n_genes % 22

    gtf_records = []
    gene_idx = 0

    for chrom in range(1, 23):
        n_genes_chrom = genes_per_chrom + (1 if chrom <= extra_genes else 0)

        # Random positions
        chrom_length = 100_000_000
        starts = sorted(np.random.randint(100000, chrom_length - 100000, size=n_genes_chrom))

        for start in starts:
            if gene_idx >= len(gene_ids):
                break

            gene_id = gene_ids[gene_idx]
            gene_name = f"GENE_{gene_idx + 1}"

            # Gene length between 1kb and 50kb
            gene_length = np.random.randint(1000, 50000)
            end = start + gene_length

            # Strand
            strand = random.choice(['+', '-'])

            # GTF format: seqname source feature start end score strand frame attributes
            attributes = f'gene_id "{gene_id}"; gene_name "{gene_name}"; gene_biotype "protein_coding";'
            record = f"chr{chrom}\tensembl\tgene\t{start}\t{end}\t.\t{strand}\t.\t{attributes}"
            gtf_records.append(record)

            # Add transcript
            transcript_id = f"{gene_id}.1"
            attributes = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; gene_name "{gene_name}"; gene_biotype "protein_coding";'
            record = f"chr{chrom}\tensembl\ttranscript\t{start}\t{end}\t.\t{strand}\t.\t{attributes}"
            gtf_records.append(record)

            gene_idx += 1

    # Write GTF (with header)
    with open(output_path, 'w') as f:
        f.write('#!genome-build GRCh38\n')
        f.write('#!genome-version GRCh38\n')
        f.write('\n'.join(gtf_records) + '\n')

    logger.info(f"GTF file written to {output_path}")


def create_covariates(output_path: str, sample_ids: list, seed: int = 42):
    """
    Create a simulated covariates CSV with age and sex.
    """
    logger.info(f"Creating covariates file with {len(sample_ids)} samples")
    np.random.seed(seed)

    n_samples = len(sample_ids)

    # Age: uniform between 20 and 80
    ages = np.random.randint(20, 81, size=n_samples)

    # Sex: roughly 50/50
    sexes = np.random.choice(['M', 'F'], size=n_samples)

    df = pd.DataFrame({
        'sample_id': sample_ids,
        'age': ages,
        'sex': sexes
    })

    df.to_csv(output_path, index=False)
    logger.info(f"Covariates file written to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Create simulated dummy dataset for multiomics preprocessing testing'
    )
    parser.add_argument(
        '--output-dir', '-o',
        default='data/raw_dummy_dataset',
        help='Output directory (default: data/raw_dummy_dataset)'
    )
    parser.add_argument(
        '--n-samples', '-n',
        type=int,
        default=200,
        help='Number of samples (default: 200)'
    )
    parser.add_argument(
        '--n-genes', '-g',
        type=int,
        default=100,
        help='Number of genes (default: 100)'
    )
    parser.add_argument(
        '--n-variants', '-v',
        type=int,
        default=100,
        help='Number of variants (default: 100)'
    )
    parser.add_argument(
        '--seed', '-s',
        type=int,
        default=42,
        help='Random seed (default: 42)'
    )

    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    logger.info(f"Output directory: {args.output_dir}")

    # Generate sample IDs
    sample_ids = create_sample_ids(args.n_samples)

    # Generate gene IDs (Ensembl-style)
    gene_ids = [f"ENSG{i:011d}" for i in range(1, args.n_genes + 1)]

    # Create all files
    vcf_path = os.path.join(args.output_dir, 'variants.vcf.gz')
    expression_path = os.path.join(args.output_dir, 'expression_counts.csv')
    gtf_path = os.path.join(args.output_dir, 'genes.gtf')
    covariates_path = os.path.join(args.output_dir, 'covariates.csv')

    create_vcf(vcf_path, sample_ids, args.n_variants, args.seed)
    create_expression_matrix(expression_path, sample_ids, gene_ids, args.seed)
    create_gtf(gtf_path, gene_ids, args.seed)
    create_covariates(covariates_path, sample_ids, args.seed)

    logger.info("=" * 50)
    logger.info("Dummy dataset creation complete!")
    logger.info(f"  Samples: {args.n_samples}")
    logger.info(f"  Genes: {args.n_genes}")
    logger.info(f"  Variants: {args.n_variants}")
    logger.info(f"  Output directory: {args.output_dir}")
    logger.info("=" * 50)


if __name__ == '__main__':
    main()
