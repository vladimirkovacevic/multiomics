"""
DIABLO Analysis using mixOmics R implementation via Docker

This module provides a Python wrapper for running DIABLO (Data Integration Analysis
for Biomarker discovery using Latent variable approaches for Omics studies) analysis
using the mixOmics R package inside a Docker container.

DIABLO performs supervised multi-omics integration to identify biomarker signatures
that discriminate between phenotypic groups (e.g., ASD vs Control).
"""

import pandas as pd
import numpy as np
import subprocess
import tempfile
import os
import shutil
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def run_diablo_docker(
    expression_indexed,
    genotypes_indexed,
    covariates_indexed,
    n_components=2,
    keepX_expr=5,
    keepX_geno=5,
    docker_image="vladimirkovacevic/diablo:latest",
    output_prefix="DIABLO"
):
    """
    Run DIABLO analysis using mixOmics R implementation in Docker container.

    DIABLO (Data Integration Analysis for Biomarker discovery using Latent variable
    approaches for Omics studies) is a supervised multi-block integration method that
    identifies multi-omics signatures discriminating between phenotypic groups.

    Parameters
    ----------
    expression_indexed : pd.DataFrame
        Expression data with samples as rows and genes as columns.
        Index should be sample IDs, columns should be gene names.
        Will be log2(x+1) transformed if values > 100 (assumes raw counts).

    genotypes_indexed : pd.DataFrame
        Genotype data with samples as rows and SNPs as columns.
        Index should be sample IDs, columns should be SNP names.
        Values should be dosage encoding (0, 1, 2).

    covariates_indexed : pd.DataFrame
        Covariates data with samples as rows.
        Must contain an 'ASD' column with values 0 (Control) or 1 (ASD).
        Index should be sample IDs matching expression and genotypes.

    n_components : int, default=2
        Number of latent components to extract.
        Typically 2-3 for interpretability.

    keepX_expr : int, default=5
        Number of genes to select per component.
        Controls sparsity of feature selection for expression data.

    keepX_geno : int, default=5
        Number of SNPs to select per component.
        Controls sparsity of feature selection for genotype data.

    docker_image : str, default="vladimirkovacevic/diablo:latest"
        Docker image containing R and mixOmics package.

    output_prefix : str, default="DIABLO"
        Prefix for output files (not including directory path).

    Returns
    -------
    results : dict
        Dictionary containing analysis results:
        - 'components': pd.DataFrame
            Sample component scores (samples × components).
            Consensus representation averaging Expression and Genotypes blocks.

        - 'loadings_expr': pd.DataFrame
            Gene loadings (selected_genes × components).
            Shows contribution of each gene to components.

        - 'loadings_geno': pd.DataFrame
            SNP loadings (selected_SNPs × components).
            Shows contribution of each SNP to components.

        - 'selected_features': pd.DataFrame
            List of selected features with columns:
            ['Component', 'Block', 'Feature', 'Loading']

        - 'performance': pd.DataFrame
            Cross-validation performance metrics:
            ['Component', 'Overall_Error', 'Control_Error', 'ASD_Error', 'BER']
            BER = Balanced Error Rate

        - 'correlations': pd.DataFrame
            Correlation between Expression and Genotypes blocks per component.

    Examples
    --------
    >>> # Load data
    >>> expression = pd.read_csv('expression.csv', index_col='sample')
    >>> genotypes = pd.read_csv('genotypes.csv', index_col='sample')
    >>> covariates = pd.read_csv('covariates.csv', index_col='sample')
    >>>
    >>> # Run DIABLO analysis
    >>> results = run_diablo_docker(
    ...     expression, genotypes, covariates,
    ...     n_components=2, keepX_expr=10, keepX_geno=15
    ... )
    >>>
    >>> # Examine selected features
    >>> print(results['selected_features'])
    >>>
    >>> # Check performance
    >>> print(results['performance'])
    >>>
    >>> # Get component scores for plotting
    >>> components = results['components']

    Notes
    -----
    - All input DataFrames must have matching sample indices
    - Sample order will be preserved in output
    - Docker must be installed and running
    - The docker image will be pulled automatically if not present
    - Analysis uses 5-fold cross-validation with 10 repeats
    - Feature selection uses sparse PLS-DA (block.splsda)
    - Design matrix uses 0.1 correlation between blocks (weak coupling)

    References
    ----------
    Singh A. et al. (2019). DIABLO: an integrative approach for identifying
    key molecular drivers from multi-omics assays. Bioinformatics, 35(17):3055-3062.
    """

    logger.info("="*80)
    logger.info("DIABLO Multi-Omics Integration Analysis")
    logger.info("="*80)
    logger.info("")

    # Validate inputs
    logger.info("Validating input data...")

    if not isinstance(expression_indexed, pd.DataFrame):
        raise TypeError("expression_indexed must be a pandas DataFrame")
    if not isinstance(genotypes_indexed, pd.DataFrame):
        raise TypeError("genotypes_indexed must be a pandas DataFrame")
    if not isinstance(covariates_indexed, pd.DataFrame):
        raise TypeError("covariates_indexed must be a pandas DataFrame")

    # Check sample alignment
    if not expression_indexed.index.equals(genotypes_indexed.index):
        raise ValueError("Expression and genotypes sample indices do not match")
    if not expression_indexed.index.equals(covariates_indexed.index):
        raise ValueError("Expression and covariates sample indices do not match")

    # Check for ASD column
    if 'ASD' not in covariates_indexed.columns:
        raise ValueError("Covariates must contain 'ASD' column (0=Control, 1=ASD)")

    # Validate ASD values
    asd_values = covariates_indexed['ASD'].unique()
    if not set(asd_values).issubset({0, 1}):
        raise ValueError("ASD column must contain only 0 (Control) or 1 (ASD)")

    n_asd = (covariates_indexed['ASD'] == 1).sum()
    n_control = (covariates_indexed['ASD'] == 0).sum()

    logger.info(f"  Expression: {expression_indexed.shape[0]} samples × {expression_indexed.shape[1]} genes")
    logger.info(f"  Genotypes:  {genotypes_indexed.shape[0]} samples × {genotypes_indexed.shape[1]} SNPs")
    logger.info(f"  Outcome:    {n_asd} ASD cases, {n_control} Controls")
    logger.info(f"  Components: {n_components}")
    logger.info(f"  KeepX:      {keepX_expr} genes, {keepX_geno} SNPs per component")
    logger.info("")

    # Create temporary directory for data exchange
    temp_dir = tempfile.mkdtemp(prefix="diablo_")
    logger.info(f"Created temporary directory: {temp_dir}")

    try:
        # Define file paths
        expression_csv = os.path.join(temp_dir, "expression.csv")
        genotypes_csv = os.path.join(temp_dir, "genotypes.csv")
        covariates_csv = os.path.join(temp_dir, "covariates.csv")
        output_temp_prefix = os.path.join(temp_dir, output_prefix)

        # Save data to CSV files
        logger.info("Saving input data to temporary files...")

        # Ensure index name is set for proper CSV export
        if expression_indexed.index.name is None:
            expression_indexed.index.name = "sample"
        if genotypes_indexed.index.name is None:
            genotypes_indexed.index.name = "sample"
        if covariates_indexed.index.name is None:
            covariates_indexed.index.name = "sample"

        expression_indexed.to_csv(expression_csv)
        genotypes_indexed.to_csv(genotypes_csv)
        covariates_indexed.to_csv(covariates_csv)

        # Build Docker command
        docker_cmd = [
            "docker", "run", "--rm",
            "-v", f"{temp_dir}:/data",    # Mount temp directory
            docker_image,
            "Rscript", "/scripts/diablo_analysis.R",
            "/data/expression.csv",
            "/data/genotypes.csv",
            "/data/covariates.csv",
            f"/data/{output_prefix}",
            str(n_components),
            str(keepX_expr),
            str(keepX_geno)
        ]

        logger.info("Running DIABLO analysis in Docker container...")
        logger.info(f"Docker image: {docker_image}")
        logger.info("")

        # Run Docker container
        result = subprocess.run(
            docker_cmd,
            capture_output=True,
            text=True,
            check=True
        )

        # Print R script output
        logger.info("R Script Output:")
        logger.info("-" * 80)
        print(result.stdout)
        logger.info("-" * 80)

        if result.stderr:
            # R messages/warnings (not necessarily errors)
            logger.warning("R Messages:")
            logger.warning(result.stderr)

        logger.info("")
        logger.info("Reading results from Docker container...")

        # Read results
        components = pd.read_csv(f"{output_temp_prefix}_components.csv")
        loadings_expr = pd.read_csv(f"{output_temp_prefix}_loadings_expression.csv")
        loadings_geno = pd.read_csv(f"{output_temp_prefix}_loadings_genotypes.csv")
        selected_features = pd.read_csv(f"{output_temp_prefix}_selected_features.csv")
        performance = pd.read_csv(f"{output_temp_prefix}_performance.csv")
        correlations = pd.read_csv(f"{output_temp_prefix}_correlations.csv")

        # Set sample index for components
        components = components.set_index('SampleID')
        components.index = expression_indexed.index  # Restore original index

        # Validate results
        if len(components) != len(expression_indexed):
            raise ValueError(
                f"Component count ({len(components)}) doesn't match "
                f"input samples ({len(expression_indexed)})"
            )

        logger.info(f"  ✓ Components: {components.shape}")
        logger.info(f"  ✓ Expression loadings: {loadings_expr.shape}")
        logger.info(f"  ✓ Genotype loadings: {loadings_geno.shape}")
        logger.info(f"  ✓ Selected features: {len(selected_features)} total")
        logger.info("")

        # Package results
        results = {
            'components': components,
            'loadings_expr': loadings_expr,
            'loadings_geno': loadings_geno,
            'selected_features': selected_features,
            'performance': performance,
            'correlations': correlations
        }

        logger.info("="*80)
        logger.info("DIABLO Analysis Completed Successfully")
        logger.info("="*80)
        logger.info("")

        return results

    except subprocess.CalledProcessError as e:
        logger.error("Docker container execution failed!")
        logger.error(f"Return code: {e.returncode}")
        logger.error(f"STDOUT: {e.stdout}")
        logger.error(f"STDERR: {e.stderr}")
        raise

    except FileNotFoundError as e:
        logger.error(f"Output file not found: {e}")
        logger.error("The R script may have failed to generate all outputs")
        raise

    finally:
        # Clean up temporary directory
        shutil.rmtree(temp_dir, ignore_errors=True)
        logger.info(f"Cleaned up temporary directory")


# For backwards compatibility and easy import
def diablo_analysis(expression_indexed, genotypes_indexed, covariates_indexed,
                   n_components=2, keepX_expr=5, keepX_geno=5):
    """
    Convenience wrapper for run_diablo_docker.

    See run_diablo_docker() for full documentation.
    """
    return run_diablo_docker(
        expression_indexed, genotypes_indexed, covariates_indexed,
        n_components=n_components, keepX_expr=keepX_expr, keepX_geno=keepX_geno
    )


if __name__ == "__main__":
    # Example usage and help
    print("=" * 80)
    print("DIABLO Docker Wrapper")
    print("=" * 80)
    print("")
    print("This module provides a Python interface to the mixOmics DIABLO implementation.")
    print("")
    print("Usage:")
    print("  from diablo_docker import run_diablo_docker")
    print("")
    print("  results = run_diablo_docker(")
    print("      expression_indexed,  # pd.DataFrame (samples × genes)")
    print("      genotypes_indexed,   # pd.DataFrame (samples × SNPs)")
    print("      covariates_indexed,  # pd.DataFrame with 'ASD' column")
    print("      n_components=2,      # Number of components")
    print("      keepX_expr=5,        # Genes to select per component")
    print("      keepX_geno=5         # SNPs to select per component")
    print("  )")
    print("")
    print("Inputs:")
    print("  - expression_indexed: pandas DataFrame (samples × genes)")
    print("  - genotypes_indexed: pandas DataFrame (samples × SNPs)")
    print("  - covariates_indexed: pandas DataFrame with 'ASD' column (0/1)")
    print("  - n_components: int, number of latent components")
    print("  - keepX_expr: int, number of genes to select per component")
    print("  - keepX_geno: int, number of SNPs to select per component")
    print("")
    print("Outputs (dict):")
    print("  - components: DataFrame with sample component scores")
    print("  - loadings_expr: DataFrame with gene loadings")
    print("  - loadings_geno: DataFrame with SNP loadings")
    print("  - selected_features: DataFrame with selected features")
    print("  - performance: DataFrame with cross-validation metrics")
    print("  - correlations: DataFrame with block correlations")
    print("")
    print("=" * 80)
