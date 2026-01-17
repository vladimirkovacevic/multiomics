"""
WGCNA Analysis using R implementation via Docker
This module provides a Python wrapper for running WGCNA analysis
using the R WGCNA package inside a Docker container.

Supports multi-omics data: expression + genotypes merged together.
"""

import pandas as pd
import numpy as np
import subprocess
import tempfile
import os
import shutil
from pathlib import Path
from sklearn.preprocessing import StandardScaler


def preprocess_expression(expression_df):
    """
    Preprocess expression data for WGCNA.

    Steps:
    1. Log2(x + 1) transformation to handle count data skewness
    2. Z-score standardization (mean=0, std=1 per gene)

    Parameters:
    -----------
    expression_df : pd.DataFrame
        Raw expression counts (samples x genes)

    Returns:
    --------
    pd.DataFrame : Preprocessed expression data
    """
    # Log-transform to handle skewed count distribution
    expr_log = np.log2(expression_df + 1)

    # Z-score standardize
    scaler = StandardScaler()
    expr_scaled = pd.DataFrame(
        scaler.fit_transform(expr_log),
        index=expression_df.index,
        columns=expression_df.columns
    )

    return expr_scaled


def preprocess_genotypes(genotypes_df):
    """
    Preprocess genotype data for WGCNA.

    Genotypes are coded as 0, 1, 2 (number of minor alleles).
    Apply z-score standardization to put on same scale as expression.

    Parameters:
    -----------
    genotypes_df : pd.DataFrame
        Genotype data (samples x SNPs), values 0/1/2

    Returns:
    --------
    pd.DataFrame : Preprocessed genotype data
    """
    scaler = StandardScaler()
    geno_scaled = pd.DataFrame(
        scaler.fit_transform(genotypes_df),
        index=genotypes_df.index,
        columns=genotypes_df.columns
    )

    return geno_scaled


def merge_multiomics(expression_scaled, genotypes_scaled):
    """
    Merge expression and genotype data into combined feature matrix.

    Adds prefixes to distinguish feature types:
    - 'expr_' for expression features
    - 'geno_' for genotype features

    Parameters:
    -----------
    expression_scaled : pd.DataFrame
        Preprocessed expression data
    genotypes_scaled : pd.DataFrame
        Preprocessed genotype data

    Returns:
    --------
    pd.DataFrame : Combined feature matrix (samples x features)
    """
    # Add prefixes to distinguish feature types
    expr_renamed = expression_scaled.copy()
    expr_renamed.columns = ['expr_' + str(col) for col in expr_renamed.columns]

    geno_renamed = genotypes_scaled.copy()
    geno_renamed.columns = ['geno_' + str(col) for col in geno_renamed.columns]

    # Merge horizontally
    combined = pd.concat([expr_renamed, geno_renamed], axis=1)

    return combined


def run_wgcna_docker(
    expression_df,
    genotypes_df=None,
    n_modules=4,
    docker_image="kkhaichau/weighted_networks",
    preprocess=True
):
    """
    Run WGCNA analysis using R implementation in Docker container.

    Supports both single-omics (expression only) and multi-omics
    (expression + genotypes) analysis.

    Parameters:
    -----------
    expression_df : pd.DataFrame
        Expression data with samples as rows and genes as columns.
        Index should be sample IDs, columns should be gene names.
    genotypes_df : pd.DataFrame, optional
        Genotype data with samples as rows and SNPs as columns.
        Values should be 0, 1, or 2 (number of minor alleles).
        If None, only expression data is used.
    n_modules : int, default=4
        Number of co-expression modules to identify.
    docker_image : str, default="kkhaichau/weighted_networks"
        Docker image containing R and WGCNA package.
    preprocess : bool, default=True
        Whether to apply preprocessing (log-transform + z-score for expression,
        z-score for genotypes). Set to False if data is already preprocessed.

    Returns:
    --------
    module_assignment : pd.DataFrame
        DataFrame with columns ['Feature', 'Module', 'Type'] mapping features to modules.
    me_df : pd.DataFrame
        DataFrame with module eigengenes (samples x modules).
        Columns are named 'ME0', 'ME1', etc.
    combined_data : pd.DataFrame
        The preprocessed combined data matrix used for WGCNA.
    """
    print("=" * 80)
    print("WGCNA Multi-Omics Analysis")
    print("=" * 80)

    # Preprocess data
    if preprocess:
        print("\nPreprocessing data...")
        print(f"  Expression: {expression_df.shape[0]} samples x {expression_df.shape[1]} genes")
        expr_processed = preprocess_expression(expression_df)
        print(f"    Applied log2(x+1) + z-score standardization")

        if genotypes_df is not None:
            print(f"  Genotypes: {genotypes_df.shape[0]} samples x {genotypes_df.shape[1]} SNPs")
            geno_processed = preprocess_genotypes(genotypes_df)
            print(f"    Applied z-score standardization")
        else:
            geno_processed = None
    else:
        expr_processed = expression_df
        geno_processed = genotypes_df

    # Merge data if genotypes provided
    if geno_processed is not None:
        combined_data = merge_multiomics(expr_processed, geno_processed)
        print(f"\nMerged multi-omics data: {combined_data.shape[0]} samples x {combined_data.shape[1]} features")
        n_expr = len([c for c in combined_data.columns if c.startswith('expr_')])
        n_geno = len([c for c in combined_data.columns if c.startswith('geno_')])
        print(f"  Expression features: {n_expr}")
        print(f"  Genotype features: {n_geno}")
    else:
        combined_data = expr_processed.copy()
        combined_data.columns = ['expr_' + str(col) for col in combined_data.columns]
        print(f"\nUsing expression data only: {combined_data.shape}")

    # Create temporary directory for data exchange
    temp_dir = tempfile.mkdtemp(prefix="wgcna_")

    try:
        # Define file paths
        input_csv = os.path.join(temp_dir, "expression_data.csv")
        module_csv = os.path.join(temp_dir, "module_assignment.csv")
        eigengene_csv = os.path.join(temp_dir, "module_eigengenes.csv")
        r_script = os.path.join(temp_dir, "wgcna_analysis.R")

        # Copy R script to temp directory
        script_dir = os.path.dirname(os.path.abspath(__file__))
        # Try multiple possible locations
        possible_paths = [
            os.path.join(script_dir, "wgcna_analysis.R"),
            os.path.join(script_dir, "..", "wgcna_analysis.R"),
            os.path.join(script_dir, "..", "scripts", "wgcna_analysis.R"),
        ]
        script_source = None
        for path in possible_paths:
            if os.path.exists(path):
                script_source = path
                break
        if script_source is None:
            raise FileNotFoundError(f"Could not find wgcna_analysis.R script. Searched: {possible_paths}")
        shutil.copy(script_source, r_script)

        # Save combined data to CSV
        if combined_data.index.name is None:
            combined_data.index.name = "SampleID"
        combined_data.to_csv(input_csv)

        # Build docker command
        docker_cmd = [
            "docker", "run", "--rm",
            "--platform", "linux/amd64",
            "-v", f"{temp_dir}:/data",
            docker_image,
            "Rscript", "/data/wgcna_analysis.R",
            "/data/expression_data.csv",
            "/data/module_assignment.csv",
            "/data/module_eigengenes.csv",
            str(n_modules)
        ]

        # Run Docker container
        print(f"\nRunning WGCNA in Docker container...")
        result = subprocess.run(
            docker_cmd,
            capture_output=True,
            text=True,
            check=True
        )

        # Print R script output
        print(result.stdout)

        if result.stderr:
            print("R Messages:", result.stderr)

        # Read results
        module_assignment_raw = pd.read_csv(module_csv)
        eigengenes_raw = pd.read_csv(eigengene_csv)

        # Process module assignment - add feature type
        module_assignment = module_assignment_raw.rename(columns={'Gene': 'Feature'})
        module_assignment['Type'] = module_assignment['Feature'].apply(
            lambda x: 'Expression' if x.startswith('expr_') else 'Genotype'
        )

        # Process eigengenes
        me_df = eigengenes_raw.set_index('SampleID')
        me_df.index = combined_data.index

        # Validate results
        if len(module_assignment) != len(combined_data.columns):
            raise ValueError(
                f"Module assignment count ({len(module_assignment)}) "
                f"doesn't match number of features ({len(combined_data.columns)})"
            )

        if me_df.shape[0] != combined_data.shape[0]:
            raise ValueError(
                f"Module eigengene samples ({me_df.shape[0]}) "
                f"doesn't match data samples ({combined_data.shape[0]})"
            )

        # Print summary
        print("\n" + "=" * 80)
        print("WGCNA COMPLETE")
        print("=" * 80)
        print(f"\nModules identified: {module_assignment['Module'].nunique()}")
        for mod in sorted(module_assignment['Module'].unique()):
            subset = module_assignment[module_assignment['Module'] == mod]
            n_expr = (subset['Type'] == 'Expression').sum()
            n_geno = (subset['Type'] == 'Genotype').sum()
            print(f"  Module {mod}: {len(subset)} features ({n_expr} expr, {n_geno} geno)")

        return module_assignment, me_df, combined_data

    except subprocess.CalledProcessError as e:
        print(f"Error running Docker container:")
        print(f"Return code: {e.returncode}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        raise

    finally:
        # Clean up temporary directory
        shutil.rmtree(temp_dir, ignore_errors=True)


def run_wgcna_python(combined_data, n_modules=5):
    """
    Run WGCNA analysis using pure Python implementation.

    This is a fallback when Docker is not available.
    Uses scipy for clustering and sklearn for PCA.

    Parameters:
    -----------
    combined_data : pd.DataFrame
        Preprocessed feature matrix (samples x features)
    n_modules : int, default=5
        Number of modules to identify

    Returns:
    --------
    module_assignment : pd.DataFrame
        DataFrame with columns ['Feature', 'Module', 'Type']
    me_df : pd.DataFrame
        Module eigengenes (samples x modules)
    """
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LinearRegression

    print("\nRunning Python-based WGCNA implementation...")

    # Step 1: Compute correlation matrix
    print("  Computing correlation matrix...")
    S = combined_data.corr(method='pearson')

    # Step 2: Select soft-thresholding power
    print("  Selecting soft-thresholding power...")
    powers = range(1, 21)
    best_power = 6  # Default

    for beta in powers:
        A = np.abs(S.values) ** beta
        k = A.sum(axis=1)
        k_unique, k_counts = np.unique(k[k > 0], return_counts=True)
        if len(k_unique) >= 3:
            log_k = np.log10(k_unique).reshape(-1, 1)
            log_pk = np.log10(k_counts)
            model = LinearRegression().fit(log_k, log_pk)
            r_sq = model.score(log_k, log_pk)
            if r_sq >= 0.8 and beta >= 4:
                best_power = beta
                break

    print(f"  Selected power: {best_power}")

    # Step 3: Build adjacency matrix
    print("  Building adjacency matrix...")
    A = np.abs(S.values) ** best_power

    # Step 4: Compute TOM
    print("  Computing TOM...")
    n = A.shape[0]
    k = A.sum(axis=1)
    L = A @ A

    TOM = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                TOM[i, j] = 1.0
            else:
                numerator = L[i, j] + A[i, j]
                denominator = min(k[i], k[j]) + 1 - A[i, j]
                TOM[i, j] = numerator / denominator if denominator > 0 else 0

    dissTOM = 1 - TOM

    # Step 5: Hierarchical clustering
    print("  Clustering features...")
    dissTOM_condensed = squareform(dissTOM, checks=False)
    tree = linkage(dissTOM_condensed, method='average')

    # Cut tree to get target number of modules
    labels = fcluster(tree, t=n_modules, criterion='maxclust')

    # Step 6: Create module assignment
    module_colors = ['grey', 'turquoise', 'blue', 'brown', 'yellow', 'green',
                     'red', 'black', 'pink', 'magenta', 'purple', 'orange']

    module_assignment = pd.DataFrame({
        'Feature': combined_data.columns.tolist(),
        'Module': labels,
        'ModuleColor': [module_colors[i % len(module_colors)] for i in labels],
        'Type': ['Expression' if c.startswith('expr_') else 'Genotype'
                 for c in combined_data.columns]
    })

    # Step 7: Compute module eigengenes
    print("  Computing module eigengenes...")
    ME_dict = {}
    for mod in sorted(module_assignment['Module'].unique()):
        features = module_assignment[module_assignment['Module'] == mod]['Feature'].tolist()
        module_data = combined_data[features]

        if len(features) == 1:
            ME_dict[f'ME{mod}'] = module_data.iloc[:, 0].values
        else:
            pca = PCA(n_components=1)
            ME_dict[f'ME{mod}'] = pca.fit_transform(module_data).flatten()

    me_df = pd.DataFrame(ME_dict, index=combined_data.index)

    return module_assignment, me_df


def load_asd_dataset(data_dir="../data/ASD_dataset/"):
    """
    Load the ASD dataset for testing.

    Parameters:
    -----------
    data_dir : str
        Path to the ASD dataset directory.

    Returns:
    --------
    expression : pd.DataFrame
        Expression data (samples x genes)
    genotypes : pd.DataFrame
        Genotype data (samples x SNPs)
    covariates : pd.DataFrame
        Clinical covariates
    """
    expression = pd.read_csv(f"{data_dir}/ASD_expression.csv").set_index('sample')
    genotypes = pd.read_csv(f"{data_dir}/ASD_genotypes.csv").set_index('sample')
    covariates = pd.read_csv(f"{data_dir}/ASD_covariates.csv").set_index('sample')

    print(f"Loaded ASD dataset:")
    print(f"  Expression: {expression.shape}")
    print(f"  Genotypes: {genotypes.shape}")
    print(f"  Covariates: {covariates.shape}")

    return expression, genotypes, covariates


# For backwards compatibility
def wgcna_analysis(expression_indexed, n_modules=4):
    """
    Convenience wrapper for run_wgcna_docker with the same interface
    as the original sklearn-based implementation.
    """
    module_assignment, me_df, _ = run_wgcna_docker(
        expression_indexed,
        genotypes_df=None,
        n_modules=n_modules
    )
    return module_assignment, me_df


def run_wgcna_multiomics(
    expression_df,
    genotypes_df=None,
    n_modules=5,
    preprocess=True,
    use_docker=True,
    docker_image="kkhaichau/weighted_networks"
):
    """
    Main entry point for multi-omics WGCNA analysis.

    Tries Docker first, falls back to Python implementation if Docker unavailable.

    Parameters:
    -----------
    expression_df : pd.DataFrame
        Expression data (samples x genes)
    genotypes_df : pd.DataFrame, optional
        Genotype data (samples x SNPs), values 0/1/2
    n_modules : int, default=5
        Number of modules to identify
    preprocess : bool, default=True
        Apply preprocessing (log2+z-score for expression, z-score for genotypes)
    use_docker : bool, default=True
        Try Docker first; if False or Docker fails, use Python implementation
    docker_image : str
        Docker image for R WGCNA

    Returns:
    --------
    module_assignment : pd.DataFrame
        Feature to module mapping
    me_df : pd.DataFrame
        Module eigengenes
    combined_data : pd.DataFrame
        Preprocessed combined data matrix
    """
    print("=" * 80)
    print("WGCNA Multi-Omics Analysis")
    print("=" * 80)

    # Preprocess data
    if preprocess:
        print("\nPreprocessing data...")
        print(f"  Expression: {expression_df.shape[0]} samples x {expression_df.shape[1]} genes")
        expr_processed = preprocess_expression(expression_df)
        print(f"    Applied log2(x+1) + z-score standardization")

        if genotypes_df is not None:
            print(f"  Genotypes: {genotypes_df.shape[0]} samples x {genotypes_df.shape[1]} SNPs")
            geno_processed = preprocess_genotypes(genotypes_df)
            print(f"    Applied z-score standardization")
        else:
            geno_processed = None
    else:
        expr_processed = expression_df
        geno_processed = genotypes_df

    # Merge data
    if geno_processed is not None:
        combined_data = merge_multiomics(expr_processed, geno_processed)
        print(f"\nMerged multi-omics data: {combined_data.shape[0]} samples x {combined_data.shape[1]} features")
    else:
        combined_data = expr_processed.copy()
        combined_data.columns = ['expr_' + str(col) for col in combined_data.columns]
        print(f"\nUsing expression data only: {combined_data.shape}")

    # Try Docker first if requested
    if use_docker:
        try:
            print("\nAttempting Docker-based R WGCNA...")
            return run_wgcna_docker(
                expression_df, genotypes_df, n_modules,
                docker_image=docker_image, preprocess=preprocess
            )
        except Exception as e:
            print(f"\nDocker unavailable: {e}")
            print("Falling back to Python implementation...")

    # Python fallback
    module_assignment, me_df = run_wgcna_python(combined_data, n_modules)

    # Print summary
    print("\n" + "=" * 80)
    print("WGCNA COMPLETE")
    print("=" * 80)
    print(f"\nModules identified: {module_assignment['Module'].nunique()}")
    for mod in sorted(module_assignment['Module'].unique()):
        subset = module_assignment[module_assignment['Module'] == mod]
        n_expr = (subset['Type'] == 'Expression').sum()
        n_geno = (subset['Type'] == 'Genotype').sum()
        print(f"  Module {mod}: {len(subset)} features ({n_expr} expr, {n_geno} geno)")

    return module_assignment, me_df, combined_data


def main():
    """Command line interface for WGCNA multi-omics analysis."""
    import sys
    import argparse

    # Determine default paths (ASD dataset)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_data_dir = os.path.join(script_dir, "..", "data", "ASD_dataset")
    default_expression = os.path.join(default_data_dir, "ASD_expression.csv")
    default_genotypes = os.path.join(default_data_dir, "ASD_genotypes.csv")
    default_covariates = os.path.join(default_data_dir, "ASD_covariates.csv")

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="WGCNA Multi-Omics Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use default ASD dataset
  python wgcna_docker.py

  # Specify custom expression file only
  python wgcna_docker.py --expression my_expression.csv

  # Specify both expression and genotypes
  python wgcna_docker.py --expression expr.csv --genotypes geno.csv

  # Specify all files with custom parameters
  python wgcna_docker.py -e expr.csv -g geno.csv -c cov.csv -n 6 --no-docker
        """
    )

    parser.add_argument(
        "-e", "--expression",
        default=default_expression,
        help=f"Expression data CSV (samples x genes). Default: ASD_expression.csv"
    )
    parser.add_argument(
        "-g", "--genotypes",
        default=default_genotypes,
        help=f"Genotypes data CSV (samples x SNPs, values 0/1/2). Default: ASD_genotypes.csv"
    )
    parser.add_argument(
        "-c", "--covariates",
        default=default_covariates,
        help=f"Covariates CSV for trait correlation. Default: ASD_covariates.csv"
    )
    parser.add_argument(
        "-n", "--n-modules",
        type=int,
        default=5,
        help="Number of modules to identify. Default: 5"
    )
    parser.add_argument(
        "--no-genotypes",
        action="store_true",
        help="Run with expression data only (no genotypes)"
    )
    parser.add_argument(
        "--no-docker",
        action="store_true",
        help="Skip Docker, use Python implementation directly"
    )
    parser.add_argument(
        "--no-preprocess",
        action="store_true",
        help="Skip preprocessing (data already log-transformed and scaled)"
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output prefix for results files (module_assignment.csv, eigengenes.csv)"
    )

    args = parser.parse_args()

    print("WGCNA Multi-Omics Analysis")
    print("=" * 80)

    try:
        # Load expression data
        if not os.path.exists(args.expression):
            print(f"Error: Expression file not found: {args.expression}")
            sys.exit(1)

        print(f"\nLoading expression data: {args.expression}")
        expression = pd.read_csv(args.expression)
        # Set first column as index if it looks like sample IDs
        if expression.columns[0].lower() in ['sample', 'sampleid', 'sample_id', 'id']:
            expression = expression.set_index(expression.columns[0])
        print(f"  Shape: {expression.shape}")

        # Load genotypes data (optional)
        genotypes = None
        if not args.no_genotypes:
            if os.path.exists(args.genotypes):
                print(f"\nLoading genotypes data: {args.genotypes}")
                genotypes = pd.read_csv(args.genotypes)
                if genotypes.columns[0].lower() in ['sample', 'sampleid', 'sample_id', 'id']:
                    genotypes = genotypes.set_index(genotypes.columns[0])
                print(f"  Shape: {genotypes.shape}")
            else:
                print(f"\nGenotypes file not found: {args.genotypes}")
                print("  Running with expression data only")

        # Load covariates (optional, for trait correlation)
        covariates = None
        if os.path.exists(args.covariates):
            print(f"\nLoading covariates: {args.covariates}")
            covariates = pd.read_csv(args.covariates)
            if covariates.columns[0].lower() in ['sample', 'sampleid', 'sample_id', 'id']:
                covariates = covariates.set_index(covariates.columns[0])
            print(f"  Shape: {covariates.shape}")

        # Run WGCNA
        print(f"\nRunning WGCNA analysis (n_modules={args.n_modules})...")
        module_assignment, me_df, combined_data = run_wgcna_multiomics(
            expression,
            genotypes_df=genotypes,
            n_modules=args.n_modules,
            preprocess=not args.no_preprocess,
            use_docker=not args.no_docker
        )

        # Save results if output prefix specified
        if args.output:
            module_file = f"{args.output}_modules.csv"
            eigengene_file = f"{args.output}_eigengenes.csv"
            module_assignment.to_csv(module_file, index=False)
            me_df.to_csv(eigengene_file)
            print(f"\nResults saved:")
            print(f"  Module assignment: {module_file}")
            print(f"  Module eigengenes: {eigengene_file}")

        # Show module-trait correlations if covariates available
        if covariates is not None:
            print("\n" + "=" * 80)
            print("Module-Trait Correlations")
            print("=" * 80)

            from scipy import stats

            # Select numeric columns for correlation
            numeric_cols = covariates.select_dtypes(include=[np.number]).columns.tolist()
            # Prioritize common trait columns
            priority_cols = ['ASD', 'Age', 'Sex', 'IQ', 'ADOS_Score']
            trait_cols = [c for c in priority_cols if c in numeric_cols]
            # Add other numeric columns (up to 10 total)
            for c in numeric_cols:
                if c not in trait_cols and len(trait_cols) < 10:
                    trait_cols.append(c)

            if trait_cols:
                traits = covariates[trait_cols].copy()

                for module in me_df.columns:
                    print(f"\n{module}:")
                    for trait in traits.columns:
                        valid = ~traits[trait].isna()
                        if valid.sum() > 0:
                            r, p = stats.pearsonr(me_df.loc[valid, module], traits.loc[valid, trait])
                            sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
                            print(f"  {trait}: r={r:.3f}, p={p:.4f} {sig}")
            else:
                print("No numeric columns found in covariates for correlation analysis")

        print("\n" + "=" * 80)
        print("Analysis complete!")
        print("=" * 80)

    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
