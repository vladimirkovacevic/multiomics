"""
WGCNA + mQTL multi-omics pipeline.

Workflow
--------
  1. Load expression / genotypes / covariates and align samples.
  2. Run WGCNA (R, in Docker) on the EXPRESSION matrix only.
  3. Treat each module eigengene as a quantitative trait.
  4. mQTL: regress every SNP against every module eigengene
     (Bonferroni + Benjamini-Hochberg FDR).
  5. Optional: correlate module eigengenes with clinical covariates.
  6. Render plots and a self-contained HTML report.

Genotypes are NEVER merged into the WGCNA input. Mixing 0/1/2 dosages
with continuous expression conflates fundamentally different signal types.

CLI
---
  python wgcna_docker.py \
      --expression  data/ASD_dataset/ASD_expression.csv \
      --genotypes   data/ASD_dataset/ASD_genotypes.csv \
      --covariates  data/ASD_dataset/ASD_covariates.csv \
      --output-dir  results/wgcna_run
"""

from __future__ import annotations

import argparse
import base64
import html
import io
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

LOG = logging.getLogger("wgcna")

DEFAULT_TRAITS = ("ASD", "Age", "Sex", "IQ", "ADOS_Score")


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
@dataclass
class WGCNAConfig:
    expression: Path
    genotypes: Path
    covariates: Optional[Path] = None
    output_dir: Path = Path("wgcna_results")
    docker_image: str = "kkhaichau/weighted_networks"
    r_script: Optional[Path] = None
    power: Optional[int] = None
    min_module_size: int = 10
    merge_cut_height: float = 0.25
    network_type: str = "unsigned"
    log_transform: bool = True
    top_var_genes: Optional[int] = None
    trait_columns: tuple = DEFAULT_TRAITS
    fdr_threshold: float = 0.10
    enrichment: bool = False
    enrichment_sets: tuple = (
        "GO_Biological_Process_2023",
        "KEGG_2021_Human",
        "MSigDB_Hallmark_2020",
    )
    enrichment_top: int = 10
    enrichment_min_module_size: int = 5


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def _read_indexed_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    first = df.columns[0]
    if first.lower() in {"sample", "sampleid", "sample_id", "id"}:
        df = df.set_index(first)
    else:
        df = df.set_index(first)
    return df


def load_inputs(cfg: WGCNAConfig):
    LOG.info("Loading inputs")
    expr = _read_indexed_csv(cfg.expression)
    geno = _read_indexed_csv(cfg.genotypes)
    cov = _read_indexed_csv(cfg.covariates) if cfg.covariates and cfg.covariates.exists() else None

    samples = expr.index.intersection(geno.index)
    if cov is not None:
        samples = samples.intersection(cov.index)

    expr = expr.loc[samples]
    geno = geno.loc[samples]
    cov = cov.loc[samples] if cov is not None else None

    poly = geno.var(axis=0) > 0
    if (~poly).any():
        LOG.info("Dropping %d monomorphic SNPs", int((~poly).sum()))
        geno = geno.loc[:, poly]

    LOG.info("Aligned: %d samples | %d genes | %d SNPs%s",
             len(samples), expr.shape[1], geno.shape[1],
             f" | {cov.shape[1]} covariates" if cov is not None else "")
    return expr, geno, cov


def preprocess_expression(expr_raw: pd.DataFrame, cfg: WGCNAConfig) -> pd.DataFrame:
    expr = expr_raw.copy()
    if cfg.log_transform:
        expr = np.log2(expr + 1)
        LOG.info("Applied log2(x + 1) to expression")
    if cfg.top_var_genes and cfg.top_var_genes < expr.shape[1]:
        keep = expr.var(axis=0).sort_values(ascending=False).head(cfg.top_var_genes).index
        expr = expr[keep]
        LOG.info("Filtered to top %d most variable genes", cfg.top_var_genes)
    return expr


# ---------------------------------------------------------------------------
# WGCNA (R inside Docker)
# ---------------------------------------------------------------------------
def _locate_r_script() -> Optional[Path]:
    here = Path(__file__).resolve().parent
    for cand in (here / "wgcna_analysis.R", here.parent / "scripts" / "wgcna_analysis.R"):
        if cand.exists():
            return cand
    return None


def _parse_params(stdout: str) -> dict:
    out = {}
    for line in stdout.splitlines():
        if "PARAMS" in line:
            tail = line.split("PARAMS", 1)[1].strip()
            for token in tail.split(","):
                if "=" in token:
                    k, v = token.split("=", 1)
                    out[k.strip()] = v.strip()
    return out


def run_wgcna_in_docker(expr_clean: pd.DataFrame, cfg: WGCNAConfig):
    r_script = cfg.r_script or _locate_r_script()
    if r_script is None or not r_script.exists():
        raise FileNotFoundError("wgcna_analysis.R not found next to wgcna_docker.py")

    with tempfile.TemporaryDirectory(prefix="wgcna_") as tmp_str:
        tmp = Path(tmp_str)
        in_csv = tmp / "expression.csv"
        modules_csv = tmp / "modules.csv"
        eigen_csv = tmp / "eigengenes.csv"
        sft_csv = tmp / "soft_threshold.csv"
        shutil.copy(r_script, tmp / "wgcna_analysis.R")

        if expr_clean.index.name is None:
            expr_clean.index.name = "SampleID"
        expr_clean.to_csv(in_csv)

        cmd = [
            "docker", "run", "--rm",
            "--platform", "linux/amd64",
            "-v", f"{tmp}:/data",
            cfg.docker_image,
            "Rscript", "/data/wgcna_analysis.R",
            "/data/expression.csv",
            "/data/modules.csv",
            "/data/eigengenes.csv",
            "/data/soft_threshold.csv",
            str(cfg.power) if cfg.power is not None else "NA",
            str(cfg.min_module_size),
            str(cfg.merge_cut_height),
            cfg.network_type,
        ]
        LOG.info("Running WGCNA in Docker (%s)", cfg.docker_image)
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if result.stdout:
            LOG.debug(result.stdout.strip())

        modules_df = pd.read_csv(modules_csv)
        eigengenes_df = pd.read_csv(eigen_csv).set_index("SampleID")
        eigengenes_df.index = expr_clean.index
        sft_df = pd.read_csv(sft_csv)

    params = _parse_params(result.stdout)
    LOG.info("WGCNA done: %s", params)
    return modules_df, eigengenes_df, sft_df, params


# ---------------------------------------------------------------------------
# Downstream: mQTL and module-trait
# ---------------------------------------------------------------------------
def run_mqtl(eigengenes: pd.DataFrame, genotypes: pd.DataFrame) -> pd.DataFrame:
    samples = eigengenes.index.intersection(genotypes.index)
    me, geno = eigengenes.loc[samples], genotypes.loc[samples]

    rows = []
    for module in me.columns:
        y = me[module].values.astype(float)
        for snp in geno.columns:
            x = geno[snp].values.astype(float)
            if np.unique(x).size < 2:
                continue
            slope, intercept, r, p, se = stats.linregress(x, y)
            rows.append({"Module": module, "SNP": snp,
                         "beta": slope, "SE": se, "r": r,
                         "pvalue": p, "n": len(x)})
    df = pd.DataFrame(rows)
    if df.empty:
        return df

    m = len(df)
    df["p_bonferroni"] = (df["pvalue"] * m).clip(upper=1.0)
    df = df.sort_values("pvalue").reset_index(drop=True)
    df["fdr_bh"] = (df["pvalue"] * m / (df.index + 1)).clip(upper=1.0)
    df["fdr_bh"] = df["fdr_bh"][::-1].cummin()[::-1]
    LOG.info("mQTL: %d tests | nominal p<0.05: %d | FDR<0.10: %d",
             m, int((df["pvalue"] < 0.05).sum()), int((df["fdr_bh"] < 0.10).sum()))
    return df


def module_trait_corr(eigengenes: pd.DataFrame,
                      covariates: Optional[pd.DataFrame],
                      traits: tuple):
    if covariates is None:
        return pd.DataFrame(), pd.DataFrame()
    avail = [t for t in traits if t in covariates.columns]
    if not avail:
        LOG.info("No matching trait columns found")
        return pd.DataFrame(), pd.DataFrame()

    r_rows, p_rows = [], []
    for module in eigengenes.columns:
        r_row, p_row = {}, {}
        for trait in avail:
            valid = ~covariates[trait].isna()
            if valid.sum() < 3:
                r_row[trait] = np.nan
                p_row[trait] = np.nan
                continue
            r, p = stats.pearsonr(eigengenes.loc[valid, module],
                                  covariates.loc[valid, trait])
            r_row[trait] = r
            p_row[trait] = p
        r_rows.append(r_row)
        p_rows.append(p_row)
    return (pd.DataFrame(r_rows, index=eigengenes.columns),
            pd.DataFrame(p_rows, index=eigengenes.columns))


# ---------------------------------------------------------------------------
# Module enrichment (Enrichr via gseapy)
# ---------------------------------------------------------------------------
def run_module_enrichment(modules_df: pd.DataFrame,
                          background_genes: list,
                          gene_sets: tuple,
                          top_n: int,
                          min_module_size: int) -> pd.DataFrame:
    """Per-module ORA against Enrichr libraries. Background = WGCNA input genes.

    Returns long-form DataFrame with one row per (Module, gene-set library, term).
    Skips module 0 (grey) and modules smaller than `min_module_size`.
    """
    try:
        import gseapy as gp
    except ImportError:
        LOG.warning("gseapy not installed; skipping --enrichment. "
                    "Install with: pip install gseapy")
        return pd.DataFrame()

    rows = []
    sizes = modules_df["Module"].value_counts()
    for module in sorted(modules_df["Module"].unique()):
        if module == 0:
            continue
        if sizes.get(module, 0) < min_module_size:
            LOG.info("Skipping module %d for enrichment (size %d < %d)",
                     module, sizes.get(module, 0), min_module_size)
            continue
        genes = modules_df.loc[modules_df["Module"] == module, "Gene"].tolist()
        try:
            enr = gp.enrichr(
                gene_list=genes,
                gene_sets=list(gene_sets),
                background=background_genes,
                organism="human",
                outdir=None,
                no_plot=True,
            )
        except Exception as exc:
            LOG.warning("Enrichr failed for module %s: %s", module, exc)
            continue

        res = enr.results
        if res is None or res.empty:
            continue
        for lib, sub in res.groupby("Gene_set"):
            top = sub.sort_values("Adjusted P-value").head(top_n)
            for _, r in top.iterrows():
                genes_str = r.get("Genes", "")
                if "Overlap" in r.index:
                    overlap = r["Overlap"]
                else:
                    n_hit = len(str(genes_str).split(";")) if genes_str else 0
                    overlap = f"{n_hit}/{len(genes)}"
                rows.append({
                    "Module":      f"ME{int(module)}",
                    "Library":     lib,
                    "Term":        r["Term"],
                    "Overlap":     overlap,
                    "P":           r["P-value"],
                    "FDR":         r["Adjusted P-value"],
                    "Odds":        r.get("Odds Ratio", np.nan),
                    "Combined":    r.get("Combined Score", np.nan),
                    "Genes":       genes_str,
                })
    df = pd.DataFrame(rows)
    if not df.empty:
        LOG.info("Enrichment: %d (module, term) hits across %d modules",
                 len(df), df["Module"].nunique())
    return df


# ---------------------------------------------------------------------------
# Plots (each saved to disk and base64-encoded for the HTML report)
# ---------------------------------------------------------------------------
def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


def _save_and_encode(fig, path: Path) -> str:
    fig.savefig(path, dpi=200, bbox_inches="tight")
    b64 = _fig_to_b64(fig)
    return b64


def make_plots(modules_df, eigengenes, sft_df, mqtl_df,
               trait_r, trait_p, genotypes, params, plot_dir: Path) -> dict:
    plot_dir.mkdir(parents=True, exist_ok=True)
    images = {}

    chosen_beta = None
    try:
        chosen_beta = int(params.get("power", ""))
    except (TypeError, ValueError):
        chosen_beta = None

    # --- 1. Soft-threshold scan ----------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    if "Power" in sft_df.columns and "SFT.R.sq" in sft_df.columns:
        ax1.plot(sft_df["Power"], sft_df["SFT.R.sq"], "o-")
        ax1.axhline(0.8, color="red", linestyle="--", label="R² = 0.8")
        if chosen_beta is not None:
            ax1.axvline(chosen_beta, color="green", linestyle="--",
                        label=f"chosen β = {chosen_beta}")
        ax1.set_xlabel("Soft-threshold power β")
        ax1.set_ylabel("Scale-free fit R²")
        ax1.set_title("Scale-free topology")
        ax1.legend()
        ax1.grid(True, alpha=0.3)
    if "Power" in sft_df.columns and "mean.k." in sft_df.columns:
        ax2.plot(sft_df["Power"], sft_df["mean.k."], "o-", color="orange")
        if chosen_beta is not None:
            ax2.axvline(chosen_beta, color="green", linestyle="--",
                        label=f"chosen β = {chosen_beta}")
            ax2.legend()
        ax2.set_xlabel("Soft-threshold power β")
        ax2.set_ylabel("Mean connectivity")
        ax2.set_title("Mean connectivity")
        ax2.grid(True, alpha=0.3)
    fig.tight_layout()
    images["soft_threshold"] = _save_and_encode(fig, plot_dir / "soft_threshold.png")

    # --- 2. Module sizes ------------------------------------------------------
    sizes = modules_df["Module"].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(sizes.index.astype(str), sizes.values, color="steelblue")
    ax.set_xlabel("Module")
    ax.set_ylabel("Number of genes")
    ax.set_title("Module sizes")
    ax.grid(True, alpha=0.3, axis="y")
    fig.tight_layout()
    images["module_sizes"] = _save_and_encode(fig, plot_dir / "module_sizes.png")

    # --- 3. Module eigengenes heatmap ----------------------------------------
    fig, ax = plt.subplots(figsize=(11, max(2.5, 0.4 * eigengenes.shape[1])))
    im = ax.imshow(eigengenes.T.values, aspect="auto", cmap="RdBu_r")
    ax.set_yticks(range(eigengenes.shape[1]))
    ax.set_yticklabels(eigengenes.columns)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Module eigengene")
    ax.set_title("Module eigengenes across samples")
    fig.colorbar(im, ax=ax, label="Eigengene value")
    fig.tight_layout()
    images["eigengenes_heatmap"] = _save_and_encode(fig, plot_dir / "eigengenes_heatmap.png")

    # --- 4. Module-trait correlations ----------------------------------------
    if not trait_r.empty:
        fig, ax = plt.subplots(figsize=(max(5, 1.2 * trait_r.shape[1]),
                                        max(2.5, 0.45 * trait_r.shape[0])))
        im = ax.imshow(trait_r.astype(float).values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
        ax.set_xticks(range(trait_r.shape[1]))
        ax.set_xticklabels(trait_r.columns)
        ax.set_yticks(range(trait_r.shape[0]))
        ax.set_yticklabels(trait_r.index)
        ax.set_title("Module–trait correlations")
        for i in range(trait_r.shape[0]):
            for j in range(trait_r.shape[1]):
                r = trait_r.iat[i, j]
                p = trait_p.iat[i, j]
                if pd.notna(r):
                    stars = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
                    ax.text(j, i, f"{r:.2f}{stars}", ha="center", va="center",
                            fontsize=8, fontweight="bold")
        fig.colorbar(im, ax=ax, label="Pearson r")
        fig.tight_layout()
        images["module_trait"] = _save_and_encode(fig, plot_dir / "module_trait.png")

    # --- 5. mQTL heatmap + Manhattan -----------------------------------------
    if not mqtl_df.empty:
        pivot = mqtl_df.pivot(index="Module", columns="SNP", values="pvalue")
        logp = -np.log10(pivot)
        fig, axes = plt.subplots(1, 2, figsize=(15, 4.5))
        im = axes[0].imshow(logp.values, cmap="viridis", aspect="auto")
        axes[0].set_yticks(range(logp.shape[0]))
        axes[0].set_yticklabels(logp.index)
        axes[0].set_xticks(range(logp.shape[1]))
        axes[0].set_xticklabels(logp.columns, rotation=90, fontsize=6)
        axes[0].set_title("mQTL  -log10(p)")
        fig.colorbar(im, ax=axes[0], label="-log10(p)")

        plot_df = mqtl_df.copy()
        plot_df["xpos"] = np.arange(len(plot_df))
        palette = {m: plt.cm.tab10(i % 10) for i, m in enumerate(plot_df["Module"].unique())}
        for module, sub in plot_df.groupby("Module"):
            axes[1].scatter(sub["xpos"], -np.log10(sub["pvalue"]),
                            s=14, color=palette[module], label=module, alpha=0.75)
        bonf = -np.log10(0.05 / len(plot_df))
        axes[1].axhline(bonf, color="red", linestyle="--", lw=0.8,
                        label=f"Bonferroni 0.05  ({bonf:.1f})")
        axes[1].axhline(-np.log10(0.05), color="gray", linestyle=":", lw=0.8,
                        label="Nominal 0.05")
        axes[1].set_xlabel("mQTL test (sorted by p)")
        axes[1].set_ylabel("-log10(p)")
        axes[1].set_title("mQTL Manhattan-style")
        axes[1].legend(loc="upper right", fontsize=7, ncol=2)
        axes[1].grid(True, alpha=0.3)
        fig.tight_layout()
        images["mqtl_overview"] = _save_and_encode(fig, plot_dir / "mqtl_overview.png")

        # Top hit boxplot
        top = mqtl_df.iloc[0]
        x = genotypes.loc[eigengenes.index, top["SNP"]].values
        y = eigengenes[top["Module"]].values
        groups = [y[x == g] for g in (0, 1, 2) if (x == g).any()]
        labels = [f"{g}\n(n={len(d)})" for g, d in zip((0, 1, 2), groups) if len(d) > 0]
        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        try:
            bp = ax.boxplot(groups, tick_labels=labels, patch_artist=True, widths=0.5)
        except TypeError:
            bp = ax.boxplot(groups, labels=labels, patch_artist=True, widths=0.5)
        for patch in bp["boxes"]:
            patch.set_facecolor("#9ecae1")
        ax.set_xlabel(f"{top['SNP']} dosage")
        ax.set_ylabel(f"{top['Module']} eigengene")
        ax.set_title(f"Top mQTL: {top['Module']} ~ {top['SNP']}  "
                     f"(β={top['beta']:.3f}, p={top['pvalue']:.2g})")
        ax.grid(True, alpha=0.3, axis="y")
        fig.tight_layout()
        images["mqtl_top_hit"] = _save_and_encode(fig, plot_dir / "mqtl_top_hit.png")

    return images


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------
_HTML_CSS = """
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
       margin: 2em auto; max-width: 1100px; color: #222; line-height: 1.45; }
h1 { border-bottom: 2px solid #333; padding-bottom: 0.3em; }
h2 { margin-top: 2em; border-bottom: 1px solid #ccc; padding-bottom: 0.2em; }
h3 { margin-top: 1.5em; color: #444; }
table { border-collapse: collapse; margin: 1em 0; font-size: 0.92em; }
th, td { border: 1px solid #ddd; padding: 4px 10px; text-align: right; }
th { background: #f4f4f4; text-align: left; }
td:first-child, th:first-child { text-align: left; }
.meta td:first-child { font-weight: bold; width: 220px; }
img { max-width: 100%; height: auto; border: 1px solid #eee; padding: 4px; }
.footer { margin-top: 3em; color: #888; font-size: 0.85em; }
.note { background: #fff8d8; border-left: 4px solid #e0c000;
        padding: 0.6em 1em; margin: 1em 0; border-radius: 4px; }
"""


def _df_to_html(df: pd.DataFrame, max_rows: int = 25) -> str:
    if df.empty:
        return "<p><i>(none)</i></p>"
    return df.head(max_rows).to_html(
        index=False, float_format=lambda v: f"{v:.4g}", escape=True,
    )


def _meta_table(rows: list[tuple[str, str]]) -> str:
    body = "".join(f"<tr><td>{html.escape(k)}</td><td>{html.escape(str(v))}</td></tr>" for k, v in rows)
    return f'<table class="meta">{body}</table>'


def render_html_report(cfg, expr, genotypes, modules_df,
                       sft_df, mqtl_df, trait_r, trait_p, enrichment_df,
                       params, images, output_path: Path) -> None:
    n_modules = modules_df["Module"].nunique()
    grey = int((modules_df["Module"] == 0).sum())
    sizes = modules_df["Module"].value_counts().sort_index()
    sizes_html = "<ul>" + "".join(
        f"<li>Module <b>{m}</b>: {n} genes{' (grey/unassigned)' if m == 0 else ''}</li>"
        for m, n in sizes.items()) + "</ul>"

    sig_mqtl = mqtl_df[mqtl_df["fdr_bh"] < cfg.fdr_threshold] if not mqtl_df.empty else pd.DataFrame()
    sig_traits = []
    if not trait_r.empty:
        for mod in trait_r.index:
            for trait in trait_r.columns:
                p = trait_p.at[mod, trait]
                r = trait_r.at[mod, trait]
                if pd.notna(p) and p < 0.05:
                    sig_traits.append((mod, trait, r, p))
    trait_rows = "".join(
        f"<tr><td>{m}</td><td>{t}</td><td>{r:+.3f}</td><td>{p:.4g}</td></tr>"
        for m, t, r, p in sig_traits) or "<tr><td colspan='4'><i>(none at p&lt;0.05)</i></td></tr>"

    def img(key, caption):
        if key not in images:
            return ""
        return (f'<figure><img alt="{caption}" '
                f'src="data:image/png;base64,{images[key]}"/>'
                f'<figcaption><i>{caption}</i></figcaption></figure>')

    # Enrichment block: per-module collapsed sections with top terms
    if enrichment_df is None or enrichment_df.empty:
        enrichment_html = (
            "<p><i>Enrichment not run "
            "(use <code>--enrichment</code>) or no significant terms.</i></p>"
        )
    else:
        cols = ["Library", "Term", "Overlap", "P", "FDR", "Genes"]
        view = enrichment_df[["Module", *cols]].copy()
        view["P"] = view["P"].map(lambda v: f"{v:.2g}")
        view["FDR"] = view["FDR"].map(lambda v: f"{v:.2g}")
        # Truncate long gene-membership strings
        view["Genes"] = view["Genes"].str.slice(0, 80) + view["Genes"].apply(
            lambda s: " ..." if len(s) > 80 else "")
        sections = []
        for mod, sub in view.groupby("Module", sort=True):
            n_fdr = (enrichment_df.loc[enrichment_df["Module"] == mod, "FDR"] < 0.05).sum()
            sections.append(
                f"<details open><summary><b>{html.escape(str(mod))}</b> "
                f"&mdash; top terms across libraries "
                f"(FDR&lt;0.05: {int(n_fdr)})</summary>"
                f"{sub[cols].to_html(index=False, escape=True)}"
                f"</details>"
            )
        enrichment_html = (
            f"<p>Per-module Over-Representation Analysis (Fisher exact, "
            f"Enrichr) against background = the genes used for WGCNA. "
            f"Showing top {cfg.enrichment_top} terms per library, sorted by FDR. "
            f"Modules with &lt;{cfg.enrichment_min_module_size} genes "
            f"and the grey/unassigned module are skipped.</p>"
            + "\n".join(sections)
        )

    default_beta = 6 if cfg.network_type == "unsigned" else 12

    if {"Power", "mean.k.", "max.k."}.issubset(sft_df.columns):
        hub_table = pd.DataFrame({
            "β":        sft_df["Power"].astype(int),
            "mean k":   sft_df["mean.k."].round(3),
            "max k":    sft_df["max.k."].round(3),
            "max/mean": (sft_df["max.k."] / sft_df["mean.k."]).round(2),
        })
    else:
        hub_table = pd.DataFrame()

    meta = _meta_table([
        ("Run timestamp", datetime.now().isoformat(timespec="seconds")),
        ("Expression file", cfg.expression),
        ("Genotypes file", cfg.genotypes),
        ("Covariates file", cfg.covariates or "—"),
        ("Samples (aligned)", expr.shape[0]),
        ("Genes used for WGCNA", expr.shape[1]),
        ("SNPs used for mQTL", genotypes.shape[1]),
        ("Docker image", cfg.docker_image),
        ("Soft-threshold β", params.get("power", "?")),
        ("Network type", params.get("network", cfg.network_type)),
        ("min module size", params.get("min_module_size", cfg.min_module_size)),
        ("merge cut height", params.get("merge_cut_height", cfg.merge_cut_height)),
        ("Modules (incl. grey)", n_modules),
        ("Genes in grey/unassigned", grey),
    ])

    body = f"""
<h1>WGCNA + mQTL Multi-Omics Report</h1>
<div class="note">
  WGCNA was run on the <b>expression matrix only</b>. Module eigengenes are
  treated as quantitative traits, then tested against every SNP in the
  genotype matrix (mQTL). This avoids merging discrete genotypes with
  continuous expression in a single similarity network.
</div>

<h2>1. Run summary</h2>
{meta}

<h2>2. Soft-threshold β selection</h2>
<p>WGCNA <code>pickSoftThreshold</code> scans β = 1..20. For each β it builds
adjacency a<sub>ij</sub> = |cor(i,j)|<sup>β</sup>, computes node connectivity
k, fits log P(k) ~ log k and reports the scale-free fit R². The
<b>smallest β with R² ≥ 0.8</b> is selected (and mean connectivity should not
be too low). If no β reaches that threshold, a default of {default_beta} is
used ({cfg.network_type} network).</p>
<p>Selected for this run: <b>β = {params.get("power", "?")}</b>
&nbsp;|&nbsp; network = <b>{params.get("network", cfg.network_type)}</b></p>
{img("soft_threshold", "Scale-free fit R² and mean connectivity vs β. Green dashed line = chosen β.")}

<h3>Hub-concentration table</h3>
<p><code>max k / mean k</code> is a rough hub-concentration index: higher β
crushes weak correlations to 0 while strong ones survive, so a few genes
accumulate disproportionate connectivity. Hubs emerge as this ratio grows.</p>
{_df_to_html(hub_table, max_rows=20)}

<details><summary>Full scale-free fit table (per β)</summary>
{_df_to_html(sft_df.round(4), max_rows=20)}
</details>

<h2>3. Modules</h2>
<p>WGCNA identified <b>{n_modules}</b> modules (including grey =
unassigned). Module sizes:</p>
{sizes_html}
{img("module_sizes", "Genes per module.")}
{img("eigengenes_heatmap", "Module eigengenes (PC1) across samples.")}

<h2>4. Module–trait correlations</h2>
<p>Significant (p &lt; 0.05) module–trait associations:</p>
<table><thead><tr><th>Module</th><th>Trait</th><th>r</th><th>p</th></tr></thead>
<tbody>{trait_rows}</tbody></table>
{img("module_trait", "Pearson correlation between module eigengenes and clinical traits.")}

<h2>5. mQTL — module eigengenes as quantitative traits</h2>
<p>Total tests: <b>{len(mqtl_df)}</b>; nominal p&lt;0.05:
<b>{int((mqtl_df['pvalue'] < 0.05).sum()) if not mqtl_df.empty else 0}</b>;
FDR&lt;{cfg.fdr_threshold:g}: <b>{len(sig_mqtl)}</b>;
Bonferroni&lt;0.05: <b>{int((mqtl_df['p_bonferroni'] < 0.05).sum()) if not mqtl_df.empty else 0}</b>.</p>

{img("mqtl_overview", "Left: -log10(p) heatmap (modules × SNPs). Right: Manhattan-style scatter.")}

<h3>Top mQTL associations</h3>
{_df_to_html(mqtl_df[['Module','SNP','beta','SE','pvalue','fdr_bh','p_bonferroni']],
             max_rows=20)}

{img("mqtl_top_hit", "Strongest single mQTL: eigengene distribution by genotype.")}

<h2>6. Module enrichment (ORA via Enrichr)</h2>
{enrichment_html}

<h2>7. Outputs on disk</h2>
<ul>
  <li><code>modules.csv</code> — gene → module assignment</li>
  <li><code>eigengenes.csv</code> — sample × module eigengene matrix</li>
  <li><code>soft_threshold.csv</code> — pickSoftThreshold fit indices</li>
  <li><code>mqtl_results.csv</code> — full mQTL table (sorted by p)</li>
  <li><code>module_trait_correlations.csv</code> — module × trait r and p</li>
  <li><code>module_enrichment.csv</code> — top GO/KEGG/Hallmark terms per module (if --enrichment)</li>
  <li><code>plots/</code> — high-resolution PNGs</li>
</ul>

<div class="footer">
  Generated by <code>scripts/wgcna_docker.py</code>.
</div>
"""

    output_path.write_text(
        f"<!doctype html><html><head><meta charset='utf-8'>"
        f"<title>WGCNA + mQTL Report</title><style>{_HTML_CSS}</style></head>"
        f"<body>{body}</body></html>"
    )
    LOG.info("Wrote HTML report: %s", output_path)


# ---------------------------------------------------------------------------
# CLI / orchestration
# ---------------------------------------------------------------------------
def _setup_logging(level: str) -> None:
    logging.basicConfig(
        level=level.upper(),
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%H:%M:%S",
    )


def _build_config(args) -> WGCNAConfig:
    return WGCNAConfig(
        expression=Path(args.expression),
        genotypes=Path(args.genotypes),
        covariates=Path(args.covariates) if args.covariates else None,
        output_dir=Path(args.output_dir),
        docker_image=args.docker_image,
        power=args.power,
        min_module_size=args.min_module_size,
        merge_cut_height=args.merge_cut_height,
        network_type=args.network_type,
        log_transform=not args.no_log_transform,
        top_var_genes=args.top_var_genes,
        trait_columns=tuple(args.traits) if args.traits else DEFAULT_TRAITS,
        fdr_threshold=args.fdr_threshold,
        enrichment=args.enrichment,
        enrichment_sets=tuple(args.enrichment_sets) if args.enrichment_sets
                        else WGCNAConfig.enrichment_sets,
        enrichment_top=args.enrichment_top,
        enrichment_min_module_size=args.enrichment_min_module_size,
    )


def _parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="WGCNA on expression (in Docker) + mQTL on module eigengenes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--expression", required=True, help="samples × genes CSV")
    p.add_argument("--genotypes",  required=True, help="samples × SNPs CSV (0/1/2)")
    p.add_argument("--covariates", default=None,  help="samples × traits CSV (optional)")
    p.add_argument("--output-dir", default="wgcna_results", help="output directory")
    p.add_argument("--docker-image", default="kkhaichau/weighted_networks")
    p.add_argument("--power", type=int, default=None,
                   help="soft-threshold β; auto-picked when omitted")
    p.add_argument("--min-module-size", type=int, default=10)
    p.add_argument("--merge-cut-height", type=float, default=0.25)
    p.add_argument("--network-type", choices=("unsigned", "signed"), default="unsigned")
    p.add_argument("--top-var-genes", type=int, default=None,
                   help="restrict to top-N most variable genes before WGCNA")
    p.add_argument("--no-log-transform", action="store_true",
                   help="skip log2(x+1) (data already on log scale)")
    p.add_argument("--traits", nargs="*", default=None,
                   help=f"trait columns from covariates (default: {' '.join(DEFAULT_TRAITS)})")
    p.add_argument("--fdr-threshold", type=float, default=0.10)
    p.add_argument("--enrichment", action="store_true",
                   help="run per-module ORA via gseapy/Enrichr (needs internet)")
    p.add_argument("--enrichment-sets", nargs="*", default=None,
                   help="Enrichr libraries to query "
                        "(default: GO_Biological_Process_2023, "
                        "KEGG_2021_Human, MSigDB_Hallmark_2020)")
    p.add_argument("--enrichment-top", type=int, default=10,
                   help="top N terms to keep per library per module")
    p.add_argument("--enrichment-min-module-size", type=int, default=5,
                   help="skip enrichment for modules smaller than this")
    p.add_argument("--log-level", default="INFO",
                   choices=("DEBUG", "INFO", "WARNING", "ERROR"))
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = _parse_args(argv)
    _setup_logging(args.log_level)
    cfg = _build_config(args)
    cfg.output_dir.mkdir(parents=True, exist_ok=True)

    expr_raw, geno, cov = load_inputs(cfg)
    expr = preprocess_expression(expr_raw, cfg)

    modules_df, eigengenes, sft_df, params = run_wgcna_in_docker(expr, cfg)
    mqtl_df = run_mqtl(eigengenes, geno)
    trait_r, trait_p = module_trait_corr(eigengenes, cov, cfg.trait_columns)

    # Persist tables
    modules_df.to_csv(cfg.output_dir / "modules.csv", index=False)
    eigengenes.to_csv(cfg.output_dir / "eigengenes.csv")
    sft_df.to_csv(cfg.output_dir / "soft_threshold.csv", index=False)
    mqtl_df.to_csv(cfg.output_dir / "mqtl_results.csv", index=False)
    if not trait_r.empty:
        merged = trait_r.add_suffix("_r").join(trait_p.add_suffix("_p"))
        merged.to_csv(cfg.output_dir / "module_trait_correlations.csv")

    enrichment_df = pd.DataFrame()
    if cfg.enrichment:
        LOG.info("Running module enrichment (Enrichr) — needs internet")
        enrichment_df = run_module_enrichment(
            modules_df,
            background_genes=expr.columns.tolist(),
            gene_sets=cfg.enrichment_sets,
            top_n=cfg.enrichment_top,
            min_module_size=cfg.enrichment_min_module_size,
        )
        if not enrichment_df.empty:
            enrichment_df.to_csv(cfg.output_dir / "module_enrichment.csv", index=False)

    images = make_plots(modules_df, eigengenes, sft_df, mqtl_df,
                        trait_r, trait_p, geno, params,
                        plot_dir=cfg.output_dir / "plots")

    render_html_report(cfg, expr, geno, modules_df,
                       sft_df, mqtl_df, trait_r, trait_p, enrichment_df,
                       params, images,
                       output_path=cfg.output_dir / "report.html")

    LOG.info("Done. Open %s", cfg.output_dir / "report.html")
    return 0


if __name__ == "__main__":
    sys.exit(main())
