"""click CLI: `mapseq2 <stage> --config config.yaml`."""

from __future__ import annotations

import logging

import click

from .config import Config
from .pipeline import (
    run_all,
    run_cluster,
    run_filter,
    run_filterplot,
    run_global_collapse,
    run_local_collapse,
    run_matrix,
    run_preprocess,
    run_rankplot,
    run_umi_collapse,
)


def _setup_logging(verbose: bool):
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )


def _config_option(f):
    return click.option(
        "-c", "--config", "config_path",
        required=True, type=click.Path(exists=True, dir_okay=False),
        help="Path to config.yaml",
    )(f)


@click.group()
@click.option("-v", "--verbose", is_flag=True)
@click.pass_context
def main(ctx, verbose):
    """mapseq2-py: config-driven MAPseq2 pipeline."""
    _setup_logging(verbose)
    ctx.ensure_object(dict)


@main.command()
@_config_option
@click.option("--engine", type=click.Choice(["bash", "python"]), default="bash",
              show_default=True,
              help="bash = pigz/awk/sort multi-core pipeline (10x faster); "
                   "python = pure-Python single-thread (portable, no shell tools needed)")
def preprocess(config_path, engine):
    """Stage 01: FASTQ → BU counts + SSI QC + UMI collapse."""
    cfg = Config.from_yaml(config_path)
    if engine == "bash":
        from .preprocess_bash import run_preprocess_bash
        run_preprocess_bash(cfg)
    else:
        run_preprocess(cfg)


@main.command()
@_config_option
@click.option("--log-x", is_flag=True, help="loglog (matches MATLAB shape=1)")
@click.option("--xlim", type=int, default=None, help="x-axis upper limit (MATLAB xlimit)")
@click.option("--annotate", is_flag=True,
              help="overlay current umi_threshold + kept_rank lines (extension, not in original)")
@click.option("--overview", is_flag=True,
              help="also emit overview.png grid of all samples (extension, not in original)")
@click.option("--suggest", is_flag=True,
              help="overlay auto-detected threshold (shoulder//2, min 2) + dump suggested_thresholds.tsv")
def rankplot(config_path, log_x, xlim, annotate, overview, suggest):
    """Per-sample BU rank plots (port of rankplot2.m).

    Default matches MATLAB: per-sample semilogy figure, simple title.
    Reads 01_preprocess/{sample}_BU.tsv produced by `preprocess`, writes
    01_preprocess/_rankplots/sample_*.png. Inspect plots, edit sample_sheet.csv
    `umi_threshold` column, then re-run `umi-collapse`.
    """
    run_rankplot(
        Config.from_yaml(config_path),
        log_x=log_x, xlim=xlim,
        annotate=annotate, overview=overview, suggest=suggest,
    )


@main.command("umi-collapse")
@_config_option
def umi_collapse_cmd(config_path):
    """Re-do UMI collapse from existing BU.tsv using current sample_sheet thresholds.

    Faster than re-running `preprocess` after a threshold edit; no FASTQ re-read.
    """
    run_umi_collapse(Config.from_yaml(config_path))


@main.command("collapse-local")
@_config_option
def collapse_local(config_path):
    """Stage 02: per-sample bowtie2 + CC."""
    run_local_collapse(Config.from_yaml(config_path))


@main.command("collapse-global")
@_config_option
def collapse_global(config_path):
    """Stage 03: pool across samples + global bowtie2 + CC."""
    run_global_collapse(Config.from_yaml(config_path))


@main.command()
@_config_option
def matrix(config_path):
    """Stage 04: barcode matrix + spike split + spike normalization."""
    run_matrix(Config.from_yaml(config_path))


@main.command()
@_config_option
@click.option("--max-thresh", default=100, type=int, show_default=True,
              help="Max threshold to scan (MATLAB xlimit; .mlx uses 100).")
@click.option("--max-targets", default=None, type=int,
              help="Cap on n_targets bin (MATLAB n; .mlx uses 19; default = number of target columns).")
@click.option("--scan", default="source",
              type=click.Choice(["source", "proj"]), show_default=True,
              help="Which threshold to scan. MATLAB original scans source only.")
@click.option("--fixed-proj", default=0, type=int, show_default=True,
              help="Fixed projection threshold during source-scan (MATLAB .mlx uses 0).")
@click.option("--fixed-source", default=0, type=int, show_default=True,
              help="Fixed source threshold during projection-scan (MATLAB .mlx uses 0).")
@click.option("--show-totals", is_flag=True,
              help="Add lower n_kept panel (extension, not in original).")
def filterplot(config_path, max_thresh, max_targets, scan, fixed_proj, fixed_source, show_totals):
    """Threshold scan (port of numproj_injthresh.m).

    Default matches MATLAB: stacked-bar of percent-by-n_targets, ylim [0,1],
    fixed_proj=0 (i.e. no projection pre-filter while scanning source).
    Reads 04_matrix/real_matrix.tsv. Writes
    04_matrix/_filterplots/{scan}_threshold_scan.{tsv,png}. Inspect plot,
    edit filtering.* in config.yaml, then run `mapseq2 filter`.
    """
    run_filterplot(
        Config.from_yaml(config_path),
        max_thresh=max_thresh,
        max_targets=max_targets,
        scan=scan,
        fixed_proj=fixed_proj,
        fixed_source=fixed_source,
        show_totals=show_totals,
    )


@main.command()
@_config_option
def filter(config_path):
    """Stage 05: threshold filtering."""
    run_filter(Config.from_yaml(config_path))


@main.command()
@_config_option
def cluster(config_path):
    """Stage 06: Ward clustering + heatmap."""
    run_cluster(Config.from_yaml(config_path))


@main.command()
@_config_option
def run(config_path):
    """Run all stages 01 → 06."""
    run_all(Config.from_yaml(config_path))


if __name__ == "__main__":
    main()
