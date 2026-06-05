"""Stage runners that glue modules into a config-driven pipeline.

Each stage reads from `<output_dir>/<prev_step>/...` and writes to its own
sub-directory, so stages can be re-run independently.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from tqdm import tqdm

from .cluster import cluster_projection_matrix
from .config import Config, Sample, load_sample_sheet, source_samples, target_samples
from .errorcorrect import (
    global_collapse,
    identify_spike_ccs,
    local_collapse_sample,
    pool_for_global,
)
from .filterplot import plot_threshold_scan, scan_thresholds
from .matrix import (
    build_barcode_matrix,
    filter_by_thresholds,
    filter_homopolymers,
    normalize_by_spikes,
    role_indices,
    save_matrix,
    split_spike_ins_by_cc,
    target_sort_indices,
)
from .preprocess import preprocess_sample, umi_collapse, write_quickout
from .rankplot import plot_overview_grid, plot_sample_rank, read_bu_counts, suggest_umi_threshold
from .viz import plot_projection_heatmap, plot_ssi_check

from collections import Counter

log = logging.getLogger(__name__)


def _stage_dir(cfg: Config, name: str) -> Path:
    p = cfg.project.output_dir / name
    p.mkdir(parents=True, exist_ok=True)
    return p


def run_preprocess(cfg: Config) -> None:
    samples = load_sample_sheet(cfg.input.sample_sheet)
    out = _stage_dir(cfg, "01_preprocess")
    for s in tqdm(samples, desc="preprocess"):
        result = preprocess_sample(s, cfg)
        bu_path = out / f"{s.sample_id}_BU.tsv"
        with bu_path.open("w") as fh:
            for bu, c in result.bu_counts.most_common():
                fh.write(f"{c}\t{bu}\n")
        write_quickout(
            umi_collapse(result.bu_counts, cfg.barcode.barcode_length, s.umi_threshold),
            out / f"{s.sample_id}_quickout.tsv",
        )
        ssi_path = out / f"{s.sample_id}_ssi_check.tsv"
        with ssi_path.open("w") as fh:
            for ssi, c in result.ssi_check.most_common():
                fh.write(f"{c}\t{ssi}\n")
        stats_path = out / f"{s.sample_id}_stats.tsv"
        with stats_path.open("w") as fh:
            fh.write("metric\tvalue\n")
            fh.write(f"reads_seen\t{result.reads_seen}\n")
            fh.write(f"reads_after_yy\t{result.reads_after_yy}\n")
            fh.write(f"reads_after_ssi\t{result.reads_after_ssi}\n")
            fh.write(f"unique_bu\t{len(result.bu_counts)}\n")
        plot_path = out / f"{s.sample_id}_ssi_check.png"
        try:
            plot_ssi_check(
                dict(result.ssi_check),
                s.sample_id,
                s.ssi_match(cfg.barcode.ssi_match_length),
                plot_path,
            )
        except Exception as e:
            log.warning("ssi-check plot failed for %s: %s", s.sample_id, e)


def _read_bu_tsv(path: Path) -> Counter[str]:
    """Reload {sample}_BU.tsv (count\\tsequence) into a Counter for re-collapse."""
    counter: Counter[str] = Counter()
    if not path.exists():
        return counter
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) < 2:
                continue
            counter[parts[1]] = int(parts[0])
    return counter


def run_umi_collapse(cfg: Config) -> None:
    """Re-run UMI collapse from existing 01_preprocess/{sample}_BU.tsv.

    Use this after editing sample_sheet.csv `umi_threshold` based on rank plots,
    so you don't re-read FASTQ files.
    """
    samples = load_sample_sheet(cfg.input.sample_sheet)
    out = _stage_dir(cfg, "01_preprocess")
    for s in tqdm(samples, desc="umi-collapse"):
        bu_path = out / f"{s.sample_id}_BU.tsv"
        bu_counts = _read_bu_tsv(bu_path)
        if not bu_counts:
            log.warning("no BU counts for sample %s, skipping", s.sample_id)
            continue
        write_quickout(
            umi_collapse(bu_counts, cfg.barcode.barcode_length, s.umi_threshold),
            out / f"{s.sample_id}_quickout.tsv",
        )


def run_rankplot(
    cfg: Config,
    log_x: bool = False,
    xlim: int | None = None,
    annotate: bool = False,
    overview: bool = False,
    suggest: bool = False,
) -> None:
    """Per-sample BU rank plots (port of rankplot2.m).

    Defaults match the MATLAB original: one semilogy figure per sample, simple
    title, no annotation, no overview. Enable extras via flags.

    suggest=True overlays the heuristic Hyopil-style threshold (shoulder // 2)
    and writes a `suggested_thresholds.tsv` next to the plots.
    """
    samples = load_sample_sheet(cfg.input.sample_sheet)
    pre_dir = cfg.project.output_dir / "01_preprocess"
    out = _stage_dir(cfg, "01_preprocess") / "_rankplots"
    out.mkdir(parents=True, exist_ok=True)

    samples_counts: list[tuple] = []
    suggestions: list[tuple[str, int, int, int]] = []
    for s in tqdm(samples, desc="rankplot"):
        bu_path = pre_dir / f"{s.sample_id}_BU.tsv"
        counts = read_bu_counts(bu_path)
        samples_counts.append((s, counts))
        if suggest:
            thr, sh, ei = suggest_umi_threshold(counts)
            suggestions.append((s.sample_id, thr, sh, ei))
        plot_sample_rank(
            counts, s, out / f"sample_{s.sample_id}_rank.png",
            log_x=log_x, xlim=xlim, annotate=annotate, suggest=suggest,
        )
    if overview:
        plot_overview_grid(
            samples_counts, out / "overview.png",
            log_x=log_x, xlim=xlim, annotate=annotate, suggest=suggest,
        )
    if suggest:
        sugg_path = out / "suggested_thresholds.tsv"
        with sugg_path.open("w") as fh:
            fh.write("sample_id\tsuggested_threshold\tshoulder_count\telbow_rank\n")
            for sid, thr, sh, ei in suggestions:
                fh.write(f"{sid}\t{thr}\t{sh}\t{ei}\n")
        log.info("suggested thresholds: %s", sugg_path)
    log.info("rank plots written to %s", out)


def run_local_collapse(cfg: Config) -> None:
    samples = load_sample_sheet(cfg.input.sample_sheet)
    pre_dir = cfg.project.output_dir / "01_preprocess"
    out = _stage_dir(cfg, "02_local_collapse")
    for s in tqdm(samples, desc="local collapse"):
        quickout = pre_dir / f"{s.sample_id}_quickout.tsv"
        if not quickout.exists():
            log.warning("missing quickout for %s, skipping", s.sample_id)
            continue
        local_collapse_sample(s.sample_id, quickout, out, cfg)


def run_global_collapse(cfg: Config) -> None:
    samples = load_sample_sheet(cfg.input.sample_sheet)
    local_dir = cfg.project.output_dir / "02_local_collapse"
    out = _stage_dir(cfg, "03_global_collapse")
    sample_collapsed_paths = {
        s.sample_id: local_dir / f"{s.sample_id}_local_collapsed.tsv" for s in samples
    }
    pooled_seqs, per_sample_counts = pool_for_global(sample_collapsed_paths)
    log.info("pooled %d unique sequences across %d samples", len(pooled_seqs), len(samples))
    global_labels, _ = global_collapse(pooled_seqs, out, cfg)

    (out / "pooled_seqs.tsv").write_text(
        "pooled_idx\tsequence\n" + "\n".join(f"{i}\t{s}" for i, s in enumerate(pooled_seqs)) + "\n"
    )
    np.save(out / "global_labels.npy", global_labels)
    psc_path = out / "per_sample_counts.tsv"
    with psc_path.open("w") as fh:
        fh.write("sample_id\tsequence\tcount\n")
        for sid, sc in per_sample_counts.items():
            for seq, c in sc.items():
                fh.write(f"{sid}\t{seq}\t{c}\n")


def _load_global_artifacts(cfg: Config) -> tuple[list[str], np.ndarray, dict[str, dict[str, int]]]:
    glob_dir = cfg.project.output_dir / "03_global_collapse"
    pooled_seqs: list[str] = []
    with (glob_dir / "pooled_seqs.tsv").open() as fh:
        next(fh)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                pooled_seqs.append(parts[1])
    global_labels = np.load(glob_dir / "global_labels.npy")
    per_sample_counts: dict[str, dict[str, int]] = {}
    with (glob_dir / "per_sample_counts.tsv").open() as fh:
        next(fh)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            sid, seq, c = parts
            per_sample_counts.setdefault(sid, {})[seq] = int(c)
    return pooled_seqs, global_labels, per_sample_counts


def run_matrix(cfg: Config) -> None:
    """Build the per-CC barcode matrix (port of produceBCmat3).

    Pipeline:
      1. aggregate per-sample counts into [n_global_cc x n_samples]
      2. tag each global CC as spike (any member ends with spike tag) or real
      3. drop CCs whose representative has a homopolymer run >= min_run
      4. split into real_matrix / spike_matrix using CC-level flag
      5. spike-normalize real_matrix (divide each column by spike column total)
    """
    samples = load_sample_sheet(cfg.input.sample_sheet)
    pooled_seqs, global_labels, per_sample_counts = _load_global_artifacts(cfg)
    out = _stage_dir(cfg, "04_matrix")

    matrix, refbarcodes, sample_ids = build_barcode_matrix(
        pooled_seqs, global_labels, per_sample_counts, samples
    )
    sample_labels = [f"{s.sample_id}:{s.region_label}:{s.role}" for s in samples]
    save_matrix(matrix, refbarcodes, sample_labels, out / "barcode_matrix.tsv")

    is_spike_cc = identify_spike_ccs(pooled_seqs, global_labels, cfg.spike_in_tag)
    n_spike_pre = int(is_spike_cc.sum())

    hp_min = cfg.quality_filter.homopolymer_min_run
    matrix, refbarcodes, hp_keep, (is_spike_cc,) = filter_homopolymers(
        matrix, refbarcodes, hp_min, extra=(is_spike_cc,)
    )
    n_dropped_hp = int((~hp_keep).sum())

    real_mat, real_bcs, spike_mat, spike_bcs = split_spike_ins_by_cc(
        matrix, refbarcodes, is_spike_cc
    )
    save_matrix(real_mat, real_bcs, sample_labels, out / "real_matrix.tsv")
    save_matrix(spike_mat, spike_bcs, sample_labels, out / "spike_matrix.tsv")

    norm_mat = normalize_by_spikes(real_mat, spike_mat, method=cfg.normalization.spike_method)
    save_matrix(norm_mat, real_bcs, sample_labels, out / "normalized_matrix.tsv")

    log.info(
        "matrix: %d CCs (spike CCs: %d, dropped by homopolymer: %d) → real %d, spike %d, %d samples",
        matrix.shape[0] + n_dropped_hp, n_spike_pre, n_dropped_hp,
        real_mat.shape[0], spike_mat.shape[0], matrix.shape[1],
    )


def run_filterplot(
    cfg: Config,
    max_thresh: int = 100,
    max_targets: int | None = None,
    scan: str = "source",
    fixed_proj: int = 0,
    fixed_source: int = 0,
    show_totals: bool = False,
) -> None:
    """Source/projection-threshold scan (port of numproj_injthresh.m).

    Reads 04_matrix/real_matrix.tsv, sweeps threshold values, writes:
      04_matrix/_filterplots/{scan}_threshold_scan.{tsv,png}

    Inspect the plot, then edit config.yaml `filtering.source_threshold_umi`
    or `proj_threshold_umi`, then run `mapseq2 filter`.
    """
    samples = load_sample_sheet(cfg.input.sample_sheet)
    mat_dir = cfg.project.output_dir / "04_matrix"
    out = mat_dir / "_filterplots"
    out.mkdir(parents=True, exist_ok=True)

    from .matrix import load_matrix

    real_mat, _real_bcs, _sample_labels = load_matrix(mat_dir / "real_matrix.tsv")
    source_idx = role_indices(samples, "source")
    target_idx = target_sort_indices(samples)
    if max_targets is None:
        max_targets = max(1, len(target_idx))

    if scan == "source":
        title = (
            f"source-UMI threshold scan "
            f"(proj_threshold fixed = {fixed_proj}, "
            f"{len(source_idx)} source / {len(target_idx)} target columns)"
        )
        scan_df = scan_thresholds(
            real_mat, source_idx, target_idx,
            max_thresh=max_thresh, max_targets=max_targets,
            fixed_proj_threshold=fixed_proj, scan="source",
        )
    elif scan == "proj":
        title = (
            f"projection-UMI threshold scan "
            f"(source_threshold fixed = {fixed_source}, "
            f"{len(source_idx)} source / {len(target_idx)} target columns)"
        )
        scan_df = scan_thresholds(
            real_mat, source_idx, target_idx,
            max_thresh=max_thresh, max_targets=max_targets,
            fixed_source_threshold=fixed_source, scan="proj",
        )
    else:
        raise ValueError(f"scan must be 'source' or 'proj', got {scan!r}")

    tsv_path = out / f"{scan}_threshold_scan.tsv"
    png_path = out / f"{scan}_threshold_scan.png"
    scan_df.to_csv(tsv_path, sep="\t", index=False)
    plot_threshold_scan(
        scan_df, png_path, max_targets=max_targets,
        scan_label=scan, show_totals=show_totals,
        title=title if show_totals else None,
    )
    log.info("filterplot (%s) → %s", scan, png_path)


def run_filter(cfg: Config) -> None:
    samples = load_sample_sheet(cfg.input.sample_sheet)
    mat_dir = cfg.project.output_dir / "04_matrix"
    out = _stage_dir(cfg, "05_filtered")

    from .matrix import load_matrix

    real_mat, real_bcs, sample_labels = load_matrix(mat_dir / "real_matrix.tsv")
    spike_mat, _, _ = load_matrix(mat_dir / "spike_matrix.tsv")

    source_idx = role_indices(samples, "source")
    target_idx = target_sort_indices(samples)

    B, Bseq, keep = filter_by_thresholds(
        real_mat, real_bcs,
        source_indices=source_idx,
        target_indices=target_idx,
        source_threshold=cfg.filtering.source_threshold_umi,
        proj_threshold=cfg.filtering.proj_threshold_umi,
        cell_body_threshold=cfg.filtering.cell_body_threshold,
    )
    save_matrix(B, Bseq, sample_labels, out / "B.tsv")

    Bnorm = normalize_by_spikes(B, spike_mat, method=cfg.normalization.spike_method)
    save_matrix(Bnorm, Bseq, sample_labels, out / "Bnorm.tsv")

    target_labels = [sample_labels[i] for i in target_idx]
    B_tar = B[:, target_idx]
    Bnorm_tar = Bnorm[:, target_idx]
    save_matrix(B_tar, Bseq, target_labels, out / "B_tar.tsv")
    save_matrix(Bnorm_tar, Bseq, target_labels, out / "Bnorm_tar.tsv")

    log.info(
        "filter: %d/%d barcodes kept (source>%g, proj>%g)",
        keep.sum(), len(real_bcs),
        cfg.filtering.source_threshold_umi, cfg.filtering.proj_threshold_umi,
    )


def run_cluster(cfg: Config) -> None:
    samples = load_sample_sheet(cfg.input.sample_sheet)
    filt_dir = cfg.project.output_dir / "05_filtered"
    out = _stage_dir(cfg, "06_cluster")

    from .matrix import load_matrix

    Bnorm_tar, bcseqs, target_labels = load_matrix(filt_dir / "Bnorm_tar.tsv")
    if Bnorm_tar.shape[0] < cfg.clustering.num_clusters:
        log.warning(
            "only %d barcodes survived filtering, fewer than num_clusters=%d; skipping clustering",
            Bnorm_tar.shape[0], cfg.clustering.num_clusters,
        )
        return

    res = cluster_projection_matrix(
        Bnorm_tar,
        num_clusters=cfg.clustering.num_clusters,
        method=cfg.clustering.method,
        metric=cfg.clustering.metric,
        log_factor=cfg.normalization.log_factor,
    )

    np.save(out / "linkage_Z.npy", res.linkage_Z)
    np.save(out / "perm.npy", res.perm)
    np.save(out / "labels.npy", res.labels)
    np.save(out / "sorted_data.npy", res.sorted_data)

    import pandas as pd

    pd.DataFrame({"barcode": bcseqs, "cluster": res.labels}).to_csv(
        out / "cluster_assignments.tsv", sep="\t", index=False,
    )

    # display region labels for heatmap
    target_idx = target_sort_indices(samples)
    region_display = [
        f"{samples[i].region}-{samples[i].side}" if samples[i].side else samples[i].region
        for i in target_idx
    ]

    plot_projection_heatmap(
        res.sorted_data,
        res.sorted_labels,
        res.linkage_Z,
        res.perm,
        region_display,
        out / "heatmap.png",
        cluster_id_to_name=cfg.clustering.cluster_labels,
    )

    log.info(
        "clustering: %d neurons → %d clusters",
        Bnorm_tar.shape[0], cfg.clustering.num_clusters,
    )


def run_all(cfg: Config) -> None:
    run_preprocess(cfg)
    run_local_collapse(cfg)
    run_global_collapse(cfg)
    run_matrix(cfg)
    run_filter(cfg)
    run_cluster(cfg)
