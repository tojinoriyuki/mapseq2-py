"""Per-sample rank plots of BU read counts (port of rankplot2.m).

Default behavior matches the original MATLAB function: one semilogy figure
per sample, title = sample id, no threshold annotation, no overview.

Opt-in extras (off by default to stay faithful to the original):
  - annotate=True  : overlay current sample-sheet threshold + kept_rank line
  - overview=True  : also emit a grid figure of all samples

Original `rankplot2(bcn, prefix, shape, xlimit)` arguments map as:
  bcn      → samples list (driven by sample sheet)
  prefix   → cfg.project.name (used in directory layout, not titles)
  shape    → log_x flag (shape=0 semilogy default; shape=1 loglog → log_x=True)
  xlimit   → xlim
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .config import Sample


def read_bu_counts(path: Path) -> np.ndarray:
    counts: list[int] = []
    if not path.exists():
        return np.array([], dtype=np.int64)
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if not parts:
                continue
            try:
                counts.append(int(parts[0]))
            except ValueError:
                continue
    return np.array(counts, dtype=np.int64)


def _annotate_threshold(ax, counts: np.ndarray, threshold: int) -> None:
    if threshold <= 0 or counts.size == 0:
        return
    below = np.where(counts < threshold)[0]
    rank = int(below[0]) if below.size > 0 else int(counts.size)
    ax.axhline(threshold, color="tab:red", linewidth=1, linestyle="--", alpha=0.7)
    ax.axvline(rank, color="tab:red", linewidth=1, linestyle=":", alpha=0.7)
    ax.text(
        0.98, 0.95, f"thresh={threshold}\nkept_rank={rank}",
        transform=ax.transAxes, ha="right", va="top",
        fontsize=8, color="tab:red",
        bbox=dict(facecolor="white", edgecolor="tab:red", alpha=0.8, boxstyle="round,pad=0.3"),
    )


def suggest_umi_threshold(
    counts: np.ndarray,
    min_threshold: int = 2,
    n_log_points: int = 100,
) -> tuple[int, int, int]:
    """Detect the shoulder of a rank plot (Kneedle algorithm, Satopaa 2011) and
    propose a UMI threshold = max(2, shoulder_count // 2).

    The shoulder = top corner of the plateau-to-drop transition (= "L" corner),
    NOT the midpoint of the steepest slope. Found as the point on the log-log
    curve at maximum perpendicular distance from the chord connecting the
    first and last sampled points.

    Method:
      1. Sample n_log_points ranks log-uniformly.
      2. Normalize both axes to [0, 1] (y inverted so plateau = 1).
      3. Chord from (0, 1) → (1, 0) is `y = 1 - x` in normalized space.
      4. For a convex decreasing curve, the elbow lies above this chord; pick
         the point with maximum `y_norm - (1 - x_norm)`.
      5. shoulder_count = counts[elbow_rank]; threshold = max(2, shoulder // 2).
    """
    counts = np.asarray(counts, dtype=float)
    n = len(counts)
    if n < 10:
        return (min_threshold, int(counts.max()) if n else 0, 0)

    rank_idx = np.unique(np.round(np.logspace(0, np.log10(n), n_log_points)).astype(int)) - 1
    rank_idx = rank_idx[(rank_idx >= 0) & (rank_idx < n)]
    if len(rank_idx) < 3:
        return (min_threshold, int(counts[0]), 0)

    log_x = np.log10(rank_idx + 1)
    log_y = np.log10(np.maximum(counts[rank_idx], 1))

    x_range = log_x[-1] - log_x[0]
    y_range = log_y[0] - log_y[-1]
    if x_range <= 0 or y_range <= 0:
        return (min_threshold, int(counts[0]), 0)

    x_norm = (log_x - log_x[0]) / x_range
    y_norm = (log_y - log_y[-1]) / y_range  # plateau ~1, tail ~0
    diff = y_norm - (1.0 - x_norm)

    elbow_seg = int(np.argmax(diff))
    elbow_rank = int(rank_idx[elbow_seg])
    shoulder_count = int(counts[elbow_rank])
    threshold = max(min_threshold, shoulder_count // 2)
    return threshold, shoulder_count, elbow_rank


def _annotate_suggestion(ax, counts: np.ndarray, threshold: int,
                         shoulder_count: int, elbow_rank: int) -> None:
    """Overlay the suggested threshold + shoulder marker on the rank plot."""
    if counts.size == 0:
        return
    # Suggested threshold horizontal line (green)
    if threshold > 0:
        ax.axhline(threshold, color="tab:green", linewidth=1.4, linestyle="--", alpha=0.8)
    # Shoulder count level (lighter, dotted)
    if shoulder_count > 0:
        ax.axhline(shoulder_count, color="tab:orange", linewidth=1, linestyle=":", alpha=0.6)
    # Mark elbow point on the curve
    if 0 <= elbow_rank < counts.size:
        ax.plot([elbow_rank + 1], [max(counts[elbow_rank], 1)],
                marker="o", color="tab:orange", markersize=8,
                markeredgecolor="black", markeredgewidth=1, alpha=0.9, zorder=10)
    # Kept rank line for the suggested threshold
    below = np.where(counts < threshold)[0]
    kept_rank = int(below[0]) if below.size > 0 else int(counts.size)
    ax.axvline(kept_rank, color="tab:green", linewidth=1, linestyle=":", alpha=0.6)

    ax.text(
        0.02, 0.05,
        f"suggested thresh = {threshold}\nshoulder = {shoulder_count}\nkept_rank = {kept_rank:,}",
        transform=ax.transAxes, ha="left", va="bottom",
        fontsize=8, color="tab:green",
        bbox=dict(facecolor="white", edgecolor="tab:green", alpha=0.85, boxstyle="round,pad=0.3"),
    )


def plot_sample_rank(
    counts: np.ndarray,
    sample: Sample,
    output_path: Path,
    log_x: bool = False,
    xlim: int | None = None,
    annotate: bool = False,
    suggest: bool = False,
) -> Path:
    """Single-sample rank plot.

    Original `rankplot2.m` produces: semilogy of counts vs rank, title = sample
    number, xlabel "Sequence rank", ylabel "Read count". With shape=1 → loglog.
    With xlimit set → xlim([0 xlimit]).
    """
    fig, ax = plt.subplots(figsize=(7, 4.5))
    if counts.size:
        ranks = np.arange(1, counts.size + 1)
        if log_x:
            ax.loglog(ranks, counts, linewidth=1)
        else:
            ax.semilogy(ranks, counts, linewidth=1)
    if xlim:
        ax.set_xlim(0, xlim)
    ax.set_title(str(sample.sample_id))
    ax.set_xlabel("Sequence rank")
    ax.set_ylabel("Read count")
    if annotate:
        _annotate_threshold(ax, counts, sample.umi_threshold)
    if suggest:
        thr, sh, ei = suggest_umi_threshold(counts)
        _annotate_suggestion(ax, counts, thr, sh, ei)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_overview_grid(
    samples_counts: list[tuple[Sample, np.ndarray]],
    output_path: Path,
    log_x: bool = False,
    xlim: int | None = None,
    ncols: int = 4,
    annotate: bool = False,
    suggest: bool = False,
) -> Path:
    """Optional grid summary (off by default; not in original rankplot2.m)."""
    n = len(samples_counts)
    if n == 0:
        return output_path
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(3.5 * ncols, 2.5 * nrows), squeeze=False
    )
    for idx, (sample, counts) in enumerate(samples_counts):
        ax = axes[idx // ncols][idx % ncols]
        if counts.size:
            ranks = np.arange(1, counts.size + 1)
            if log_x:
                ax.loglog(ranks, counts, linewidth=0.8)
            else:
                ax.semilogy(ranks, counts, linewidth=0.8)
        if xlim:
            ax.set_xlim(0, xlim)
        ax.set_title(str(sample.sample_id), fontsize=9)
        ax.tick_params(labelsize=7)
        if annotate:
            _annotate_threshold(ax, counts, sample.umi_threshold)
        if suggest:
            thr, sh, ei = suggest_umi_threshold(counts)
            _annotate_suggestion(ax, counts, thr, sh, ei)
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    fig.supxlabel("Sequence rank", fontsize=10)
    fig.supylabel("Read count", fontsize=10)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return output_path
