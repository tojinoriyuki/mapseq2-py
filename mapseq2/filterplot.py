"""Source/projection-UMI-threshold scan plot (port of numproj_injthresh.m).

Default behavior matches the MATLAB original: a single stacked-bar plot of
fraction-of-kept-barcodes by `n target regions projected to`, parula/viridis
colormap, ylim [0,1]. Returns (percent, total) arrays — emitted as TSV.

Opt-in extras (off by default):
  - show_totals=True : add a lower panel of n_kept vs threshold (semilog)

The MATLAB function only scans source threshold. We accept a `scan` parameter
so the projection-threshold variant (analogous to numproj_projthresh.m which
also exists in the helperfunctions directory) is available with the same code.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def scan_thresholds(
    real_mat: np.ndarray,
    source_indices: list[int],
    target_indices: list[int],
    max_thresh: int,
    max_targets: int,
    fixed_proj_threshold: int = 0,
    fixed_source_threshold: int = 0,
    scan: Literal["source", "proj"] = "source",
) -> pd.DataFrame:
    """Numerical scan returning a tidy DataFrame.

    Columns: threshold, n_kept, k_1, k_2, ..., k_<max_targets>
    where k_i = #kept barcodes detected in exactly i target columns
    (UMI > 0); >max_targets binned into k_<max_targets>.
    """
    if real_mat.shape[0] == 0:
        return pd.DataFrame(
            columns=["threshold", "n_kept", *(f"k_{k}" for k in range(1, max_targets + 1))]
        )

    src_max = (
        real_mat[:, source_indices].max(axis=1)
        if source_indices else np.zeros(real_mat.shape[0])
    )
    tgt_max = (
        real_mat[:, target_indices].max(axis=1)
        if target_indices else np.zeros(real_mat.shape[0])
    )
    tgt_presence = (
        (real_mat[:, target_indices] > 0).sum(axis=1)
        if target_indices else np.zeros(real_mat.shape[0], dtype=np.int64)
    )
    tgt_presence_capped = np.minimum(tgt_presence, max_targets).astype(np.int64)

    rows = []
    for thresh in range(1, max_thresh + 1):
        if scan == "source":
            keep = (src_max > thresh) & (tgt_max > fixed_proj_threshold)
        else:
            keep = (tgt_max > thresh) & (src_max > fixed_source_threshold)
        total = int(keep.sum())
        if total == 0:
            counts = np.zeros(max_targets + 1, dtype=np.int64)
        else:
            counts = np.bincount(tgt_presence_capped[keep], minlength=max_targets + 1)
        row = {"threshold": thresh, "n_kept": total}
        for k in range(1, max_targets + 1):
            row[f"k_{k}"] = int(counts[k])
        rows.append(row)
    return pd.DataFrame(rows)


def plot_threshold_scan(
    scan_df: pd.DataFrame,
    output_path: Path,
    max_targets: int,
    scan_label: str = "source",
    show_totals: bool = False,
    title: str | None = None,
) -> Path:
    """Stacked bar of fractions, ylim [0,1], parula/viridis colormap.

    Matches numproj_injthresh.m's `bar(percent,'stacked',...)` plot.
    With show_totals=True, adds a lower semilog panel of n_kept (extension).
    """
    k_cols = [f"k_{k}" for k in range(1, max_targets + 1)]
    counts = scan_df[k_cols].to_numpy().astype(float)
    totals = counts.sum(axis=1)
    safe_totals = np.where(totals > 0, totals, 1.0)
    fractions = counts / safe_totals[:, None]

    cmap = plt.get_cmap("viridis", max_targets)
    colors = [cmap(i) for i in range(max_targets)]

    x = scan_df["threshold"].to_numpy()

    if show_totals:
        fig, (ax_bar, ax_total) = plt.subplots(
            2, 1, figsize=(11, 7), sharex=True,
            gridspec_kw=dict(height_ratios=[3, 1], hspace=0.08),
        )
    else:
        fig, ax_bar = plt.subplots(figsize=(10, 5))
        ax_total = None

    bottom = np.zeros(len(scan_df))
    for k in range(max_targets):
        ax_bar.bar(
            x, fractions[:, k], bottom=bottom, color=colors[k],
            width=1.0, edgecolor="none",
        )
        bottom += fractions[:, k]
    ax_bar.set_ylim(0, 1)
    if title:
        ax_bar.set_title(title)

    if ax_total is not None:
        ax_total.semilogy(x, scan_df["n_kept"], color="black", linewidth=1.2)
        ax_total.set_ylabel("# kept barcodes")
        ax_total.set_xlabel(f"{scan_label} UMI threshold")
        ax_total.set_xlim(x.min() - 0.5, x.max() + 0.5)
        ax_total.grid(True, which="both", alpha=0.3)
    else:
        ax_bar.set_xlim(x.min() - 0.5, x.max() + 0.5)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return output_path
