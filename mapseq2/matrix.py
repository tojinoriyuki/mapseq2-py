"""Barcode-by-region matrix construction, spike-in split, normalization, filtering.

Replaces MATLAB `produceBCmat3` (matrix build + spike split) and `normBCmat2`
(spike-in normalization + threshold filtering).
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .config import Sample


def build_barcode_matrix(
    pooled_seqs: list[str],
    global_labels: np.ndarray,
    per_sample_counts: dict[str, dict[str, int]],
    samples: list[Sample],
) -> tuple[np.ndarray, list[str], list[str]]:
    """Aggregate per-sample counts into a [n_global_barcodes × n_samples] matrix.

    Returns:
      matrix: int64 array
      refbarcodes: representative sequence per global ID (highest total-count member)
      sample_ids: column labels (order matches `samples`)
    """
    seq_to_global: dict[str, int] = {
        seq: int(global_labels[i]) for i, seq in enumerate(pooled_seqs)
    }
    if len(global_labels) == 0:
        n_global = 0
    else:
        n_global = int(global_labels.max()) + 1
    sample_ids = [s.sample_id for s in samples]
    n_samples = len(samples)

    matrix = np.zeros((n_global, n_samples), dtype=np.int64)
    member_totals: dict[tuple[int, str], int] = defaultdict(int)

    for col_idx, s in enumerate(samples):
        sample_counts = per_sample_counts.get(s.sample_id, {})
        for seq, count in sample_counts.items():
            global_id = seq_to_global.get(seq)
            if global_id is None:
                continue
            matrix[global_id, col_idx] += count
            member_totals[(global_id, seq)] += count

    repr_best: dict[int, tuple[str, int]] = {}
    for (gid, seq), total in member_totals.items():
        cur = repr_best.get(gid)
        if cur is None or total > cur[1]:
            repr_best[gid] = (seq, total)
    refbarcodes = [repr_best.get(i, ("", 0))[0] for i in range(n_global)]

    return matrix, refbarcodes, sample_ids


def split_spike_ins(
    matrix: np.ndarray,
    refbarcodes: list[str],
    spike_tag: str,
) -> tuple[np.ndarray, list[str], np.ndarray, list[str]]:
    """Partition rows by representative-end check (simple/legacy variant).

    For CC-level spike detection that matches Hyopil's `produceBCmat3`, use
    `split_spike_ins_by_cc` together with `errorcorrect.identify_spike_ccs`.
    """
    is_spike = np.array([bc.endswith(spike_tag) for bc in refbarcodes], dtype=bool)
    return _split_by_mask(matrix, refbarcodes, is_spike)


def split_spike_ins_by_cc(
    matrix: np.ndarray,
    refbarcodes: list[str],
    is_spike_cc: np.ndarray,
) -> tuple[np.ndarray, list[str], np.ndarray, list[str]]:
    """Partition rows by per-CC spike flag (port of produceBCmat3 logic).

    `is_spike_cc[k]` must be True iff the global CC with label k contains at
    least one sequence whose tail equals the spike tag.
    """
    if is_spike_cc.shape[0] != matrix.shape[0]:
        raise ValueError(
            f"is_spike_cc length {is_spike_cc.shape[0]} != matrix rows {matrix.shape[0]}"
        )
    return _split_by_mask(matrix, refbarcodes, is_spike_cc.astype(bool))


def _split_by_mask(matrix, refbarcodes, is_spike):
    real_mat = matrix[~is_spike]
    real_bcs = [bc for bc, sp in zip(refbarcodes, is_spike) if not sp]
    spike_mat = matrix[is_spike]
    spike_bcs = [bc for bc, sp in zip(refbarcodes, is_spike) if sp]
    return real_mat, real_bcs, spike_mat, spike_bcs


def has_homopolymer_run(seq: str, min_run: int) -> bool:
    """True if seq contains any base repeated `min_run` or more times in a row.

    Port of `findhomopolymers.m` (Hyopil): MATLAB uses minrunlength=7 by default.
    """
    if min_run <= 1 or not seq:
        return False
    run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            run += 1
            if run >= min_run:
                return True
        else:
            run = 1
    return False


def filter_homopolymers(
    matrix: np.ndarray,
    refbarcodes: list[str],
    min_run: int,
    extra: tuple[np.ndarray, ...] = (),
) -> tuple[np.ndarray, list[str], np.ndarray, tuple[np.ndarray, ...]]:
    """Drop rows whose representative contains a homopolymer run of >= min_run.

    `extra` lets you keep auxiliary arrays (e.g. is_spike_cc) row-aligned to
    the surviving rows.
    Returns (matrix, refbarcodes, keep_mask, extra_filtered).
    """
    if min_run <= 1:
        keep = np.ones(len(refbarcodes), dtype=bool)
    else:
        keep = np.array(
            [not has_homopolymer_run(bc, min_run) for bc in refbarcodes], dtype=bool
        )
    filtered_bcs = [bc for bc, k in zip(refbarcodes, keep) if k]
    filtered_extra = tuple(arr[keep] for arr in extra)
    return matrix[keep], filtered_bcs, keep, filtered_extra


def normalize_by_spikes(
    real_mat: np.ndarray,
    spike_mat: np.ndarray,
    method: str = "unique_count",
) -> np.ndarray:
    """Per-column normalize by spike-in load per sample.

    method:
      "unique_count" (default, matches Hyopil normBCmat2.m): divide by the
          number of unique spike CCs detected in each sample, i.e.
          (spike_mat > 0).sum(axis=0). Equivalent to MATLAB
          `size(spikes(i).counts2u, 1)`.
      "umi_sum" (legacy / original Kebschull MAPseq): divide by total spike
          UMI counts in each sample, i.e. spike_mat.sum(axis=0).
    """
    if method == "unique_count":
        x = (spike_mat > 0).sum(axis=0).astype(np.float64)
    elif method == "umi_sum":
        x = spike_mat.sum(axis=0).astype(np.float64)
    else:
        raise ValueError(
            f"normalize_by_spikes: unknown method {method!r} "
            f"(expected 'unique_count' or 'umi_sum')"
        )
    x = np.where(x == 0, 1.0, x)
    return real_mat.astype(np.float64) / x[None, :]


def role_indices(samples: list[Sample], role: str) -> list[int]:
    return [i for i, s in enumerate(samples) if s.role == role]


def target_sort_indices(samples: list[Sample]) -> list[int]:
    """Indices of target samples in display sort order (sort_order ascending)."""
    targets = [(i, s) for i, s in enumerate(samples) if s.role == "target"]
    targets.sort(key=lambda t: (
        t[1].sort_order if t[1].sort_order is not None else 1_000_000,
        t[1].sample_id,
    ))
    return [i for i, _ in targets]


def filter_by_thresholds(
    real_mat: np.ndarray,
    real_bcs: list[str],
    source_indices: list[int],
    target_indices: list[int],
    source_threshold: float,
    proj_threshold: float,
    cell_body_threshold: Optional[float] = None,
) -> tuple[np.ndarray, list[str], np.ndarray]:
    """normBCmat2-style filtering on raw UMI matrix.

    Keep rows where max-over-source > source_threshold AND max-over-target > proj_threshold.
    Optionally drop rows where any target column exceeds cell_body_threshold.
    """
    n = real_mat.shape[0]
    src_max = (
        real_mat[:, source_indices].max(axis=1) if source_indices else np.zeros(n)
    )
    tgt_max = (
        real_mat[:, target_indices].max(axis=1) if target_indices else np.zeros(n)
    )
    keep = (src_max > source_threshold) & (tgt_max > proj_threshold)

    if cell_body_threshold is not None:
        sub_max = real_mat[:, target_indices].max(axis=1) if target_indices else np.zeros(n)
        keep &= sub_max < cell_body_threshold

    B = real_mat[keep]
    Bseq = [bc for i, bc in enumerate(real_bcs) if keep[i]]
    return B, Bseq, keep


def save_matrix(
    matrix: np.ndarray,
    refbarcodes: list[str],
    sample_labels: list[str],
    path: Path,
) -> None:
    df = pd.DataFrame(matrix, index=refbarcodes, columns=sample_labels)
    df.index.name = "barcode"
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".parquet":
        df.reset_index().to_parquet(path, index=False)
    else:
        df.to_csv(path, sep="\t")


def load_matrix(path: Path) -> tuple[np.ndarray, list[str], list[str]]:
    if path.suffix == ".parquet":
        df = pd.read_parquet(path)
        df = df.set_index("barcode")
    else:
        df = pd.read_csv(path, sep="\t", index_col="barcode")
    return df.to_numpy(), df.index.tolist(), df.columns.tolist()
