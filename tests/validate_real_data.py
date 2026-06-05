"""Validation against Hyopil's published .mat files.

Strategy (option A in our discussion):
  1. Load barcodematrixMAPseq2-{1,2}.mat (post-CC, post-spike-split barcode matrix)
  2. Load spikesMAPseq2-{1,2}.mat (per-sample spike CC counts)
  3. Apply my filter + normalize stages directly
  4. Compare to Hyopil's MAPseq2_1 / MAPseq2_2 (the published filtered+normalized matrices)
  5. Concatenate, log-norm, cluster, compare cluster labels to Hyopil's T_mapseq2

Hyopil's normBCmat2 call:
    normBCmat2(barcodematrix, refbarcodes, spikes, 20, 10, 3, [1 2 4:27], sorting2)
  sourcethresh = 20
  projthresh   = 10
  sourcesite   = 3 (1-indexed MATLAB → index 2 in 0-indexed Python)
  projsite     = [1 2 4..27] (1-indexed) → all but the source
  sorting2 reorders the 26-target output

Note: normBCmat2.m overwrites B_tar with `B(:, sorting)` after first setting
B_tar = B(:, projsite). Since sorting2 has values 1..26 but B has 27 cols,
this is ambiguous. The intended interpretation per targets26 in the .mlx is that
sorting2 indexes into the post-projsite list. We follow that.
"""

import sys
from pathlib import Path

PROJ_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ_ROOT))

import h5py
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.io import loadmat

from mapseq2.matrix import filter_by_thresholds

REPO = Path("/home/noriyuki/MAPseq/MAPseq2-and-POINTseq/MAPseq2/MATLAB")

# 1-indexed → 0-indexed: source at MATLAB col 3 → Python col 2
SOURCE_COL_PY = 2
TARGET_COLS_PY = [i for i in range(27) if i != SOURCE_COL_PY]
SOURCE_THRESH = 20
PROJ_THRESH = 10

# sorting2 (1-indexed, MATLAB) indexes into the 26-target ordered list (post-projsite)
SORTING2_MATLAB_1IDX = [1, 10, 2, 3, 11, 4, 6, 8, 5, 7, 9, 14, 12, 18, 16,
                       15, 13, 19, 17, 20, 21, 22, 23, 25, 24, 26]
SORTING2_PY = [i - 1 for i in SORTING2_MATLAB_1IDX]


def load_barcodematrix(prefix):
    """Returns (matrix [n_bc, n_samples], refbarcodes [n_bc])."""
    with h5py.File(REPO / f"barcodematrix{prefix}.mat", "r") as h:
        bm = np.array(h["barcodematrix"])  # h5py: (27, n_bc)
        rb = np.array(h["refbarcodes"])    # h5py: (32, n_bc), int8 ASCII
    matrix = bm.T  # → (n_bc, 27)
    seq_chars = rb.T  # → (n_bc, 32)
    refbarcodes = ["".join(chr(b) for b in row) for row in seq_chars]
    return matrix, refbarcodes


def load_spike_totals(prefix):
    """Returns per-sample spike normalizer. Hyopil normBCmat2.m uses
    `size(spikes(i).counts2u, 1)` = count of UNIQUE spike CCs detected
    in sample i (not total UMI)."""
    d = loadmat(REPO / f"spikes{prefix}.mat")
    spikes = d["spikes"]
    norm_x = np.zeros(spikes.shape[1], dtype=np.float64)
    for i in range(spikes.shape[1]):
        c2u = np.asarray(spikes[0, i]["counts2u"])
        norm_x[i] = float(c2u.shape[0]) if c2u.size > 0 else 1.0
    return norm_x


def filter_normalize_sort(matrix, spike_x):
    """Reproduce normBCmat2(.,20,10,3,[1 2 4..27], sorting2).

    `spike_x` is the per-sample normalizer (unique-spike count per Hyopil).
    """
    real_mat, _, keep = filter_by_thresholds(
        matrix,
        [""] * matrix.shape[0],
        source_indices=[SOURCE_COL_PY],
        target_indices=TARGET_COLS_PY,
        source_threshold=SOURCE_THRESH,
        proj_threshold=PROJ_THRESH,
        cell_body_threshold=None,
    )
    safe_x = np.where(spike_x == 0, 1.0, spike_x)
    Bnorm = real_mat.astype(np.float64) / safe_x[None, :]
    Bnorm_targets_raw_order = Bnorm[:, TARGET_COLS_PY]
    Bnorm_tar_sorted = Bnorm_targets_raw_order[:, SORTING2_PY]
    return real_mat, Bnorm, Bnorm_tar_sorted, keep


def compare_matrices(my, hyopil, label):
    my = np.asarray(my)
    hy = np.asarray(hyopil)
    print(f"\n[{label}]  shapes: mine={my.shape} hyopil={hy.shape}")
    if my.shape != hy.shape:
        print(f"  ✗ SHAPE MISMATCH")
        return False
    abs_diff = np.abs(my - hy)
    print(f"  max abs diff = {abs_diff.max():.6g}")
    print(f"  mean abs diff = {abs_diff.mean():.6g}")
    nonzero_my = (my > 0).sum()
    nonzero_hy = (hy > 0).sum()
    print(f"  nonzero entries: mine={nonzero_my} hyopil={nonzero_hy}")
    if abs_diff.max() < 1e-6:
        print(f"  ✓ EXACT MATCH")
        return True
    if abs_diff.max() < 1e-4:
        print(f"  ~ NEAR MATCH (tiny numerical diff)")
        return True
    print(f"  → row sums: mine_max={my.sum(axis=1).max():.4f}, hyopil_max={hy.sum(axis=1).max():.4f}")
    return False


def main():
    expected = loadmat(REPO / "MAPseq1andMAPseq2clustering.mat")
    print("== expected keys (shapes) ==")
    for k in ["MAPseq2_1", "MAPseq2_2", "MAPseq2",
              "T_mapseq2", "Tsorted_mapseq2", "Z_mapseq2",
              "perm_mapseq2", "data_mapseq2", "randorder_mapseq2",
              "num_clusters_mapseq2"]:
        if k in expected:
            print(f"  {k}: shape={expected[k].shape} dtype={expected[k].dtype}")

    print("\n== loading MAPseq2-1 ==")
    bm1, refs1 = load_barcodematrix("MAPseq2-1")
    spike_totals_1 = load_spike_totals("MAPseq2-1")
    print(f"barcodematrix: {bm1.shape}, refs: {len(refs1)}")
    print(f"spike_totals_1 (first 5): {spike_totals_1[:5]}")
    real1, Bnorm1, my_Bnorm_tar_1, keep1 = filter_normalize_sort(bm1, spike_totals_1)
    print(f"after filter: {real1.shape[0]} / {bm1.shape[0]} barcodes kept ({keep1.sum()})")

    print("\n== loading MAPseq2-2 ==")
    bm2, refs2 = load_barcodematrix("MAPseq2-2")
    spike_totals_2 = load_spike_totals("MAPseq2-2")
    print(f"barcodematrix: {bm2.shape}, refs: {len(refs2)}")
    real2, Bnorm2, my_Bnorm_tar_2, keep2 = filter_normalize_sort(bm2, spike_totals_2)
    print(f"after filter: {real2.shape[0]} / {bm2.shape[0]} barcodes kept ({keep2.sum()})")

    print("\n== compare to Hyopil's MAPseq2_1, MAPseq2_2 ==")
    compare_matrices(my_Bnorm_tar_1, expected["MAPseq2_1"], "MAPseq2_1 (target,sorted)")
    compare_matrices(my_Bnorm_tar_2, expected["MAPseq2_2"], "MAPseq2_2 (target,sorted)")

    print("\n== compare unsorted (just target subset, raw order) ==")
    my_unsorted_1 = Bnorm1[:, TARGET_COLS_PY]
    my_unsorted_2 = Bnorm2[:, TARGET_COLS_PY]

    print("\n== concat + cluster (replicating .mlx) ==")
    MAPseq2 = np.vstack([my_Bnorm_tar_1, my_Bnorm_tar_2])
    print(f"MAPseq2 (concat): {MAPseq2.shape}")
    print(f"Hyopil MAPseq2:   {expected['MAPseq2'].shape}")

    randorder = expected["randorder_mapseq2"].ravel().astype(int) - 1
    print(f"randorder length: {len(randorder)}")
    if MAPseq2.shape[0] == len(randorder):
        inputdata = MAPseq2[randorder]
        row_sums = inputdata.sum(axis=1, keepdims=True)
        row_sums = np.where(row_sums == 0, 1.0, row_sums)
        datanorm = inputdata / row_sums
        logdata = np.log(1 + datanorm * 100)
        scaled = logdata / logdata.max() if logdata.max() > 0 else logdata
        compare_matrices(scaled, expected["data_mapseq2"], "data_mapseq2 (scaled)")

        Z = linkage(scaled, method="ward", metric="euclidean")
        labels = fcluster(Z, t=int(expected["num_clusters_mapseq2"].item()),
                          criterion="maxclust")
        hyopil_labels = expected["T_mapseq2"].ravel().astype(int)
        print(f"\n[cluster labels]")
        print(f"  mine:   {pd.Series(labels).value_counts().sort_index().to_dict()}")
        print(f"  hyopil: {pd.Series(hyopil_labels).value_counts().sort_index().to_dict()}")
        try:
            from sklearn.metrics import adjusted_rand_score
            ari = adjusted_rand_score(hyopil_labels, labels)
            print(f"  adjusted Rand Index: {ari:.4f} (1.0 = identical clustering)")
        except ImportError:
            print("  (sklearn not available, skipping ARI)")
    else:
        print(f"  shape mismatch with randorder, skipping cluster comparison")


if __name__ == "__main__":
    main()
