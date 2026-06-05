"""Per-row log normalization + Ward hierarchical clustering.

Replaces the clustering portion of MAPseq1andMAPseq2clustering.mlx.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage


@dataclass
class ClusterResult:
    data_scaled: np.ndarray
    linkage_Z: np.ndarray
    perm: np.ndarray
    labels: np.ndarray
    sorted_data: np.ndarray
    sorted_labels: np.ndarray


def log_norm_scale(matrix: np.ndarray, log_factor: float = 100) -> np.ndarray:
    """Per-row L1 normalize → log(1 + x*factor) → scale by global max to [0,1]."""
    matrix = matrix.astype(np.float64)
    row_sums = matrix.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0, 1.0, row_sums)
    data_norm = matrix / row_sums
    log_data = np.log(1.0 + data_norm * log_factor)
    max_val = log_data.max()
    if max_val == 0:
        return log_data
    return log_data / max_val


def cluster_projection_matrix(
    matrix: np.ndarray,
    num_clusters: int,
    method: str = "ward",
    metric: str = "euclidean",
    log_factor: float = 100,
    rng_seed: int = 0,
) -> ClusterResult:
    """Cluster a projection matrix (rows = neurons, cols = target regions).

    Matches the MATLAB pipeline: randomize row order, log-norm, Ward, fcluster.
    Returns sorted data and labels in dendrogram leaf order.
    """
    rng = np.random.default_rng(rng_seed)
    perm_input = rng.permutation(matrix.shape[0])
    shuffled = matrix[perm_input]

    scaled = log_norm_scale(shuffled, log_factor=log_factor)

    Z = linkage(scaled, method=method, metric=metric)
    labels = fcluster(Z, t=num_clusters, criterion="maxclust").astype(np.int64)

    dd = dendrogram(Z, no_plot=True)
    perm = np.asarray(dd["leaves"], dtype=np.int64)

    return ClusterResult(
        data_scaled=scaled,
        linkage_Z=Z,
        perm=perm,
        labels=labels,
        sorted_data=scaled[perm],
        sorted_labels=labels[perm],
    )


def cluster_size_summary(labels: np.ndarray) -> dict[int, int]:
    unique, counts = np.unique(labels, return_counts=True)
    return {int(u): int(c) for u, c in zip(unique, counts)}
