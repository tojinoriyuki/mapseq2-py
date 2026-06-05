"""Projection-matrix visualization (dendrogram + heatmap + cluster colorbar)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram


def _default_cluster_palette(num_clusters: int) -> np.ndarray:
    cmap = plt.get_cmap("tab10" if num_clusters <= 10 else "tab20")
    return np.array([cmap(i % cmap.N)[:3] for i in range(num_clusters)])


def plot_projection_heatmap(
    sorted_data: np.ndarray,
    sorted_labels: np.ndarray,
    Z: np.ndarray,
    perm: np.ndarray,
    region_labels: list[str],
    output_path: Path,
    cluster_id_to_name: dict[int, str] | None = None,
    cmap: str = "viridis",
) -> Path:
    unique_clusters = sorted(set(int(l) for l in sorted_labels))
    palette = _default_cluster_palette(len(unique_clusters))
    cluster_to_color = {cid: palette[i] for i, cid in enumerate(unique_clusters)}

    cluster_bar = np.array(
        [cluster_to_color[int(l)] for l in sorted_labels]
    ).reshape(-1, 1, 3)

    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(1, 100, wspace=0.05)

    ax_dendro = fig.add_subplot(gs[0, :10])
    dendrogram(
        Z,
        ax=ax_dendro,
        orientation="left",
        no_labels=True,
        link_color_func=lambda _k: "black",
        above_threshold_color="black",
    )
    ax_dendro.invert_yaxis()
    ax_dendro.axis("off")

    ax_cbar = fig.add_subplot(gs[0, 11:13])
    ax_cbar.imshow(cluster_bar, aspect="auto", interpolation="nearest")
    ax_cbar.set_xticks([])
    ax_cbar.set_yticks([])

    ax_heat = fig.add_subplot(gs[0, 14:])
    im = ax_heat.imshow(sorted_data, aspect="auto", cmap=cmap, interpolation="nearest")
    ax_heat.set_xticks(range(len(region_labels)))
    ax_heat.set_xticklabels(region_labels, rotation=90, fontsize=9)
    ax_heat.set_yticks([])
    ax_heat.yaxis.tick_right()

    cbar = fig.colorbar(im, ax=ax_heat, fraction=0.025, pad=0.02)
    cbar.set_label("scaled log proj")

    legend_handles = []
    for cid in unique_clusters:
        name = cluster_id_to_name.get(cid, f"cluster {cid}") if cluster_id_to_name else f"cluster {cid}"
        legend_handles.append(
            plt.Rectangle((0, 0), 1, 1, fc=cluster_to_color[cid], label=name)
        )
    if legend_handles:
        ax_heat.legend(
            handles=legend_handles,
            loc="upper left",
            bbox_to_anchor=(1.05, 1.0),
            fontsize=9,
            frameon=False,
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_ssi_check(
    ssi_counts: dict[str, int],
    sample_id: str,
    expected_ssi: str,
    output_path: Path,
    top_n: int = 10,
) -> Path:
    items = sorted(ssi_counts.items(), key=lambda kv: -kv[1])[:top_n]
    labels = [s for s, _ in items]
    values = [v for _, v in items]
    colors = ["tab:green" if s == expected_ssi else "tab:gray" for s in labels]

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(range(len(labels)), values, color=colors)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("read count")
    ax.set_title(f"sample {sample_id} SSI distribution (expected: {expected_ssi})")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return output_path
