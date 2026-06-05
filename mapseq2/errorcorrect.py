"""bowtie2 all-vs-all + connected components for local + global barcode collapse.

This replaces the implied MATLAB functions `produceBCmat2_nospikes`
(per-sample CC + pool) and the global-collapse portion of `produceBCmat3`.

The Bash pipeline writes per-sample bowtie2 SAM files; we collapse here and
emit the same per-sample collapsed counts plus a global CC label map.
"""

from __future__ import annotations

import logging
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from .config import Bowtie2Config, Config

log = logging.getLogger(__name__)


def run_bowtie2_all_vs_all(
    sequences: list[str],
    work_dir: Path,
    bowtie2_cfg: Bowtie2Config,
    score_min: str,
    name_prefix: str,
) -> Path:
    """Build a bowtie2 index from `sequences` and align them against itself.

    Returns the path to the produced SAM file. FASTA / index live in `work_dir`.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = work_dir / f"{name_prefix}.fasta"
    with fasta_path.open("w") as fh:
        for i, seq in enumerate(sequences):
            fh.write(f">{i}\n{seq}\n")

    index_prefix = work_dir / f"{name_prefix}_idx"
    build_log = work_dir / f"{name_prefix}_build.log"
    with build_log.open("w") as logf:
        subprocess.run(
            ["bowtie2-build", "-q", str(fasta_path), str(index_prefix)],
            check=True,
            stdout=logf,
            stderr=subprocess.STDOUT,
        )

    sam_path = work_dir / f"{name_prefix}.sam"
    align_log = work_dir / f"{name_prefix}_align.log"
    args = [
        "bowtie2",
        *bowtie2_cfg.args(score_min),
        "-x", str(index_prefix),
        "-f", "-U", str(fasta_path),
        "--no-head",
        "-S", str(sam_path),
    ]
    with align_log.open("w") as logf:
        subprocess.run(args, check=True, stdout=logf, stderr=subprocess.STDOUT)

    return sam_path


def parse_sam_edges(sam_path: Path) -> list[tuple[int, int]]:
    """Extract (query_idx, ref_idx) edges from a bowtie2 SAM (no header)."""
    edges: list[tuple[int, int]] = []
    with sam_path.open() as fh:
        for line in fh:
            if not line or line.startswith("@"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            qname, _flag, rname = parts[0], parts[1], parts[2]
            if rname == "*":
                continue
            edges.append((int(qname), int(rname)))
    return edges


def find_connected_components(n: int, edges: Iterable[tuple[int, int]]) -> np.ndarray:
    """Return CC label per node. Isolated nodes get unique labels."""
    edges = list(edges)
    if not edges:
        return np.arange(n, dtype=np.int64)
    rows, cols = zip(*edges)
    data = np.ones(len(edges), dtype=np.uint8)
    adj = csr_matrix((data, (rows, cols)), shape=(n, n))
    _, labels = connected_components(adj, directed=False)
    return labels.astype(np.int64)


def collapse_by_components(
    sequences: list[str],
    counts: list[int],
    labels: np.ndarray,
) -> list[tuple[str, int]]:
    """Group sequences by CC label; representative = member with highest count,
    aggregated count = sum within CC. Output sorted by aggregated count desc."""
    agg_count: dict[int, int] = defaultdict(int)
    repr_idx: dict[int, int] = {}
    for i, lab in enumerate(labels):
        lab = int(lab)
        agg_count[lab] += counts[i]
        if lab not in repr_idx or counts[i] > counts[repr_idx[lab]]:
            repr_idx[lab] = i
    out: list[tuple[str, int]] = []
    for lab, total in agg_count.items():
        out.append((sequences[repr_idx[lab]], total))
    out.sort(key=lambda kv: -kv[1])
    return out


def read_quickout(path: Path) -> tuple[list[int], list[str]]:
    counts: list[int] = []
    seqs: list[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) < 2:
                continue
            counts.append(int(parts[0]))
            seqs.append(parts[1])
    return counts, seqs


def local_collapse_sample(
    sample_id: str,
    quickout_path: Path,
    out_dir: Path,
    config: Config,
) -> Path:
    """Run bowtie2 all-vs-all on a per-sample quickout, find CC,
    write collapsed (count, sequence) TSV. Returns the path."""
    counts, seqs = read_quickout(quickout_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    collapsed_path = out_dir / f"{sample_id}_local_collapsed.tsv"

    if not seqs:
        collapsed_path.write_text("")
        return collapsed_path

    if len(seqs) == 1:
        collapsed_path.write_text(f"{counts[0]}\t{seqs[0]}\n")
        return collapsed_path

    sample_dir = out_dir / f"_workdir_{sample_id}"
    sam_path = run_bowtie2_all_vs_all(
        seqs, sample_dir, config.bowtie2,
        config.bowtie2.score_min_local,
        name_prefix=f"local_{sample_id}",
    )
    edges = parse_sam_edges(sam_path)
    labels = find_connected_components(len(seqs), edges)
    collapsed = collapse_by_components(seqs, counts, labels)

    with collapsed_path.open("w") as fh:
        for seq, c in collapsed:
            fh.write(f"{c}\t{seq}\n")

    cc_path = out_dir / f"{sample_id}_local_cc_labels.tsv"
    with cc_path.open("w") as fh:
        fh.write("idx\tcc_label\tsequence\tcount\n")
        for i, lab in enumerate(labels):
            fh.write(f"{i}\t{int(lab)}\t{seqs[i]}\t{counts[i]}\n")

    return collapsed_path


def pool_for_global(
    sample_collapsed_paths: dict[str, Path],
) -> tuple[list[str], dict[str, dict[str, int]]]:
    """Pool locally-collapsed sequences across samples.

    Returns:
      pooled_seqs: deterministic-ordered unique sequences across samples
      per_sample_counts: {sample_id: {seq: count}}
    """
    per_sample_counts: dict[str, dict[str, int]] = {}
    all_seqs: set[str] = set()
    for sample_id, path in sample_collapsed_paths.items():
        sc: dict[str, int] = {}
        if path.exists():
            with path.open() as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 2:
                        continue
                    sc[parts[1]] = int(parts[0])
        per_sample_counts[sample_id] = sc
        all_seqs.update(sc.keys())
    pooled_seqs = sorted(all_seqs)
    return pooled_seqs, per_sample_counts


def identify_spike_ccs(
    pooled_seqs: list[str],
    global_labels: np.ndarray,
    spike_tag: str,
) -> np.ndarray:
    """CC-level spike-in flag (port of produceBCmat3 spike detection).

    A global CC is flagged as a spike CC if ANY of its member pooled sequences
    ends with `spike_tag`. This catches sequencing variants of the spike whose
    representative might not retain the exact tag.

    Returns a boolean array indexed by global CC label.
    """
    if len(global_labels) == 0:
        return np.zeros(0, dtype=bool)
    n_cc = int(global_labels.max()) + 1
    is_spike = np.zeros(n_cc, dtype=bool)
    for i, seq in enumerate(pooled_seqs):
        if seq.endswith(spike_tag):
            is_spike[int(global_labels[i])] = True
    return is_spike


def global_collapse(
    pooled_seqs: list[str],
    out_dir: Path,
    config: Config,
) -> tuple[np.ndarray, Path]:
    """Run bowtie2 all-vs-all on pooled sequences. Returns (cc_labels, sam_path)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    if len(pooled_seqs) == 0:
        return np.array([], dtype=np.int64), out_dir / "global.sam"
    if len(pooled_seqs) == 1:
        return np.array([0], dtype=np.int64), out_dir / "global.sam"

    work = out_dir / "_workdir_global"
    sam_path = run_bowtie2_all_vs_all(
        pooled_seqs, work, config.bowtie2,
        config.bowtie2.score_min_global,
        name_prefix="global",
    )
    edges = parse_sam_edges(sam_path)
    labels = find_connected_components(len(pooled_seqs), edges)
    cc_path = out_dir / "global_cc_labels.tsv"
    with cc_path.open("w") as fh:
        fh.write("pooled_idx\tcc_label\tsequence\n")
        for i, lab in enumerate(labels):
            fh.write(f"{i}\t{int(lab)}\t{pooled_seqs[i]}\n")
    return labels, sam_path
