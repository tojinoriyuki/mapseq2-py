"""FASTQ extraction + per-sample SSI demultiplex + UMI collapse.

This module replaces the first half of bash_processing_v2.txt:
  - Step 1: extract barcode (R1[:32]) + UMI+SSI (R2[:20]), apply library YY filter
  - Step 2: demultiplex by per-sample SSI (last `ssi_match_length` bases)
  - Step 3: quick UMI collapse using per-sample threshold

All hardcoded parameters from the bash (cut lengths, YY regex, spike tag, thresholds)
come from the Config + sample sheet.
"""

from __future__ import annotations

import gzip
import logging
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Optional

from .config import BarcodeConfig, Config, LibraryFilter, Sample

log = logging.getLogger(__name__)


def _open_fastq(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")


def _iter_seq_lines(path: Path) -> Iterator[str]:
    """Yield sequence lines (every 4th line starting from line 2) from a FASTQ file."""
    with _open_fastq(path) as f:
        for idx, line in enumerate(f):
            if idx % 4 == 1:
                yield line.rstrip("\n")


def iter_paired_combined(
    r1_paths: list[Path],
    r2_paths: list[Path],
    barcode_len: int,
    r2_extract_len: int,
) -> Iterator[str]:
    """Yield concatenated reads: R1[:barcode_len] + R2[:r2_extract_len] per read pair."""
    if len(r1_paths) != len(r2_paths):
        raise ValueError(f"R1/R2 file count mismatch: {len(r1_paths)} vs {len(r2_paths)}")
    for r1_path, r2_path in zip(r1_paths, r2_paths):
        if not r1_path.exists():
            raise FileNotFoundError(r1_path)
        if not r2_path.exists():
            raise FileNotFoundError(r2_path)
        r1_iter = _iter_seq_lines(r1_path)
        r2_iter = _iter_seq_lines(r2_path)
        for r1_seq, r2_seq in zip(r1_iter, r2_iter):
            yield r1_seq[:barcode_len] + r2_seq[:r2_extract_len]


def find_fastq_files(
    raw_dir: Path, pattern: str, prefix: str, sample_id: str
) -> list[Path]:
    """Resolve glob pattern with prefix/sample_id substitution and return sorted file list."""
    formatted = pattern.format(prefix=prefix, sample_id=sample_id)
    paths = sorted(raw_dir.glob(formatted))
    return paths


@dataclass
class PreprocessOutputs:
    bu_counts: Counter[str]
    ssi_check: Counter[str]
    reads_seen: int
    reads_after_yy: int
    reads_after_ssi: int


def preprocess_sample(
    sample: Sample,
    config: Config,
    progress_every: int = 0,
) -> PreprocessOutputs:
    """Process one sample's FASTQs end-to-end through library filter and SSI demux.

    Returns BU-mer (barcode+UMI) counts and SSI distribution for QC.
    """
    r1_paths = find_fastq_files(
        config.input.raw_fastq_dir,
        config.input.r1_pattern,
        config.project.name,
        sample.sample_id,
    )
    r2_paths = find_fastq_files(
        config.input.raw_fastq_dir,
        config.input.r2_pattern,
        config.project.name,
        sample.sample_id,
    )
    if not r1_paths:
        raise FileNotFoundError(
            f"no R1 fastqs found for sample {sample.sample_id}: pattern={config.input.r1_pattern}"
        )

    bc = config.barcode
    yy = config.library_filter.compile()
    ssi_match = sample.ssi_match(bc.ssi_match_length)
    ssi_slice = slice(bc.bu_length, bc.bu_length + bc.ssi_match_length)

    bu_counts: Counter[str] = Counter()
    ssi_check: Counter[str] = Counter()
    reads_seen = 0
    reads_after_yy = 0
    reads_after_ssi = 0

    for combined in iter_paired_combined(
        r1_paths, r2_paths, bc.barcode_length, bc.r2_extract_length
    ):
        reads_seen += 1
        if progress_every and reads_seen % progress_every == 0:
            log.info("sample %s: %d reads processed", sample.sample_id, reads_seen)

        if "N" in combined:
            continue
        if yy is not None and not yy.match(combined):
            continue
        reads_after_yy += 1

        observed_ssi = combined[ssi_slice]
        ssi_check[observed_ssi] += 1
        if observed_ssi != ssi_match:
            continue
        reads_after_ssi += 1

        bu_counts[combined[: bc.bu_length]] += 1

    return PreprocessOutputs(
        bu_counts=bu_counts,
        ssi_check=ssi_check,
        reads_seen=reads_seen,
        reads_after_yy=reads_after_yy,
        reads_after_ssi=reads_after_ssi,
    )


def umi_collapse(
    bu_counts: Counter[str],
    barcode_length: int,
    threshold: int,
) -> Counter[str]:
    """Collapse (BC+UMI) counts to per-barcode unique-UMI counts above read threshold.

    Mirrors the bash quick-collapse step. For each BU-mer sorted by read count desc,
    keep entries with read_count >= threshold (if threshold != 1, otherwise keep all),
    then count distinct UMIs per 32-mer barcode.

    The bash logic: find first rank where count < threshold, take top-N entries,
    cut to barcode, sort | uniq -c.
    """
    if threshold == 1 or threshold <= 0:
        kept_bus = list(bu_counts.keys())
    else:
        sorted_bus = sorted(bu_counts.items(), key=lambda kv: -kv[1])
        kept_bus = [bu for bu, c in sorted_bus if c >= threshold]

    barcode_umi_counts: Counter[str] = Counter()
    for bu in kept_bus:
        barcode_umi_counts[bu[:barcode_length]] += 1
    return barcode_umi_counts


def write_quickout(
    quickout: Counter[str],
    out_path: Path,
) -> None:
    """Write quickout in `<count> <sequence>` lines sorted by count desc.

    Matches bash output of `${Prefix}_BC${i}_quickout.txt`.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        for seq, count in sorted(quickout.items(), key=lambda kv: -kv[1]):
            f.write(f"{count}\t{seq}\n")
