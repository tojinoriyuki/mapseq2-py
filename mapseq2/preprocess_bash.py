"""Bash-backed preprocess: wraps Hyopil's bash pipeline (pigz/awk/cut/grep/sort)
for ~10x speedup over the pure-Python implementation in preprocess.py.

The pure-Python `preprocess.py` is preserved as the reference / portable fallback.
Both produce the same per-sample outputs (`{sample_id}_BU.tsv`, `_quickout.tsv`,
`_ssi_check.tsv`, `_stats.tsv`, `_ssi_check.png`).
"""

from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path

from .config import Config, Sample, load_sample_sheet
from .viz import plot_ssi_check

log = logging.getLogger(__name__)

BASH_SCRIPT = Path(__file__).parent / "bash_wrappers" / "preprocess.sh"


def run_preprocess_bash(cfg: Config) -> None:
    """Run the bash preprocess pipeline for all samples in the sheet.

    Generates SSI.txt + threshold array from sample_sheet.csv + config.yaml,
    invokes the bash script, then converts Hyopil-format outputs (`<count> <seq>`)
    to TSV (`count\\tseq`) and writes stats + SSI plots.
    """
    samples = load_sample_sheet(cfg.input.sample_sheet)
    out = cfg.project.output_dir / "01_preprocess"
    out.mkdir(parents=True, exist_ok=True)

    # 1) Generate SSI.txt: line N = ssi for sample_id N (1-indexed)
    by_id = {int(s.sample_id): s for s in samples}
    max_id = max(by_id.keys())
    ssi_file = out / "SSI.txt"
    with ssi_file.open("w") as f:
        for i in range(1, max_id + 1):
            s = by_id.get(i)
            f.write((s.ssi_full_sequence if s else "N" * cfg.barcode.ssi_match_length) + "\n")

    # 2) Threshold array: leading 0 placeholder, then per-sample umi_threshold
    threshold_vals = ["0"]
    for i in range(1, max_id + 1):
        s = by_id.get(i)
        threshold_vals.append(str(s.umi_threshold) if s else "1")
    threshold_str = " ".join(threshold_vals)

    # 3) Build env for bash script
    env = os.environ.copy()
    env.update({
        "PREFIX": cfg.project.name,
        "RAW_FASTQ_DIR": str(cfg.input.raw_fastq_dir),
        "SSI_FILE": str(ssi_file),
        "OUT_DIR": str(out),
        "YY_REGEX": (cfg.library_filter.regex or "") if cfg.library_filter.enabled else "",
        "THRESHOLD_ARRAY": threshold_str,
        "NSAMPLES": str(max_id),
        "BARCODE_LEN": str(cfg.barcode.barcode_length),
        "UMI_LEN": str(cfg.barcode.umi_length),
        "SSI_LEN": str(cfg.barcode.ssi_length),
        "SSI_MATCH_LEN": str(cfg.barcode.ssi_match_length),
        "FASTQ_R1_PATTERN": cfg.input.r1_pattern,
        "FASTQ_R2_PATTERN": cfg.input.r2_pattern,
        "SORT_PARALLEL": env.get("MAPSEQ2_SORT_PARALLEL", "8"),
        "SORT_BUF": env.get("MAPSEQ2_SORT_BUF", "2G"),
        "PIGZ_THREADS": env.get("MAPSEQ2_PIGZ_THREADS", "4"),
        "TMPDIR": env.get("MAPSEQ2_TMPDIR", str(Path.home() / "scratch_mapseq2_sort")),
    })
    Path(env["TMPDIR"]).mkdir(parents=True, exist_ok=True)

    log.info("running bash preprocess: %s", BASH_SCRIPT)
    log.info("  PREFIX=%s NSAMPLES=%s YY=%s", env["PREFIX"], env["NSAMPLES"], env["YY_REGEX"])
    log.info("  threshold=%s", env["THRESHOLD_ARRAY"])
    subprocess.run(["bash", str(BASH_SCRIPT)], env=env, check=True)

    # 4) Convert Hyopil-format outputs to TSV + stats + plots
    intermediate = out / "intermediate"
    thresholds_dir = out / "thresholds"

    for s in samples:
        i = int(s.sample_id)
        prefix = cfg.project.name

        bu_src = intermediate / f"{prefix}_BC{i}_BU.txt"
        bu_dst = out / f"{s.sample_id}_BU.tsv"
        _convert_count_to_tsv(bu_src, bu_dst)

        qo_src = thresholds_dir / f"{prefix}_BC{i}_quickout.txt"
        qo_dst = out / f"{s.sample_id}_quickout.tsv"
        _convert_count_to_tsv(qo_src, qo_dst)

        ssi_src = intermediate / f"{prefix}_BC{i}_SSIcheck.txt"
        ssi_dst = out / f"{s.sample_id}_ssi_check.tsv"
        _convert_count_to_tsv(ssi_src, ssi_dst)

        if bu_dst.exists() and ssi_dst.exists():
            _write_stats(s, bu_dst, ssi_dst, cfg, out)
            try:
                _plot_ssi_check(s, ssi_dst, cfg, out)
            except Exception as e:
                log.warning("ssi-check plot failed for %s: %s", s.sample_id, e)

    log.info("bash preprocess done; outputs in %s", out)


def _convert_count_to_tsv(src: Path, dst: Path) -> None:
    """Convert Hyopil-format `<leading_ws><count> <seq>` to TSV `count\\tseq`."""
    if not src.exists():
        return
    with src.open() as fin, dst.open("w") as fout:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            parts = line.split(None, 1)
            if len(parts) >= 2:
                fout.write(f"{parts[0]}\t{parts[1]}\n")


def _write_stats(sample: Sample, bu_path: Path, ssi_path: Path,
                 cfg: Config, out: Path) -> None:
    expected_ssi = sample.ssi_match(cfg.barcode.ssi_match_length)
    ssi_total = 0
    ssi_matched = 0
    if ssi_path.exists():
        with ssi_path.open() as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    c = int(parts[0])
                    ssi_total += c
                    if parts[1] == expected_ssi:
                        ssi_matched += c

    bu_unique = 0
    bu_total_reads = 0
    if bu_path.exists():
        with bu_path.open() as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    bu_unique += 1
                    bu_total_reads += int(parts[0])

    stats_path = out / f"{sample.sample_id}_stats.tsv"
    with stats_path.open("w") as f:
        f.write("metric\tvalue\n")
        f.write(f"reads_after_yy\t{ssi_total}\n")
        f.write(f"reads_after_ssi\t{ssi_matched}\n")
        f.write(f"unique_bu\t{bu_unique}\n")
        f.write(f"total_reads_kept\t{bu_total_reads}\n")


def _plot_ssi_check(sample: Sample, ssi_path: Path, cfg: Config, out: Path) -> None:
    counts: dict[str, int] = {}
    with ssi_path.open() as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                counts[parts[1]] = int(parts[0])
    if not counts:
        return
    plot_ssi_check(
        counts, sample.sample_id,
        sample.ssi_match(cfg.barcode.ssi_match_length),
        out / f"{sample.sample_id}_ssi_check.png",
    )
