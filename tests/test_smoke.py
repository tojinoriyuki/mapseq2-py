"""Smoke test: config load + preprocess on synthetic FASTQ fixture.

Constructs a minimal synthetic experiment with a known number of reads matching
SSI #1 and verifies the preprocess pipeline returns expected BU counts.
"""

import gzip
import sys
from pathlib import Path
from textwrap import dedent

PROJ_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ_ROOT))

from mapseq2.config import Config, load_sample_sheet  # noqa: E402
from mapseq2.preprocess import preprocess_sample, umi_collapse  # noqa: E402


def make_fastq(path: Path, records: list[tuple[str, str]]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as f:
        for i, (seq, qual) in enumerate(records):
            f.write(f"@read{i}\n{seq}\n+\n{qual}\n")


def main():
    work = PROJ_ROOT / "tests" / "_smoke_work"
    raw = work / "raw"
    if work.exists():
        import shutil

        shutil.rmtree(work)
    work.mkdir(parents=True)
    raw.mkdir(parents=True)

    bc_len = 32
    umi_len = 12
    ssi_match_len = 8

    bc_ok = "A" * 30 + "TT"
    bc_ok_alt = "A" * 29 + "TCG"
    bc_bad = "A" * 30 + "AA"

    ssi1 = "ACGTACGT"
    ssi2 = "TGCATGCA"

    def mkread(bc: str, umi: str, ssi: str) -> str:
        assert len(bc) == bc_len
        assert len(umi) == umi_len
        assert len(ssi) == ssi_match_len
        return bc + umi + ssi

    full_records: list[tuple[str, str]] = []

    def add(seq, n=1):
        for _ in range(n):
            full_records.append((seq, "I" * len(seq)))

    add(mkread(bc_ok, "C" * umi_len, ssi1), n=5)
    add(mkread(bc_ok, "T" * umi_len, ssi1), n=3)
    add(mkread(bc_ok_alt, "G" * umi_len, ssi1), n=4)
    add(mkread(bc_bad, "C" * umi_len, ssi1), n=2)
    add(mkread(bc_ok, "A" * umi_len, ssi2), n=6)
    add(mkread(bc_ok, "C" * umi_len, ssi1).replace("A", "N", 1), n=1)

    r1_records = [(r[:bc_len], q[:bc_len]) for r, q in full_records]
    r2_records = [(r[bc_len:], q[bc_len:]) for r, q in full_records]

    make_fastq(raw / "TEST-1_S1_L001_R1_001.fastq.gz", r1_records)
    make_fastq(raw / "TEST-1_S1_L001_R2_001.fastq.gz", r2_records)
    make_fastq(raw / "TEST-2_S2_L001_R1_001.fastq.gz", [])
    make_fastq(raw / "TEST-2_S2_L001_R2_001.fastq.gz", [])

    sample_sheet = work / "sample_sheet.csv"
    sample_sheet.write_text(
        dedent(
            f"""\
            sample_id,region,side,role,ssi_full_sequence,umi_threshold,sort_order
            1,RegionA,Ip,target,{ssi1}NNNNNNNNNNNN,3,1
            2,RegionB,Con,source,{ssi2}NNNNNNNNNNNN,3,
            """
        )
    )

    cfg_path = work / "config.yaml"
    cfg_path.write_text(
        dedent(
            f"""\
            project:
              name: TEST
              output_dir: ./out
            input:
              raw_fastq_dir: ./raw
              r1_pattern: "{{prefix}}-{{sample_id}}_*_L*_R1_001.fastq.gz"
              r2_pattern: "{{prefix}}-{{sample_id}}_*_L*_R2_001.fastq.gz"
              sample_sheet: ./sample_sheet.csv
            barcode:
              barcode_length: {bc_len}
              umi_length: {umi_len}
              ssi_length: {ssi_match_len}
              ssi_match_length: {ssi_match_len}
            library_filter:
              enabled: true
              regex: '^(.{{30}}[TC][TC]|.{{29}}[TC][TC]G)'
            spike_in:
              tag_sequence: CGTCAGTC
            bowtie2:
              N: 1
              L: 16
              D: 20
              R: 3
              i: "S,1,1.15"
              score_min_local: "L,0,-0.7"
              score_min_global: "L,0,-0.6"
              threads: 4
            filtering:
              source_threshold_umi: 20
              proj_threshold_umi: 10
              cell_body_threshold: null
            normalization:
              method: rowsum_log_scale
              log_factor: 100
            clustering:
              method: ward
              metric: euclidean
              num_clusters: 2
              cluster_labels: {{1: A, 2: B}}
            """
        )
    )

    cfg = Config.from_yaml(cfg_path)
    samples = load_sample_sheet(cfg.input.sample_sheet)
    print(f"config loaded: project={cfg.project.name}, {len(samples)} samples")

    sample1 = next(s for s in samples if s.sample_id == "1")
    out = preprocess_sample(sample1, cfg)
    print(
        f"sample {sample1.sample_id}: "
        f"reads_seen={out.reads_seen}, after_yy={out.reads_after_yy}, "
        f"after_ssi={out.reads_after_ssi}, bu_unique={len(out.bu_counts)}"
    )
    print(f"SSI check top entries: {out.ssi_check.most_common(3)}")
    print(f"BU counts: {dict(out.bu_counts)}")

    assert out.reads_seen == 21, f"expected 21 total reads, got {out.reads_seen}"
    assert out.reads_after_yy == 18, f"expected 18 after YY+N filter, got {out.reads_after_yy}"
    assert out.reads_after_ssi == 12, f"expected 12 after SSI match, got {out.reads_after_ssi}"
    assert len(out.bu_counts) == 3, f"expected 3 unique BU mers, got {len(out.bu_counts)}"

    quick = umi_collapse(out.bu_counts, cfg.barcode.barcode_length, threshold=3)
    print(f"UMI collapse (thresh=3): {dict(quick)}")
    assert quick[bc_ok] == 2
    assert quick.get(bc_ok_alt, 0) == 1

    sample2 = next(s for s in samples if s.sample_id == "2")
    out2 = preprocess_sample(sample2, cfg)
    print(
        f"sample {sample2.sample_id} (empty fastq): "
        f"reads_seen={out2.reads_seen}, bu_unique={len(out2.bu_counts)}"
    )
    assert out2.reads_seen == 0

    print("\nALL SMOKE TESTS PASSED")


if __name__ == "__main__":
    main()
