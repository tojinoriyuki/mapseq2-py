"""Test rankplot with PCR-amplification-simulated FASTQ.

Generates BU 44-mers where some have many reads (true sequences with high PCR
amplification) and many have few reads (singletons / errors), so the rank plot
shows the characteristic elbow.
"""

import gzip
import random
import shutil
import sys
from pathlib import Path
from textwrap import dedent

PROJ_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ_ROOT))

from mapseq2.config import Config
from mapseq2.pipeline import run_preprocess, run_rankplot, run_umi_collapse


BC_LEN = 32
UMI_LEN = 12
SSI_LEN = 8


def rand_dna(rng, n):
    return "".join(rng.choices("ACGT", k=n))


def make_yy_compliant(rng):
    """Generate a 32-mer ending in [TC][TC] or [TC][TC]G."""
    body = rand_dna(rng, BC_LEN - 2)
    tail = rng.choice(["TT", "TC", "CT", "CC"])
    return body + tail


def main():
    rng = random.Random(42)
    work = PROJ_ROOT / "tests" / "_rankplot_work"
    if work.exists():
        shutil.rmtree(work)
    raw = work / "raw"
    raw.mkdir(parents=True)
    out = work / "out"

    n_real = 30
    real_bcs = [make_yy_compliant(rng) for _ in range(n_real)]

    ssi = "ACGTACGT"

    r1_records: list[tuple[str, str]] = []
    r2_records: list[tuple[str, str]] = []

    for idx, bc in enumerate(real_bcs):
        n_umis = 10
        for u in range(n_umis):
            umi = rand_dna(rng, UMI_LEN)
            pcr_count = rng.randint(20, 200)
            for _ in range(pcr_count):
                r1_records.append((bc, "I" * BC_LEN))
                r2_records.append((umi + ssi, "I" * (UMI_LEN + SSI_LEN)))

    for _ in range(300):
        bc = make_yy_compliant(rng)
        umi = rand_dna(rng, UMI_LEN)
        r1_records.append((bc, "I" * BC_LEN))
        r2_records.append((umi + ssi, "I" * (UMI_LEN + SSI_LEN)))

    rng.shuffle(r1_records)
    rng.shuffle(r2_records)
    print(f"generated {len(r1_records)} reads")

    def write_fq(path, recs):
        with gzip.open(path, "wt") as f:
            for i, (s, q) in enumerate(recs):
                f.write(f"@r{i}\n{s}\n+\n{q}\n")

    write_fq(raw / "TEST-1_S1_L001_R1_001.fastq.gz", r1_records)
    write_fq(raw / "TEST-1_S1_L001_R2_001.fastq.gz", r2_records)

    (work / "sample_sheet.csv").write_text(dedent(
        f"""\
        sample_id,region,side,role,ssi_full_sequence,umi_threshold,sort_order
        1,Inj,Ip,source,{ssi}NNNNNNNNNNNN,1,
        """
    ))

    (work / "config.yaml").write_text(dedent(
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
          barcode_length: {BC_LEN}
          umi_length: {UMI_LEN}
          ssi_length: {SSI_LEN}
          ssi_match_length: {SSI_LEN}
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
          threads: 2
        filtering:
          source_threshold_umi: 1
          proj_threshold_umi: 1
          cell_body_threshold: null
        quality_filter:
          homopolymer_min_run: 7
        normalization:
          method: rowsum_log_scale
          log_factor: 100
        clustering:
          method: ward
          metric: euclidean
          num_clusters: 2
          cluster_labels: {{1: A, 2: B}}
        """
    ))

    cfg = Config.from_yaml(work / "config.yaml")

    print("\n=== preprocess (initial threshold=1) ===")
    run_preprocess(cfg)

    print("\n=== rankplot (MATLAB-aligned default) ===")
    run_rankplot(cfg, log_x=True)
    rp_dir = out / "01_preprocess" / "_rankplots"
    assert (rp_dir / "sample_1_rank.png").exists()
    assert not (rp_dir / "overview.png").exists(), "overview should be off by default"
    print(f"per-sample plot: {(rp_dir / 'sample_1_rank.png').name}")

    print("\n=== rankplot --overview --annotate (opt-in extras) ===")
    run_rankplot(cfg, log_x=True, overview=True, annotate=True)
    assert (rp_dir / "overview.png").exists(), "overview should be created when overview=True"
    print(f"overview + annotated plots created")

    print("\n=== umi-collapse re-run with threshold=10 ===")
    (work / "sample_sheet.csv").write_text(dedent(
        f"""\
        sample_id,region,side,role,ssi_full_sequence,umi_threshold,sort_order
        1,Inj,Ip,source,{ssi}NNNNNNNNNNNN,10,
        """
    ))
    cfg2 = Config.from_yaml(work / "config.yaml")
    run_umi_collapse(cfg2)

    quickout = out / "01_preprocess" / "1_quickout.tsv"
    lines = [l for l in quickout.read_text().splitlines() if l.strip()]
    print(f"after re-collapse (threshold=10): {len(lines)} barcodes in quickout")
    assert len(lines) >= 1
    print("\nrankplot + umi-collapse test PASSED")


if __name__ == "__main__":
    main()
