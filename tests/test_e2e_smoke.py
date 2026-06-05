"""End-to-end smoke test: synthetic FASTQ → all 6 pipeline stages.

Builds a tiny experiment with 3 samples + ~10 barcodes + 1 near-variant pair
(to exercise bowtie2 global CC) + 1 spike-in. Verifies each stage's outputs.
"""

import shutil
import sys
from pathlib import Path

PROJ_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ_ROOT))

import gzip
import random
from textwrap import dedent

from mapseq2.config import Config, load_sample_sheet
from mapseq2.pipeline import (
    run_cluster,
    run_filter,
    run_global_collapse,
    run_local_collapse,
    run_matrix,
    run_preprocess,
)


BC_LEN = 32
UMI_LEN = 12
SSI_LEN = 8

PREFIXES = [
    "GCAGCAGCAGCAGCAG",
    "TGCATGCATGCATGCA",
    "AATTCCGGAATTCCGG",
    "GGCCAAGGCCAAGGCC",
    "ACGTACGTACGTACGT",
    "TTAACCGGTTAACCGG",
    "CCGGAATTCCGGAATT",
    "AGCTAGCTAGCTAGCT",
    "TCGATCGATCGATCGA",
    "CATGCATGCATGCATG",
    "GTACGTACGTACGTAC",
    "TAGCTAGCTAGCTAGC",
]

def make_bc(prefix: str, body: str, tail: str = "TT") -> str:
    pad = BC_LEN - len(prefix) - len(tail)
    bc = prefix + body[:pad].ljust(pad, "A") + tail
    assert len(bc) == BC_LEN
    return bc


REAL_BARCODES = [make_bc(p, "AAAAAAAAAAAAAA") for p in PREFIXES[:10]]
# near-variant of bc[0]: same 16bp seed, 1 mismatch at position 20
BC0_VARIANT = REAL_BARCODES[0][:19] + "C" + REAL_BARCODES[0][20:]
assert BC0_VARIANT != REAL_BARCODES[0]
assert BC0_VARIANT[:16] == REAL_BARCODES[0][:16]

SPIKE_BC = make_bc(PREFIXES[11], "AAAAAAAA", tail="CGTCAGTC")  # 16+8+8 = 32, ends TC

SSI_BY_SAMPLE = {
    "1": "ACGTACGT",
    "2": "TGCATGCA",
    "3": "GATCGATC",
}


def random_umi(rng: random.Random) -> str:
    return "".join(rng.choices("ACGT", k=UMI_LEN))


def make_fastq(path: Path, records: list[tuple[str, str]]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as f:
        for i, (seq, qual) in enumerate(records):
            f.write(f"@r{i}\n{seq}\n+\n{qual}\n")


def build_reads(rng: random.Random, sample_id: str, bc_counts: dict[str, int]):
    r1_records: list[tuple[str, str]] = []
    r2_records: list[tuple[str, str]] = []
    ssi = SSI_BY_SAMPLE[sample_id]
    for bc, count in bc_counts.items():
        for _ in range(count):
            umi = random_umi(rng)
            r1 = bc
            r2 = umi + ssi
            r1_records.append((r1, "I" * len(r1)))
            r2_records.append((r2, "I" * len(r2)))
    rng.shuffle(r1_records)
    rng.shuffle(r2_records)
    return r1_records, r2_records


def main():
    rng = random.Random(0)
    work = PROJ_ROOT / "tests" / "_e2e_work"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    raw = work / "raw"
    out = work / "out"

    # Sample 1 (source): all 10 real + variant + spike
    # Sample 2 (target, sort=1): some real + spike
    # Sample 3 (target, sort=2): some real + variant + spike
    plans = {
        "1": {bc: 30 for bc in REAL_BARCODES}
            | {BC0_VARIANT: 8, SPIKE_BC: 20},
        "2": {
            REAL_BARCODES[0]: 25, REAL_BARCODES[1]: 22, REAL_BARCODES[2]: 18,
            REAL_BARCODES[3]: 14, REAL_BARCODES[4]: 10, REAL_BARCODES[5]: 7,
            SPIKE_BC: 30,
        },
        "3": {
            REAL_BARCODES[2]: 20, REAL_BARCODES[3]: 18, REAL_BARCODES[4]: 16,
            REAL_BARCODES[5]: 12, REAL_BARCODES[6]: 10, REAL_BARCODES[7]: 8,
            REAL_BARCODES[8]: 6, REAL_BARCODES[9]: 5,
            BC0_VARIANT: 15, SPIKE_BC: 25,
        },
    }
    for sid, bc_counts in plans.items():
        r1, r2 = build_reads(rng, sid, bc_counts)
        make_fastq(raw / f"TEST-{sid}_S{sid}_L001_R1_001.fastq.gz", r1)
        make_fastq(raw / f"TEST-{sid}_S{sid}_L001_R2_001.fastq.gz", r2)

    sample_sheet = work / "sample_sheet.csv"
    sample_sheet.write_text(dedent(
        f"""\
        sample_id,region,side,role,ssi_full_sequence,umi_threshold,sort_order
        1,RegionInj,Ip,source,{SSI_BY_SAMPLE['1']}NNNNNNNNNNNN,1,
        2,RegionA,Ip,target,{SSI_BY_SAMPLE['2']}NNNNNNNNNNNN,1,1
        3,RegionB,Ip,target,{SSI_BY_SAMPLE['3']}NNNNNNNNNNNN,1,2
        """
    ))

    cfg_path = work / "config.yaml"
    cfg_path.write_text(dedent(
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
          extra_flags: ["-a"]
        filtering:
          source_threshold_umi: 2
          proj_threshold_umi: 2
          cell_body_threshold: null
        quality_filter:
          homopolymer_min_run: 0
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

    cfg = Config.from_yaml(cfg_path)
    samples = load_sample_sheet(cfg.input.sample_sheet)
    print(f"Loaded config + {len(samples)} samples")

    print("\n=== Stage 01: preprocess ===")
    run_preprocess(cfg)
    s1_quick = out / "01_preprocess" / "1_quickout.tsv"
    assert s1_quick.exists()
    lines = [l for l in s1_quick.read_text().splitlines() if l.strip()]
    print(f"sample 1 quickout: {len(lines)} barcodes")
    assert len(lines) >= 11, f"expected >=11 (10 real + variant + spike), got {len(lines)}"

    print("\n=== Stage 02: local collapse ===")
    run_local_collapse(cfg)
    s1_local = out / "02_local_collapse" / "1_local_collapsed.tsv"
    assert s1_local.exists()
    lines = [l for l in s1_local.read_text().splitlines() if l.strip()]
    print(f"sample 1 local-collapsed: {len(lines)} barcodes")

    print("\n=== Stage 03: global collapse ===")
    run_global_collapse(cfg)
    gl = out / "03_global_collapse"
    assert (gl / "pooled_seqs.tsv").exists()
    assert (gl / "global_labels.npy").exists()
    import numpy as np
    labels = np.load(gl / "global_labels.npy")
    pooled = (gl / "pooled_seqs.tsv").read_text().strip().splitlines()[1:]
    print(f"global pool: {len(pooled)} sequences → {len(set(labels))} CCs")
    bc0_idx = next(
        i for i, line in enumerate(pooled) if line.split("\t")[1] == REAL_BARCODES[0]
    )
    bcv_idx = next(
        i for i, line in enumerate(pooled) if line.split("\t")[1] == BC0_VARIANT
    )
    assert labels[bc0_idx] == labels[bcv_idx], (
        f"BC0 and BC0_variant should be in same global CC (labels: "
        f"{labels[bc0_idx]} vs {labels[bcv_idx]})"
    )
    print(f"  BC0 + BC0_variant collapsed to same global CC: ✓")

    print("\n=== Stage 04: matrix ===")
    run_matrix(cfg)
    m4 = out / "04_matrix"
    for fn in ["barcode_matrix.tsv", "real_matrix.tsv", "spike_matrix.tsv", "normalized_matrix.tsv"]:
        assert (m4 / fn).exists(), f"missing {fn}"
    import pandas as pd
    real_df = pd.read_csv(m4 / "real_matrix.tsv", sep="\t", index_col="barcode")
    spike_df = pd.read_csv(m4 / "spike_matrix.tsv", sep="\t", index_col="barcode")
    print(f"real matrix: {real_df.shape}, spike: {spike_df.shape}")
    assert real_df.shape[1] == 3
    assert spike_df.shape[0] >= 1, "expected at least one spike-in row"

    print("\n=== Stage 05: filter ===")
    run_filter(cfg)
    m5 = out / "05_filtered"
    for fn in ["B.tsv", "Bnorm.tsv", "B_tar.tsv", "Bnorm_tar.tsv"]:
        assert (m5 / fn).exists(), f"missing {fn}"
    btar = pd.read_csv(m5 / "B_tar.tsv", sep="\t", index_col="barcode")
    print(f"B_tar after filter: {btar.shape}")

    print("\n=== Stage 06: cluster ===")
    run_cluster(cfg)
    m6 = out / "06_cluster"
    if (m6 / "cluster_assignments.tsv").exists():
        ca = pd.read_csv(m6 / "cluster_assignments.tsv", sep="\t")
        print(f"clustering: {len(ca)} barcodes → {ca['cluster'].nunique()} clusters")
        assert (m6 / "heatmap.png").exists()
        print(f"heatmap saved to {m6 / 'heatmap.png'}")
    else:
        print(f"clustering skipped (insufficient barcodes); stage dir exists: {m6.exists()}")

    print("\nALL E2E SMOKE TESTS PASSED")


if __name__ == "__main__":
    main()
