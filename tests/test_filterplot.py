"""Test filterplot scan + plot on a synthetic barcode matrix.

Builds a real_matrix.tsv with neurons at varying source-UMI counts and varying
projection breadths, runs scan_thresholds + plot, verifies outputs.
"""

import shutil
import sys
from pathlib import Path
from textwrap import dedent

PROJ_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJ_ROOT))

import numpy as np
import pandas as pd

from mapseq2.config import Config
from mapseq2.matrix import save_matrix
from mapseq2.pipeline import run_filterplot


def main():
    work = PROJ_ROOT / "tests" / "_filterplot_work"
    if work.exists():
        shutil.rmtree(work)
    out = work / "out"
    mat_dir = out / "04_matrix"
    mat_dir.mkdir(parents=True)

    rng = np.random.default_rng(0)
    n_neurons = 500
    n_samples = 6  # 1 source + 5 targets

    src_umi = rng.integers(0, 200, size=n_neurons)
    n_targets_per_neuron = np.clip(
        np.round((src_umi / 200) * 5 + rng.normal(0, 0.5, n_neurons)),
        0, 5,
    ).astype(int)

    real_mat = np.zeros((n_neurons, n_samples), dtype=np.int64)
    real_mat[:, 0] = src_umi
    for i in range(n_neurons):
        if n_targets_per_neuron[i] == 0:
            continue
        target_choice = rng.choice(np.arange(1, n_samples), size=n_targets_per_neuron[i], replace=False)
        for t in target_choice:
            real_mat[i, t] = rng.integers(1, 80)

    bcs = [f"BC{i:04d}" + "A" * 28 for i in range(n_neurons)]
    sample_labels = ["1:Inj-Ip:source"] + [f"{i+2}:Tgt{i+1}-Ip:target" for i in range(5)]
    save_matrix(real_mat, bcs, sample_labels, mat_dir / "real_matrix.tsv")

    sample_sheet = work / "sample_sheet.csv"
    rows = ["sample_id,region,side,role,ssi_full_sequence,umi_threshold,sort_order"]
    rows.append("1,Inj,Ip,source,ACGTACGTNNNNNNNNNNNN,1,")
    for i in range(5):
        rows.append(f"{i+2},Tgt{i+1},Ip,target,AAAAAAAA{i}NNNNNNNNNNN,1,{i+1}")
    sample_sheet.write_text("\n".join(rows) + "\n")

    (work / "config.yaml").write_text(dedent(
        f"""\
        project:
          name: TEST
          output_dir: ./out
        input:
          raw_fastq_dir: ./raw
          r1_pattern: "x"
          r2_pattern: "y"
          sample_sheet: ./sample_sheet.csv
        barcode:
          barcode_length: 32
          umi_length: 12
          ssi_length: 8
          ssi_match_length: 8
        library_filter:
          enabled: false
          regex: null
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
          source_threshold_umi: 20
          proj_threshold_umi: 10
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

    cfg = Config.from_yaml(work / "config.yaml")
    run_filterplot(cfg, max_thresh=50, max_targets=5, scan="source", fixed_proj=0)
    run_filterplot(cfg, max_thresh=30, max_targets=5, scan="proj", fixed_source=0)

    fp_dir = mat_dir / "_filterplots"
    assert (fp_dir / "source_threshold_scan.png").exists()
    assert (fp_dir / "source_threshold_scan.tsv").exists()
    assert (fp_dir / "proj_threshold_scan.png").exists()

    df = pd.read_csv(fp_dir / "source_threshold_scan.tsv", sep="\t")
    print(df.head(8).to_string(index=False))
    print("...")
    print(df.tail(3).to_string(index=False))

    assert df["n_kept"].iloc[0] >= df["n_kept"].iloc[-1], "n_kept should monotonically decrease as threshold rises"
    assert df.shape[0] == 50
    print("\nfilterplot test PASSED")


if __name__ == "__main__":
    main()
