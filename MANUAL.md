# mapseq2-py User Manual

A config-driven, Python port of Hyopil Kim's MAPseq2 pipeline. Takes
demultiplexed Illumina FASTQs of barcoded brain regions and produces a
spike-normalized projection matrix + Ward-clustered heatmap.

For installation, see [INSTALL.md](INSTALL.md).

## Conceptual overview

```
FASTQ                                    ──┐
  └─ R1: 32 nt barcode + 12 nt UMI         │ stage 01: preprocess
  └─ R2: 12 nt UMI continuation + 8 nt SSI │   FASTQ -> 44-mer BU counts
                                           │   -> 32-mer quickout (UMI ≥ j)
                                         ──┘
quickout (per sample)                    ──┐
  └─ bowtie2 all-vs-all + connected-comp   │ stage 02: collapse-local
  └─ merge near-identical barcodes/sample  │   sequencing-error correction
                                         ──┘
pooled barcodes (all samples)            ──┐
  └─ bowtie2 all-vs-all + CC across pool   │ stage 03: collapse-global
  └─ assign each barcode to a global label │   cross-sample agreement
                                         ──┘
matrix (barcode × sample, raw UMI)       ──┐
  └─ split spike-in CCs (tag CGTCAGTC)     │ stage 04: matrix
  └─ normalize real by unique spike count  │   real_matrix + spike_matrix
                                           │   + normalized_matrix
                                         ──┘
filter on source ≥ N, projection ≥ M     ──┐ stage 05: filter
clean projection matrix B / Bnorm        ──┘

Ward clustering + dendrogram + heatmap   ──┐ stage 06: cluster
                                         ──┘
```

## Quick start (3 commands)

```bash
conda activate mapseq

# 1. Copy and edit the example configs for your experiment
cp configs/example_config.yaml my_config.yaml
cp configs/example_samplesheet.csv my_samples.csv
# (edit paths, sample_id↔SSI mapping, umi_threshold)

# 2. Run end-to-end (preprocess → cluster)
mapseq2 -v run --config my_config.yaml
```

The heatmap lands in `<output_dir>/06_cluster/heatmap.png`.

## Inputs

### 1. FASTQ files

* Paired-end, gzipped (`*_R1_001.fastq.gz` / `*_R2_001.fastq.gz`).
* R1 must be ≥ 32 nt (the random barcode).
* R2 must be ≥ 20 nt: 12 nt UMI continuation + 8 nt sample-specific SSI.
* Filename pattern is configurable (see `input.r1_pattern` /
  `r2_pattern`); placeholders `{prefix}` (= `project.name`) and
  `{sample_id}` are substituted.

### 2. Sample sheet (CSV)

One row per FASTQ sample. Required columns:

| column | meaning |
|---|---|
| `sample_id` | integer 1, 2, 3, … matches `{sample_id}` in FASTQ filename pattern |
| `region` | brain region label (e.g. `MOp`, `STR-r`) — used in column headers |
| `side` | `Ip` (ipsilateral) / `Con` (contralateral) — used in headers |
| `role` | `source` (injection site, ≥1 required) or `target` |
| `ssi_full_sequence` | 8-letter SSI used to demultiplex this sample's reads |
| `umi_threshold` | per-sample UMI count threshold (see "Choosing thresholds" below) |
| `sort_order` | column order in the final heatmap; leave blank for `source` rows |

A complete example is in [configs/example_samplesheet.csv](configs/example_samplesheet.csv).

### 3. Config YAML

Annotated example in [configs/example_config.yaml](configs/example_config.yaml).
The fields you typically edit per experiment:

* `project.name`: the `{prefix}` token in FASTQ filenames (e.g.
  `MAPseq2-1`).
* `project.output_dir`: where pipeline outputs go.
* `input.raw_fastq_dir`: where R1/R2 FASTQs live.
* `input.sample_sheet`: path to your CSV.
* `library_filter.enabled` and `regex`: see **Library anchor filter**.
* `filtering.source_threshold_umi`, `proj_threshold_umi`: hard UMI
  cutoffs at the matrix-filter stage (Hyopil defaults: 20 source, 10
  projection).
* `clustering.num_clusters`: number of Ward-clustering classes (5 in
  the MAPseq2 paper).

### Library anchor filter

The random barcode in MAPseq virus libraries ends with a 2- or 3-letter
**library anchor** at positions 31-32 (or 30-32) of the 32-mer.
Hyopil's "v2" pipeline (DAT-1 library) uses YY = `[TC][TC]`. Some older
or different libraries use a different anchor. The `library_filter`
block controls this:

```yaml
library_filter:
  enabled: true
  regex: '^(.{30}[TC][TC]|.{29}[TC][TC]G)'   # v2 / DAT-1 YY filter
```

* `enabled: false` disables the filter entirely (matches the v1 bash
  that generated the published `MAPseq2-1.mat`).
* Change `regex` if your library uses a different anchor. The regex is
  matched against the 32-mer R1 barcode; anything that doesn't match
  is dropped *before* UMI counting.

## Outputs (per stage)

```
output_dir/
├── 01_preprocess/
│   ├── <id>_BU.tsv         44-mer (32 barcode + 12 UMI) raw counts
│   ├── <id>_quickout.tsv   32-mer after UMI threshold j
│   ├── <id>_ssi_check.png  SSI distribution QC plot
│   ├── <id>_stats.tsv      read-throughput metrics
│   └── _rankplots/         (after `rankplot`) per-sample BU rank plots
├── 02_local_collapse/
│   └── <id>_local_collapsed.tsv   sample's barcodes after CC merge
├── 03_global_collapse/
│   ├── pooled_seqs.tsv     all samples' barcodes pooled + globally CC'd
│   ├── per_sample_counts.tsv  per-(sample,barcode) UMI count
│   └── global_labels.npy
├── 04_matrix/
│   ├── barcode_matrix.tsv     all barcodes × samples (UMI)
│   ├── spike_matrix.tsv       spike-tag barcodes only
│   ├── real_matrix.tsv        non-spike barcodes
│   └── normalized_matrix.tsv  real / unique spike count per sample
├── 05_filtered/
│   ├── B.tsv          source ≥ s_thr, projection ≥ p_thr (all samples)
│   ├── B_tar.tsv      same, but source columns dropped (target-only)
│   ├── Bnorm.tsv      spike-normalized version of B
│   └── Bnorm_tar.tsv  spike-normalized B_tar — **input to clustering**
└── 06_cluster/
    ├── linkage_Z.npy          scipy Ward linkage matrix
    ├── perm.npy               dendrogram leaf order
    ├── labels.npy             cluster assignment per barcode
    ├── sorted_data.npy        Bnorm_tar reordered by `perm`
    ├── cluster_assignments.tsv  barcode → cluster_id mapping
    └── heatmap.png            dendro + cluster bar + heatmap
```

## Stage-by-stage invocation

`mapseq2 run` does all six stages, but you can invoke them
individually — useful when iterating on thresholds:

```bash
mapseq2 preprocess        --config my_config.yaml   # stage 01
mapseq2 rankplot          --config my_config.yaml --suggest --annotate
mapseq2 umi-collapse      --config my_config.yaml   # re-run UMI step only
mapseq2 collapse-local    --config my_config.yaml   # stage 02
mapseq2 collapse-global   --config my_config.yaml   # stage 03
mapseq2 matrix            --config my_config.yaml   # stage 04
mapseq2 filterplot        --config my_config.yaml   # scan source/proj thresholds
mapseq2 filter            --config my_config.yaml   # stage 05
mapseq2 cluster           --config my_config.yaml   # stage 06
```

### preprocess engine

```bash
mapseq2 preprocess --config my_config.yaml --engine bash      # default; fast (pigz/awk/sort, multi-core)
mapseq2 preprocess --config my_config.yaml --engine python    # portable; slow (single-thread Python)
```

The bash engine requires `pigz`, `awk`, `sort`, `grep`, `paste`, `cut`,
`dd` on PATH (all included in the conda env). The Python engine works
anywhere Python does, but is ~30× slower per sample.

## Choosing UMI thresholds (`umi_threshold` column)

For each sample, the UMI threshold *j* keeps only 44-mer BU sequences
that appear *j* or more times. Setting *j* too low keeps PCR/sequencing
noise; too high discards real barcodes.

1. Run `mapseq2 preprocess` once with any threshold (e.g. `j=1`).
2. Run `mapseq2 rankplot --config … --suggest`.
   This produces `01_preprocess/_rankplots/sample_*.png` and a
   `suggested_thresholds.tsv` based on the Kneedle algorithm
   (perpendicular distance from the chord on log-log axes).
3. Inspect the plots — pick the *L*-shaped knee, where the curve
   transitions from a plateau to a steep drop.
4. Edit `umi_threshold` in your sample sheet.
5. Re-run `mapseq2 umi-collapse` (does **not** re-read FASTQ; reuses
   `BU.tsv` from preprocess), then everything downstream.

Hyopil's MAPseq2-1 published thresholds, as reference:

```
sample 1..27:
  15 15 2 15 15 15 15 15 15 10 15 9 15 15 15 15 15 15 15 15 15 15 15 20 20 15 15
```

## Choosing matrix filter thresholds

After `matrix`, scan how `source_threshold` (UMI count at the
injection site) affects the kept-barcode population:

```bash
mapseq2 filterplot --config my_config.yaml --scan source --max-thresh 100
```

Produces `04_matrix/_filterplots/source_threshold_scan.png`. Edit
`filtering.source_threshold_umi` and `proj_threshold_umi` in your
config, then `mapseq2 filter && mapseq2 cluster`.

## Re-using preprocess output between threshold experiments

Stage 01 is the slow stage (FASTQ I/O + decompression). To avoid
re-reading FASTQs when changing **only** `umi_threshold`:

```bash
# … edit umi_threshold in sample sheet …
mapseq2 umi-collapse  --config my_config.yaml   # cheap; reads BU.tsv
mapseq2 collapse-local --config my_config.yaml   # must rerun (input changed)
mapseq2 collapse-global --config my_config.yaml
mapseq2 matrix --config my_config.yaml
mapseq2 filter --config my_config.yaml
mapseq2 cluster --config my_config.yaml
```

To change **only** the filter thresholds (after matrix exists):

```bash
mapseq2 filter --config my_config.yaml
mapseq2 cluster --config my_config.yaml
```

## Validating against Hyopil's `.mat` files

Hyopil's published `barcodematrixMAPseq2-1.mat` (in his
`MAPseq2-and-POINTseq` repo) holds the reference barcode set for the
v1 pipeline. If you want bit-exact agreement:

* Set `library_filter.enabled: false` (v1 has no YY filter).
* Use Hyopil's published `umi_threshold` values.
* Compare your `04_matrix/real_matrix.tsv` against the .mat's
  `refbarcodes` field (see `tests/validate_real_data.py` for the
  loader pattern).

Note: v2 (DAT-1) and v1 (MAPseq2-1) use **different virus libraries**
with different anchor designs. Applying v2's YY filter to MAPseq2-1
reads drops ~60% of valid barcodes (we measured this — see
`tests/validate_real_data.py`).

## Tests

```bash
conda activate mapseq
cd path/to/mapseq2-py
pytest tests/test_smoke.py -x              # end-to-end on tiny fixture
pytest tests/test_spike_homopolymer.py -x  # spike split + HP filter
pytest tests/test_rankplot.py -x           # rankplot matches MATLAB shape
pytest tests/test_filterplot.py -x         # filterplot matches MATLAB shape
```

## Known performance notes

| | typical wall-clock |
|---|---|
| preprocess (bash) for 27 × 4 GB FASTQs | 30–80 min |
| local CC (j=15) | 10–20 min |
| global CC | 5–10 min |
| matrix + filter + cluster | < 2 min |

* Multi-core sort + pigz dominate preprocess time. `bowtie2.threads`
  in `config.yaml` is the main lever for CC stages.
* WSL drvfs (`/mnt/c`, `/mnt/d`) is ~10× slower than WSL ext4 for
  many-small-read workloads. The bash preprocess stages each
  drvfs-resident FASTQ to a local ext4 scratch dir via
  `dd bs=64M`; downstream stages do not benefit from drvfs.

## License / attribution

The reference MATLAB and bash code is by Hyopil Kim
(`kebschulllab/MAPseq2-and-POINTseq`). This package re-implements the
v2 bash workflow + MATLAB post-processing in Python, keeping the
sequencing-error-correction logic numerically identical (scipy
`connected_components` ≡ MATLAB `graphconncomp`, Adjusted Rand Index
1.0 against the published `.mat` for the filter + cluster stages).
