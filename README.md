# mapseq2-py

Config-driven Python port of Hyopil Kim's MAPseq2 pipeline. End-to-end
from demultiplexed Illumina FASTQs to a Ward-clustered projection
heatmap.

Reference MATLAB + bash code:
[kebschulllab/MAPseq2-and-POINTseq](https://github.com/kebschulllab/MAPseq2-and-POINTseq).
This package re-implements the v2 bash workflow + MATLAB post-processing
in Python (scipy `connected_components` ≡ MATLAB `graphconncomp`,
validated against the published `.mat` files with `AdjustedRandIndex = 1.0`
on filter + cluster stages).

## Quick links

* **[INSTALL.md](INSTALL.md)** — Windows (WSL2) / Linux / macOS install
* **[MANUAL.md](MANUAL.md)** — full user guide, config reference, threshold tuning
* `configs/example_config.yaml` — annotated config template
* `configs/example_samplesheet.csv` — sample-sheet template
* `configs/samplesheet_mapseq2-1.csv` — sample sheet for SRA SAMN53373349 (MAPseq2-1)

## Install (3 commands)

```bash
# Linux / WSL2 / macOS
git clone https://github.com/<your-username>/mapseq2-py.git
cd mapseq2-py
bash install.sh
```

```powershell
# Windows (PowerShell; assumes WSL2 + Ubuntu)
git clone https://github.com/<your-username>/mapseq2-py.git
cd mapseq2-py
.\install.ps1
```

If you don't have WSL2 yet, see [INSTALL.md](INSTALL.md#a-1-一度きり-wsl2--ubuntu-のセットアップ).

## Run (3 commands)

```bash
conda activate mapseq
cp configs/example_config.yaml my_config.yaml
cp configs/example_samplesheet.csv my_samples.csv
# edit my_config.yaml + my_samples.csv to point at your FASTQs
mapseq2 -v run --config my_config.yaml
```

Heatmap drops at `<output_dir>/06_cluster/heatmap.png`.

## Stages

```
FASTQ ─► preprocess ─► collapse-local ─► collapse-global
                                              │
                                              ▼
heatmap ◄── cluster ◄── filter ◄── matrix
```

Each stage is invocable on its own — see [MANUAL.md](MANUAL.md#stage-by-stage-invocation).

## License

MIT. See [LICENSE](LICENSE).
