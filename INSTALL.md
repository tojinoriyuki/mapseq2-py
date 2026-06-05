# mapseq2-py インストールマニュアル

`mapseq2-py` のインストール方法を OS 別に説明します。
使い方は [MANUAL.md](MANUAL.md)、設定ファイルの書き方は
`configs/example_config.yaml` を参照してください。

---

## 0. 動作環境

| OS | 動く？ | 方法 |
|---|---|---|
| **Windows 10/11 ネイティブ** | × | `bowtie2` が無いので不可 |
| **Windows 10/11 + WSL2 (Ubuntu)** | ◯ | `install.ps1` |
| **Linux (Ubuntu / Debian / RHEL / Fedora)** | ◯ | `install.sh` |
| **macOS (Intel / Apple Silicon)** | ◯ | `install.sh` |

### Windows ユーザーへの注意

ネイティブ Windows では動きません。`mapseq2` の内部処理は

- bash + pigz / awk / sort (preprocess 段階)
- bowtie2 (CC 段階; Linux/macOS のみ提供)

を必要とするため、Windows では必ず **WSL2 (Windows Subsystem for
Linux)** の中の Ubuntu で動かします。`install.ps1` はそのセットアップ
を自動化しているだけで、走っているのは Linux です。

```
Windows OS
  └─ WSL2 (Linux カーネル)
       └─ Ubuntu
            └─ conda env "mapseq"
                 ├─ bowtie2
                 ├─ pigz
                 └─ Python + mapseq2-py    ← ここで動く
```

ファイルアクセス: WSL 側からは Windows の `C:\Users\you\...` を
`/mnt/c/Users/you/...` として読めます。逆に WSL 内の `~/...` は
エクスプローラから `\\wsl$\Ubuntu\home\you\...` で読めます。
ただし `/mnt/c`, `/mnt/d` は drvfs 経由で遅い(~10×遅延)ので、
大きな FASTQ は WSL 内 ext4 (`~/mapseq_data/...`) にコピーした方が
preprocess は速くなります。

---

## 1. ハードウェア要件

| 項目 | 目安 |
|---|---|
| ディスク | conda 環境 ~3 GB + FASTQ 1サンプルあたり ~4 GB 作業領域 |
| RAM | 16 GB 推奨（8 GB でも 27サンプル run は通る） |
| CPU | コア数が多いほど preprocess + bowtie2 が速い (8 コア以上推奨) |
| ネット | 初回インストールで Miniforge ~80 MB + conda パッケージ ~2.5 GB DL |

---

## 2. 配布物の準備（送り側）

別 PC にインストールする場合の事前準備:

```powershell
# Windows (送り側) — zip にまとめる
Compress-Archive -Path C:\Users\rocke\MAPseq\mapseq2-py `
                 -DestinationPath mapseq2-py.zip
```

```bash
# Linux/macOS (送り側) — tar にまとめる
tar -czf mapseq2-py.tar.gz -C ~/MAPseq mapseq2-py
```

zip / tar に含まれる必要があるもの:

```
mapseq2-py/
├── mapseq2/                       # Pythonパッケージ本体
│   ├── bash_wrappers/preprocess.sh
│   └── *.py
├── configs/
│   ├── example_config.yaml         # configテンプレ
│   ├── example_samplesheet.csv
│   └── samplesheet_mapseq2-1.csv   # MAPseq2-1 SRA用 (参考)
├── tests/                           # 単体テスト・検証スクリプト
├── reference/                       # Hyopil原版 bash の参考
├── pyproject.toml
├── environment.yml
├── install.sh
├── install.ps1
├── README.md
├── INSTALL.md  (このファイル)
└── MANUAL.md
```

含めなくていいもの: `out/`, `raw/`, `__pycache__/`, `.egg-info/`,
個別の解析結果フォルダ。

---

## 3. インストール手順

### A. Windows (WSL2)

#### A-1. 一度きり: WSL2 + Ubuntu のセットアップ

**管理者 PowerShell** で:

```powershell
wsl --install -d Ubuntu
```

→ 再起動 → Ubuntu 初回起動でユーザー名 / パスワードを設定。
(設定したパスワードは `sudo` で使うので忘れないこと)

#### A-2. mapseq2-py インストール

zip を展開してから、**通常の PowerShell** (管理者不要) で:

```powershell
# 受け取った zip を展開 (例: Documents 配下に置く)
Expand-Archive $HOME\Downloads\mapseq2-py.zip -DestinationPath $HOME\MAPseq
cd $HOME\MAPseq\mapseq2-py
.\install.ps1
```

`install.ps1` が中で何をしているか:

1. WSL2 と Ubuntu が入っているか確認
2. WSL の `install.sh` を起動
3. conda が無ければ Miniforge を `~/miniforge3` に自動 DL/展開
4. `environment.yml` から conda env `mapseq` を作成
   (bowtie2 / pigz / Python 依存をまとめて入れる)
5. `mapseq2` CLI を editable モードでインストール
6. smoke test (`bowtie2 --version`, `pigz --version`, `mapseq2 --help`)

所要時間: **5〜15 分** (ネット速度依存、初回は約 2.5 GB DL)

#### A-3. 動作確認

スタートメニューから **Ubuntu** アプリを開いて:

```bash
conda activate mapseq
mapseq2 --help
bowtie2 --version
```

`mapseq2: command not found` が出る場合は、新しい Ubuntu シェルを
開き直してから再度試してください (conda init が PATH を更新します)。

#### A-4. カスタム env 名にする場合

```powershell
$env:ENV_NAME = "mapseq_v2"; .\install.ps1
```

---

### B. Linux / macOS

```bash
# tar を展開
tar -xzf mapseq2-py.tar.gz
cd mapseq2-py

# インストール実行
bash install.sh
```

`install.sh` が中で何をしているか:

1. conda が PATH にあるか確認 (無ければ Miniforge を `$HOME/miniforge3` に自動 DL)
2. `environment.yml` から conda env `mapseq` を作成 (既存なら `update --prune`)
3. smoke test

#### B-1. 動作確認

```bash
conda activate mapseq
mapseq2 --help
bowtie2 --version
```

#### B-2. カスタム env 名

```bash
ENV_NAME=mapseq_v2 bash install.sh
```

---

## 4. インストール後の日常ワークフロー

### Windows ユーザー

新しいシェルを開くたびに WSL → conda activate が必要:

```
1. スタートメニュー → "Ubuntu" を開く
2. プロンプトで:
       conda activate mapseq
       cd /mnt/c/Users/you/MAPseq/my_experiment   # またはWSL内のディレクトリ
       mapseq2 -v run --config config.yaml
```

VS Code から使う場合:
- 拡張機能 `Remote - WSL` をインストール
- VS Code 左下の `><` → "Connect to WSL"
- WSL 内のフォルダを開く → 統合ターミナルが Ubuntu bash になる

### Linux / macOS ユーザー

```bash
conda activate mapseq
cd ~/my_experiment
mapseq2 -v run --config config.yaml
```

### 最小実行例 (どの OS でも共通)

```bash
conda activate mapseq

# 解析ごとに作業ディレクトリを作る
mkdir ~/my_experiment && cd ~/my_experiment

# テンプレを取ってきて編集
cp ~/MAPseq/mapseq2-py/configs/example_config.yaml      config.yaml
cp ~/MAPseq/mapseq2-py/configs/example_samplesheet.csv sample_sheet.csv

# config.yaml で編集する項目:
#   project.name:        FASTQファイル名の prefix (例: "MyExp-A")
#   input.raw_fastq_dir: FASTQの場所 (例: "./raw")
#   input.sample_sheet:  "./sample_sheet.csv"
#   library_filter.regex: ライブラリのアンカー配列に合わせる
#   filtering.source_threshold_umi, proj_threshold_umi: UMIカットオフ

# sample_sheet.csv で各サンプルの:
#   region, side, role, ssi_full_sequence, umi_threshold を記入

# 実行
mapseq2 -v run --config config.yaml
```

出力は `<output_dir>/06_cluster/heatmap.png`。

詳しいパラメータ調整は [MANUAL.md](MANUAL.md) 参照。

---

## 5. アップデート

`mapseq2-py` の最新版を入手したとき:

```bash
conda activate mapseq
cd path/to/mapseq2-py
pip install -e .                                     # コード変更を反映
# environment.yml が更新されている場合のみ:
conda env update -n mapseq -f environment.yml --prune
```

---

## 6. アンインストール

```bash
conda env remove -n mapseq
rm -rf path/to/mapseq2-py
```

WSL ごと消したい場合 (Windows):

```powershell
wsl --unregister Ubuntu
```

(注: WSL 内のすべてのファイルが消えます)

---

## 7. トラブルシュート

| 症状 | 原因 / 対処 |
|---|---|
| `wsl: command not found` (Windows) | WSL2 未インストール。**管理者** PowerShell で `wsl --install -d Ubuntu` → 再起動 |
| `Ubuntu` がスタートメニューに無い | `wsl --install -d Ubuntu` 後、初回起動が必要 (Microsoft Store からも入手可) |
| `conda: command not found` | install.sh が conda init した直後は古いシェルから見えない。**新しいシェルを開き直す** |
| `mapseq2: command not found` | `conda activate mapseq` を忘れている |
| `bowtie2: command not found` | env が壊れている。`conda env remove -n mapseq && bash install.sh` で作り直し |
| `Permission denied` (`.sh`) | `chmod +x install.sh` か `bash install.sh` で実行 |
| C: 容量不足 | PowerShell の `Get-PSDrive C` で実容量を確認 (**WSL の `df` は嘘** — 仮想 VHDX 容量を返す) |
| WSL が動かない (BIOS で仮想化無効) | BIOS で Intel VT-x / AMD-V を有効化。Hyper-V 機能も `Windows の機能の有効化` で ON |
| install.ps1 が `実行ポリシー` で止まる | `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned` で許可 (1回だけ) |
| プロキシ環境で conda が DL 失敗 | `~/.condarc` に `proxy_servers:` 設定を追加 |
| FASTQ が `/mnt/c` 配下で preprocess が極端に遅い | drvfs オーバーヘッド。FASTQ を WSL 内 ext4 (`~/mapseq_data/raw/`) にコピーしてから実行 |
| `WSL の df は 1TB free と言うのに install が "no space" で失敗` | WSL の df は VHDX 仮想容量。`Get-PSDrive C` で Windows 側の実 free を確認 |

---

## 8. まとめ (1ページ要約)

```
# Windows ───────────────────────────────────────────
# 一度だけ (管理者 PowerShell):
wsl --install -d Ubuntu                       # 再起動 + Ubuntu 初回設定

# zip を展開して通常の PowerShell で:
cd path\to\mapseq2-py
.\install.ps1                                 # 5-15 分

# 以後、Ubuntu 起動 → 解析:
conda activate mapseq
mapseq2 -v run --config config.yaml

# Linux / macOS ─────────────────────────────────────
cd path/to/mapseq2-py
bash install.sh                               # 5-15 分

conda activate mapseq
mapseq2 -v run --config config.yaml
```

詳しい使い方は [MANUAL.md](MANUAL.md) を見てください。
