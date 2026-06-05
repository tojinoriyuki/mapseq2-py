#!/bin/bash
# Bash preprocess for MAPseq2-py: ports Hyopil bash_processing_v2.txt stages 1-3.
# Closely mirrors Hyopil's structure (sequential pigz reads, intermediate ext4 files)
# to avoid the parallel-process-substitution back-pressure issues that surface on drvfs.
#
# Required env vars (set by mapseq2/preprocess_bash.py):
#   PREFIX            experiment prefix (matches FASTQ filenames)
#   NSAMPLES          highest sample id (loop runs 1..NSAMPLES)
#   OUT_DIR           output directory ($OUT_DIR/intermediate, $OUT_DIR/thresholds created)
#   RAW_FASTQ_DIR     directory containing FASTQ files
#   SSI_FILE          path to SSI.txt (one SSI per line, line N = sample N)
#   YY_REGEX          library YY filter regex
#   THRESHOLD_ARRAY   space-separated thresholds, leading 0 placeholder
#   BARCODE_LEN       barcode length (typ 32)
#   UMI_LEN           UMI length (typ 12)
#   SSI_LEN           SSI length (typ 8)
#   SSI_MATCH_LEN     SSI match length (typ 8)
#   FASTQ_R1_PATTERN  pattern with {prefix} and {sample_id} placeholders
#   FASTQ_R2_PATTERN  same
#
# Optional:
#   SORT_PARALLEL     threads for sort (default 8)
#   SORT_BUF          sort -S buffer size (default 2G)
#   PIGZ_THREADS      pigz threads (default 8; sequential single-pigz-at-a-time)
#   TMPDIR            sort temp dir (default /tmp)
#   STAGE_DIR         ext4 staging dir for intermediate per-sample files
#                     (default $HOME/scratch_mapseq2_stage)

set -eu
shopt -s nullglob

: "${SORT_PARALLEL:=8}"
: "${SORT_BUF:=2G}"
: "${PIGZ_THREADS:=8}"
: "${TMPDIR:=/tmp}"
: "${STAGE_DIR:=$HOME/scratch_mapseq2_stage}"

mkdir -p "$OUT_DIR/intermediate" "$OUT_DIR/thresholds" "$TMPDIR" "$STAGE_DIR"

threshold=($THRESHOLD_ARRAY)

bu_len=$((BARCODE_LEN + UMI_LEN))
r2_extract=$((UMI_LEN + SSI_LEN))
ssi_pos_start=$((bu_len + 1))
ssi_pos_end=$((bu_len + SSI_MATCH_LEN))

echo "== mapseq2 bash preprocess =="
echo "  PREFIX=$PREFIX, NSAMPLES=$NSAMPLES"
echo "  BU_LEN=$bu_len (barcode $BARCODE_LEN + UMI $UMI_LEN)"
echo "  SSI position in 52-mer: $ssi_pos_start-$ssi_pos_end"
echo "  YY regex: $YY_REGEX"
echo "  SORT_PARALLEL=$SORT_PARALLEL, SORT_BUF=$SORT_BUF, PIGZ_THREADS=$PIGZ_THREADS"
echo "  TMPDIR=$TMPDIR  STAGE_DIR=$STAGE_DIR"
echo "  raw_fastq_dir=$RAW_FASTQ_DIR"
echo

T_GLOBAL=$(date +%s)

for i in $(seq 1 $NSAMPLES); do
  R1_PAT=$(echo "$FASTQ_R1_PATTERN" | sed "s|{prefix}|$PREFIX|g; s|{sample_id}|$i|g")
  R2_PAT=$(echo "$FASTQ_R2_PATTERN" | sed "s|{prefix}|$PREFIX|g; s|{sample_id}|$i|g")
  R1_ARR=($RAW_FASTQ_DIR/$R1_PAT)
  R2_ARR=($RAW_FASTQ_DIR/$R2_PAT)
  R1_FILES="${R1_ARR[*]:-}"
  R2_FILES="${R2_ARR[*]:-}"

  if [ -z "$R1_FILES" ] || [ -z "$R2_FILES" ]; then
    echo "[$i/$NSAMPLES] $(date +%H:%M:%S) NO FASTQ, skipping"
    continue
  fi

  ssi=$(awk -v n=$i 'NR==n {print; exit}' "$SSI_FILE" | tr -d ' \n')
  if [ -z "$ssi" ]; then
    echo "[$i/$NSAMPLES] $(date +%H:%M:%S) NO SSI on line $i, skipping"
    continue
  fi
  ssi_prefix=$(echo "$ssi" | cut -b 1-$SSI_MATCH_LEN)

  T0=$(date +%s)
  echo "[$i/$NSAMPLES] $(date +%H:%M:%S) start (SSI=$ssi_prefix)"

  # === Stage 0: copy FASTQ from drvfs (/mnt/d) to local ext4 stage_dir ===
  # Reading drvfs through small-block tools (cat default 128KB, pigz 16KB) is
  # extremely slow due to drvfs per-syscall overhead. dd with bs=64M reaches
  # ~180 MB/s on the same drive. So front-load a bulk dd copy here.
  R1_LOCAL="$STAGE_DIR/${PREFIX}_${i}_R1.fastq.gz"
  R2_LOCAL="$STAGE_DIR/${PREFIX}_${i}_R2.fastq.gz"
  rm -f "$R1_LOCAL" "$R2_LOCAL"
  for f in $R1_FILES; do dd if="$f" bs=64M status=none; done > "$R1_LOCAL"
  for f in $R2_FILES; do dd if="$f" bs=64M status=none; done > "$R2_LOCAL"
  T_CP=$(date +%s)
  echo "[$i/$NSAMPLES] $(date +%H:%M:%S)   cp to ext4 in $((T_CP-T0))s ($(du -h $R1_LOCAL | cut -f1) + $(du -h $R2_LOCAL | cut -f1))"

  R1_STRIP="$STAGE_DIR/${PREFIX}_${i}_R1_stripped.txt"
  PE_FILE="$STAGE_DIR/${PREFIX}_BC${i}_PE.txt"
  BU_OUT="$OUT_DIR/intermediate/${PREFIX}_BC${i}_BU.txt"
  SSI_CHK="$OUT_DIR/intermediate/${PREFIX}_BC${i}_SSIcheck.txt"
  QO="$OUT_DIR/thresholds/${PREFIX}_BC${i}_quickout.txt"

  # === Stage 1: extract R1 32-mers (now from ext4, fast pigz) ===
  pigz -p $PIGZ_THREADS -dc "$R1_LOCAL" \
    | awk 'NR%4==2' \
    | cut -b 1-$BARCODE_LEN \
    > "$R1_STRIP"
  rm -f "$R1_LOCAL"
  T1=$(date +%s)
  echo "[$i/$NSAMPLES] $(date +%H:%M:%S)   R1 stripped in $((T1-T_CP))s ($(wc -l < $R1_STRIP) reads)"

  # === Stage 1b: paste with R2 20-mer + N filter + (optional) YY filter ===
  if [ -n "$YY_REGEX" ]; then
    pigz -p $PIGZ_THREADS -dc "$R2_LOCAL" \
      | awk 'NR%4==2' \
      | cut -b 1-$r2_extract \
      | paste -d '' "$R1_STRIP" - \
      | grep -v N \
      | grep -E "$YY_REGEX" \
      > "$PE_FILE"
  else
    pigz -p $PIGZ_THREADS -dc "$R2_LOCAL" \
      | awk 'NR%4==2' \
      | cut -b 1-$r2_extract \
      | paste -d '' "$R1_STRIP" - \
      | grep -v N \
      > "$PE_FILE"
  fi
  rm -f "$R1_STRIP" "$R2_LOCAL"
  T2=$(date +%s)
  echo "[$i/$NSAMPLES] $(date +%H:%M:%S)   PE built in $((T2-T1))s ($(wc -l < $PE_FILE) reads)"

  # === Stage 2: SSI distribution ===
  cut -b ${ssi_pos_start}-${ssi_pos_end} "$PE_FILE" \
    | sort -S $SORT_BUF --parallel=$SORT_PARALLEL -T "$TMPDIR" \
    | uniq -c \
    | sort -nr -S $SORT_BUF --parallel=$SORT_PARALLEL -T "$TMPDIR" \
    > "$SSI_CHK"
  T3=$(date +%s)
  echo "[$i/$NSAMPLES] $(date +%H:%M:%S)   SSI check in $((T3-T2))s ($(wc -l < $SSI_CHK) unique SSIs)"

  # === Stage 3: SSI demux + BU collapse (mirrors Hyopil step 2) ===
  grep "${ssi_prefix}\$" "$PE_FILE" \
    | cut -b 1-$bu_len \
    | sort -S $SORT_BUF --parallel=$SORT_PARALLEL -T "$TMPDIR" \
    | uniq -c \
    | sort -nr -S $SORT_BUF --parallel=$SORT_PARALLEL -T "$TMPDIR" \
    > "$BU_OUT"
  rm -f "$PE_FILE"
  T4=$(date +%s)
  echo "[$i/$NSAMPLES] $(date +%H:%M:%S)   BU built in $((T4-T3))s ($(wc -l < $BU_OUT) unique BU mers)"

  # === Stage 4: UMI threshold + 32mer collapse (mirrors Hyopil step 3) ===
  j=${threshold[$i]:-1}
  if [ "$j" != "1" ]; then
    awk -v j=$j -v L=$BARCODE_LEN '$1 >= j { print substr($2, 1, L) }' "$BU_OUT" \
      | sort -S $SORT_BUF --parallel=$SORT_PARALLEL -T "$TMPDIR" \
      | uniq -c \
      | sort -nr -S $SORT_BUF --parallel=$SORT_PARALLEL -T "$TMPDIR" \
      > "$QO"
  else
    awk -v L=$BARCODE_LEN '{ print substr($2, 1, L) }' "$BU_OUT" \
      | sort -S $SORT_BUF --parallel=$SORT_PARALLEL -T "$TMPDIR" \
      | uniq -c \
      | sort -nr -S $SORT_BUF --parallel=$SORT_PARALLEL -T "$TMPDIR" \
      > "$QO"
  fi
  T5=$(date +%s)
  echo "[$i/$NSAMPLES] $(date +%H:%M:%S) DONE total $((T5-T0))s (j=$j, quickout=$(wc -l < $QO))"
done

T_END=$(date +%s)
echo
echo "== preprocess all samples done in $((T_END-T_GLOBAL))s =="
