#!/usr/bin/env bash
# WGS batch runner — process chr2 through chrY sequentially.
# Follows CLAUDE.md §5 (disconnect-resilient). Launched in tmux.
#
# For each chromosome:
#   1. Build position labels (if not present)
#   2. Forward pass on all windows
#   3. Aggregate chunks → per-context summary
#   4. Update status file
#
# Resumable: skips chromosomes with an existing chr{N}_context_summary.csv.

set -uo pipefail

TDIG=/NHNHOME/WORKSPACE/0526040123_A/darejinn/tdig
CHROMS=(chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
STAMP=$(date +%Y%m%d_%H%M%S)
BATCH_LOG=$TDIG/wgs/logs/batch_runner_${STAMP}.log
STATUS_FILE=$TDIG/wgs/logs/batch_status.txt
mkdir -p "$(dirname "$BATCH_LOG")"

log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" | tee -a "$BATCH_LOG"; }

echo "STARTED $(date)" > $STATUS_FILE

source /home/yuhs_seun/miniforge3/etc/profile.d/conda.sh
conda activate gdtr
export HF_HOME=$TDIG/env/hf_cache
# GPU 0 is occupied by user's CLM/DeepSets train.py; use GPU 1 only.
export CUDA_VISIBLE_DEVICES=1

log "== WGS batch runner started =="
log "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES (GPU 1 only)"
log "chromosomes: ${CHROMS[*]}"

for CHR in "${CHROMS[@]}"; do
  SUMMARY=$TDIG/wgs/results/$CHR/${CHR}_context_summary.csv
  if [ -f "$SUMMARY" ]; then
    log "[$CHR] already aggregated ($SUMMARY exists) — skip"
    echo "$CHR SKIPPED $(date)" >> $STATUS_FILE
    continue
  fi

  log "[$CHR] STARTING"
  echo "$CHR STARTED $(date)" >> $STATUS_FILE
  T0=$(date +%s)

  # 1. labels (use v2 with UTR fix)
  LAB_FILE=$TDIG/wgs/data/labels/${CHR}_position_labels.npy
  if [ ! -f "$LAB_FILE" ]; then
    log "[$CHR] building labels …"
    python $TDIG/wgs/scripts/build_chr_position_labels_v2.py --chrom $CHR 2>&1 | tail -20 | tee -a $BATCH_LOG
    if [ ! -f "$LAB_FILE" ]; then
      log "[$CHR] LABEL BUILD FAILED — skip chromosome"
      echo "$CHR LABEL_FAIL $(date)" >> $STATUS_FILE
      continue
    fi
  fi

  # 2. forward
  FWD_LOG=$TDIG/wgs/logs/e1_00_${CHR}_${STAMP}.log
  log "[$CHR] forward pass → $FWD_LOG"
  python $TDIG/wgs/scripts/e1_00_chr_forward.py --chrom $CHR > $FWD_LOG 2>&1
  RC=$?
  if [ $RC -ne 0 ]; then
    log "[$CHR] FORWARD FAILED rc=$RC — inspect $FWD_LOG"
    echo "$CHR FORWARD_FAIL $(date)" >> $STATUS_FILE
    continue
  fi

  # 3. aggregate
  log "[$CHR] aggregate"
  python $TDIG/wgs/scripts/e1_10_aggregate_chr.py --chrom $CHR 2>&1 | tail -30 | tee -a $BATCH_LOG

  T1=$(date +%s)
  DT=$((T1 - T0))
  log "[$CHR] DONE in $((DT / 60)) min"
  echo "$CHR DONE $(date) $((DT / 60))min" >> $STATUS_FILE
done

log "== WGS batch runner COMPLETE =="
echo "COMPLETED $(date)" >> $STATUS_FILE
touch $TDIG/wgs/logs/batch_runner.done
