#!/usr/bin/env bash
# =============================================================================================
# PASSO 5 - bigWig normalizzati in CPM per IP e Input (deepTools bamCoverage 3.5.6).
#   --normalizeUsing CPM   --binSize 10   --extendReads 200
# I bigWig IP - Input usati per i metagene si ottengono poi con 06_bigwig_IPminusInput.sh.
#
# Uso:  BAM_DIR=results/chipseq/bam OUT=results/chipseq/bigwig_cpm THREADS=4 bash 05_bigwig_cpm.sh
# =============================================================================================
set -euo pipefail

BAM_DIR="${BAM_DIR:-results/chipseq/bam}"
OUT="${OUT:-results/chipseq/bigwig_cpm}"
THREADS="${THREADS:-4}"
mkdir -p "$OUT"

for bam in "$BAM_DIR"/*.sorted.bam; do
  s=$(basename "$bam" .sorted.bam)
  bw="$OUT/${s}.bw"
  [[ -s "$bw" ]] && { echo "gia' fatto: $s"; continue; }
  echo ">>> $s"
  bamCoverage -b "$bam" -o "$bw" -p "$THREADS" \
    --normalizeUsing CPM --binSize 10 --extendReads 200
done
