#!/usr/bin/env bash
# =============================================================================================
# PASSO 7 - Peak calling con MACS3 (v3.0.2) in modalita' paired-end (BAMPE), IP contro il proprio Input.
#   narrow : -f BAMPE -g hs -q 0.01
#   broad  : -f BAMPE -g hs --broad --broad-cutoff 0.1
# Le stesse impostazioni sono state usate per il dataset pubblico UT/KCl (GSE152062).
#
# L'Input di ogni IP si trova dal nome: <ID>_IP_<cond>_<rep>_<S> -> *_Input_<cond>_<rep>_*.sorted.bam
#
# Uso:  BAM_DIR=results/chipseq/bam OUT=results/chipseq/peaks bash 07_macs3_peak_calling.sh
# Output: $OUT/narrow/<IP>_narrow_peaks.narrowPeak, $OUT/broad/<IP>_broad_peaks.broadPeak
# =============================================================================================
set -euo pipefail

BAM_DIR="${BAM_DIR:-results/chipseq/bam}"
OUT="${OUT:-results/chipseq/peaks}"
mkdir -p "$OUT/narrow" "$OUT/broad"

for ip in "$BAM_DIR"/*_IP_*.sorted.bam; do
  base=$(basename "$ip" .sorted.bam)                               # es. S61454_IP_Amr_1_S33
  key=$(echo "$base" | sed -E 's/^S[0-9]+_IP_(.+)_S[0-9]+$/\1/')   # es. Amr_1
  ctrl=$(ls "$BAM_DIR"/*_Input_"${key}"_*.sorted.bam 2>/dev/null | head -1 || true)
  [[ -n "$ctrl" ]] || { echo "ATTENZIONE: nessun Input per $base, salto"; continue; }
  echo ">>> $base  (controllo: $(basename "$ctrl"))"

  [[ -s "$OUT/narrow/${base}_narrow_peaks.narrowPeak" ]] || \
    macs3 callpeak -t "$ip" -c "$ctrl" -f BAMPE -g hs -n "${base}_narrow" \
      --outdir "$OUT/narrow" -q 0.01

  [[ -s "$OUT/broad/${base}_broad_peaks.broadPeak" ]] || \
    macs3 callpeak -t "$ip" -c "$ctrl" -f BAMPE -g hs -n "${base}_broad" \
      --outdir "$OUT/broad" --broad --broad-cutoff 0.1
done

wc -l "$OUT"/narrow/*_peaks.narrowPeak "$OUT"/broad/*_peaks.broadPeak 2>/dev/null || true
