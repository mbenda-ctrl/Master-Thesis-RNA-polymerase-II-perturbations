#!/usr/bin/env bash
# =============================================================================================
# PASSO 1 - Controllo qualita' delle read grezze (paired-end): FastQC su R1 e R2 + MultiQC.
#
# Per ogni campione FastQC analizza separatamente R1 e R2. MultiQC aggrega i report in tre file:
# tutti i campioni, solo IP e solo Input. Il file multiqc_config.yaml unisce R1 e R2 dello stesso
# campione in una sola riga della tabella generale (modalita' paired-end).
#
# Opzionale (RUN_FASTP=1): report fastp in sola lettura (nessun filtro o trimming), come nel
# controllo qualita' originale; anche questi report finiscono in MultiQC.
#
# Input : FASTQ_DIR con file <campione>_R1_001.fastq.gz e <campione>_R2_001.fastq.gz
#         (cerca anche nelle sottocartelle, es. Sample_<campione>/)
# Output: $OUT/fastqc/*, $OUT/fastp/* (opzionale), $OUT/multiqc_{all,IP,Input}.html
#
# Uso:  FASTQ_DIR=data/fastq OUT=results/chipseq/qc THREADS=4 bash 01_fastqc_multiqc.sh
# =============================================================================================
set -euo pipefail

FASTQ_DIR="${FASTQ_DIR:-data/fastq}"
OUT="${OUT:-results/chipseq/qc}"
THREADS="${THREADS:-4}"
RUN_FASTP="${RUN_FASTP:-0}"
SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$OUT/fastqc" "$OUT/fastp"

mapfile -t R1S < <(find "$FASTQ_DIR" -name "*_R1_001.fastq.gz" | sort)
[[ ${#R1S[@]} -gt 0 ]] || { echo "ERRORE: nessun *_R1_001.fastq.gz in $FASTQ_DIR"; exit 1; }
echo "Campioni trovati: ${#R1S[@]}"

for r1 in "${R1S[@]}"; do
  r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
  s=$(basename "$r1" _R1_001.fastq.gz)
  [[ -s "$r2" ]] || { echo "ATTENZIONE: manca R2 per $s, salto"; continue; }

  # FastQC: R1 e R2 dello stesso campione
  if [[ ! -s "$OUT/fastqc/${s}_R2_001_fastqc.zip" ]]; then
    echo ">>> FastQC $s"
    fastqc -t "$THREADS" -o "$OUT/fastqc" "$r1" "$r2"
  fi

  # fastp in modalita' paired-end, solo report (come nel QC originale)
  if [[ "$RUN_FASTP" == 1 && ! -s "$OUT/fastp/${s}.json" ]]; then
    echo ">>> fastp $s"
    fastp -i "$r1" -I "$r2" -w "$THREADS" \
          --disable_quality_filtering --disable_length_filtering --disable_adapter_trimming \
          -j "$OUT/fastp/${s}.json" -h "$OUT/fastp/${s}.html" -o /dev/null -O /dev/null
  fi
done

# MultiQC: tutti i campioni, poi IP e Input separati (meno memoria e report piu' leggibili)
CONF="$SCRIPTS/multiqc_config.yaml"
multiqc -f -c "$CONF" -n multiqc_all   -o "$OUT" "$OUT/fastqc" "$OUT/fastp"
multiqc -f -c "$CONF" -n multiqc_IP    -o "$OUT" --ignore "*_Input_*" "$OUT/fastqc" "$OUT/fastp"
multiqc -f -c "$CONF" -n multiqc_Input -o "$OUT" --ignore "*_IP_*"    "$OUT/fastqc" "$OUT/fastp"
echo "Report: $OUT/multiqc_all.html, multiqc_IP.html, multiqc_Input.html"
