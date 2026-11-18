#!/usr/bin/env bash
# =============================================================================================
# PASSO 3 - Allineamento paired-end su GRCh38 con Bowtie2, conversione in BAM, ordinamento e indice.
#
# Parametri ricavati dalle intestazioni (@PG) dei BAM originali:
#   bowtie2 (v2.4.5) -x <indice> -1 R1 -2 R2  -> parametri di default (end-to-end, --sensitive)
#   samtools view -bS                          -> SAM -> BAM, nessun filtro
#   samtools sort                              -> ordinamento per coordinate (samtools v1.19.2)
# Nessuna rimozione dei duplicati PCR.
# Qui i tre passaggi sono collegati con una pipe (stesso risultato, senza il SAM intermedio).
#
# Uso:  FASTQ_DIR=data/fastq BT2_INDEX=reference/hg38/GRCh38_no_alt OUT=results/chipseq/bam \
#       THREADS=8 bash 03_align_bowtie2_pe.sh
# Output: $OUT/<campione>.sorted.bam (+ .bai) e $OUT/logs/<campione>.bowtie2.log
# =============================================================================================
set -euo pipefail

FASTQ_DIR="${FASTQ_DIR:-data/fastq}"
BT2_INDEX="${BT2_INDEX:-reference/hg38/GRCh38_no_alt}"
OUT="${OUT:-results/chipseq/bam}"
THREADS="${THREADS:-8}"

mkdir -p "$OUT/logs"
[[ -s "${BT2_INDEX}.1.bt2" ]] || { echo "ERRORE: indice Bowtie2 non trovato ($BT2_INDEX)"; exit 1; }

mapfile -t R1S < <(find "$FASTQ_DIR" -name "*_R1_001.fastq.gz" | sort)
for r1 in "${R1S[@]}"; do
  r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
  s=$(basename "$r1" _R1_001.fastq.gz)
  bam="$OUT/${s}.sorted.bam"
  [[ -s "$bam.bai" ]] && { echo "gia' fatto: $s"; continue; }
  echo ">>> $s"

  bowtie2 -p "$THREADS" -x "$BT2_INDEX" -1 "$r1" -2 "$r2" 2> "$OUT/logs/${s}.bowtie2.log" \
    | samtools view -@ "$THREADS" -bS - \
    | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index -@ "$THREADS" "$bam"
  tail -1 "$OUT/logs/${s}.bowtie2.log"     # tasso di allineamento complessivo
done
