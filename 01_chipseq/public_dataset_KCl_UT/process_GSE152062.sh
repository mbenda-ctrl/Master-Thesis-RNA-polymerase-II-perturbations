#!/usr/bin/env bash
# =============================================================================================
# Dataset pubblico GSE152062 (HEK293T, anti-Pol II ChIP-seq, UT vs KCl, paired-end).
# Usato come riferimento per lo stress osmotico (UT* / KCl).
#
# Passaggi (dalle intestazioni @PG dei BAM):
#   1. download dei FASTQ da SRA (fasterq-dump)
#   2. FastQC + MultiQC                                    (come il passo 1)
#   3. bowtie2 su GRCh38 no_alt (default, paired-end)      (indice del passo 2)
#   4. samtools view -q 30 -F 2820 -> sort -> index        (MAPQ >= 30; esclude non mappate,
#      secondarie, QC-fail, duplicati marcati e supplementari)
#   5. bigWig CPM, MACS3 narrow/broad BAMPE                (passi 5 e 7, stessi parametri)
#
# Input: config/public_GSE152062_runs.tsv (colonne run, sample)
# Uso:   RUNS=config/public_GSE152062_runs.tsv BT2_INDEX=reference/hg38/GRCh38_no_alt \
#        OUT=results/public_KCl_UT THREADS=8 bash process_GSE152062.sh
# =============================================================================================
set -euo pipefail

RUNS="${RUNS:-config/public_GSE152062_runs.tsv}"
BT2_INDEX="${BT2_INDEX:-reference/hg38/GRCh38_no_alt}"
OUT="${OUT:-results/public_KCl_UT}"
THREADS="${THREADS:-8}"
mkdir -p "$OUT"/{fastq,bam,logs,qc}

grep -v '^#' "$RUNS" | tail -n +2 | while IFS=$'\t' read -r run sample; do
  [[ -n "$run" ]] || continue
  echo ">>> $sample ($run)"
  if [[ ! -s "$OUT/fastq/${sample}_R1_001.fastq.gz" ]]; then
    fasterq-dump --split-files -e "$THREADS" -O "$OUT/fastq" "$run"
    gzip -c "$OUT/fastq/${run}_1.fastq" > "$OUT/fastq/${sample}_R1_001.fastq.gz"
    gzip -c "$OUT/fastq/${run}_2.fastq" > "$OUT/fastq/${sample}_R2_001.fastq.gz"
    rm -f "$OUT/fastq/${run}_1.fastq" "$OUT/fastq/${run}_2.fastq"
  fi
  bam="$OUT/bam/${sample}.q30.sorted.bam"
  [[ -s "$bam.bai" ]] && continue
  bowtie2 -p "$THREADS" -x "$BT2_INDEX" --rg-id "$sample" --rg "SM:$sample" \
          -1 "$OUT/fastq/${sample}_R1_001.fastq.gz" -2 "$OUT/fastq/${sample}_R2_001.fastq.gz" \
          2> "$OUT/logs/${sample}.bowtie2.log" \
    | samtools view -@ "$THREADS" -bh -q 30 -F 2820 - \
    | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index "$bam"
done

SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FASTQ_DIR="$OUT/fastq" OUT="$OUT/qc" bash "$SCRIPTS/01_fastqc_multiqc.sh"
echo "Per bigWig e picchi: lanciare 05_bigwig_cpm.sh e 07_macs3_peak_calling.sh con BAM_DIR=$OUT/bam"
echo "(nomi dei BAM: *_PolII_* = IP, *_Input_* = Input; adattare il riconoscimento IP/Input se necessario)"
