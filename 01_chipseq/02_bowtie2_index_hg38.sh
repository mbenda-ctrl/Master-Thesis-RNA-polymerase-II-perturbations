#!/usr/bin/env bash
# =============================================================================================
# PASSO 2 - Indice Bowtie2 del genoma umano GRCh38 (hg38).
#
# Di default scarica il "no_alt analysis set" di GRCh38 (NCBI, nomi dei cromosomi in stile UCSC:
# chr1, chr2, ...), la versione consigliata per l'allineamento perche' esclude i contig alternativi
# che renderebbero multimappate molte read. Si puo' indicare un FASTA gia' presente con GENOME_FA.
#
# Uso:  OUTDIR=reference/hg38 THREADS=8 bash 02_bowtie2_index_hg38.sh
#       GENOME_FA=/percorso/hg38.fa OUTDIR=reference/hg38 bash 02_bowtie2_index_hg38.sh
# =============================================================================================
set -euo pipefail

OUTDIR="${OUTDIR:-reference/hg38}"
THREADS="${THREADS:-8}"
URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
PREFIX="$OUTDIR/GRCh38_no_alt"

mkdir -p "$OUTDIR"

if [[ -z "${GENOME_FA:-}" ]]; then
  GENOME_FA="$OUTDIR/GRCh38_no_alt_analysis_set.fna"
  if [[ ! -s "$GENOME_FA" ]]; then
    echo ">>> Scarico GRCh38 no_alt analysis set"
    wget -c -O "$GENOME_FA.gz" "$URL"
    gunzip -f "$GENOME_FA.gz"
  fi
fi

echo ">>> bowtie2-build ($GENOME_FA)"
bowtie2-build --threads "$THREADS" "$GENOME_FA" "$PREFIX" > "$OUTDIR/bowtie2-build.log" 2>&1
samtools faidx "$GENOME_FA"
cut -f1,2 "$GENOME_FA.fai" > "$OUTDIR/chrom.sizes"
echo "Indice pronto: $PREFIX  (usa BT2_INDEX=$PREFIX nel passo 3)"
