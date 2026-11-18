#!/usr/bin/env bash
# =============================================================================================
# PASSO 12b - Conte Pol II per gene (IP e Input) e dimensioni di libreria, per integrazione e heatmap.
#   1. geni: protein-coding senza sovrapposizioni (passo 12a) intersecati PER NOME con il set consensus
#      del passo 8 (geni con picchi); regione = gene intero. Nomi doppi: si tiene la prima copia.
#      -> geni_intersezione_MANE.bed  (BED dei geni su cui si assegnano le regioni DiffBind)
#   2. featureCounts su 10 IP e 10 Input: paired-end = si contano i frammenti (-p --countReadPairs)
#   3. read mappate di ogni BAM (samtools idxstats) -> librerie.tsv
# I campioni (etichetta, BAM IP, BAM Input) sono in config/chipseq_samples.tsv.
#
# Uso:  NONOVERLAP=results/chipseq/genes/genes_pc_nonoverlap.bed \
#       CONSENSUS=results/chipseq/consensus/genes_pc_union_MANE_n2_min1kb.bed \
#       BAM_DIR=results/chipseq/bam OUT=results/chipseq/gene_counts bash 12b_count_chip_genes.sh
# Output: $OUT/geni_intersezione_MANE.bed, $OUT/conteggi_gene.tab (gene, chr, start, end, <IP>, <IP>_input),
#         $OUT/librerie.tsv (campione, bam, read_mappate)
# =============================================================================================
set -euo pipefail

NONOVERLAP="${NONOVERLAP:?indica il BED del passo 12a}"
CONSENSUS="${CONSENSUS:?indica il BED consensus del passo 8 (nome del gene in colonna 4)}"
BAM_DIR="${BAM_DIR:-results/chipseq/bam}"
SAMPLES="${SAMPLES:-config/chipseq_samples.tsv}"
OUT="${OUT:-results/chipseq/gene_counts}"
THREADS="${THREADS:-4}"

command -v featureCounts >/dev/null || { echo "ERRORE: featureCounts non trovato (conda install -c bioconda subread)"; exit 1; }
mkdir -p "$OUT"

# ---------------- 1. BAM e librerie: prima i 10 IP, poi i 10 Input ----------------
B=(); L=()
echo -e "campione\tbam\tread_mappate" > "$OUT/librerie.tsv"
for col in ip_bam input_bam; do
  while IFS=$'\t' read -r label cond rep ip inp; do
    [[ "$label" == label ]] && continue
    if [[ $col == ip_bam ]]; then bam="$BAM_DIR/$ip"; lab="$label"; else bam="$BAM_DIR/$inp"; lab="${label}_input"; fi
    [[ -s "$bam" ]] || { echo "ERRORE: manca $bam"; exit 1; }
    B+=("$bam"); L+=("$lab")
    n=$(samtools idxstats "$bam" | awk '{s += $3} END {print s}')
    echo -e "$lab\t$(basename "$bam")\t$n" >> "$OUT/librerie.tsv"
  done < "$SAMPLES"
done
cat "$OUT/librerie.tsv"

# ---------------- 2. geni: senza sovrapposizioni ∩ consensus (per nome) ----------------
PRIMO=$(samtools idxstats "${B[0]}" | awk 'NR == 1 {print $1}')
[[ "$PRIMO" == chr* ]] && PRE="chr" || PRE=""
awk -v OFS='\t' -v P="$PRE" '
  NR == FNR { nostri[$4] = 1; next }
  ($4 in nostri) { c = $1; sub(/^chr/, "", c); if (seen[$4]++) next; print P c, $2, $3, $4 }
' "$CONSENSUS" "$NONOVERLAP" | sort -k1,1 -k2,2n > "$OUT/geni_intersezione_MANE.bed"
echo "Geni senza sovrapposizioni: $(wc -l < "$NONOVERLAP") | consensus: $(wc -l < "$CONSENSUS") | in comune: $(wc -l < "$OUT/geni_intersezione_MANE.bed")"
awk -v OFS='\t' 'BEGIN {print "GeneID", "Chr", "Start", "End", "Strand"} {print $4, $1, $2 + 1, $3, "+"}' \
  "$OUT/geni_intersezione_MANE.bed" > "$OUT/geni_intersezione_MANE.saf"

# ---------------- 3. featureCounts (frammenti paired-end, senza filamento) ----------------
featureCounts -F SAF -a "$OUT/geni_intersezione_MANE.saf" -o "$OUT/featureCounts_gene.txt" \
  -T "$THREADS" -s 0 -p --countReadPairs "${B[@]}"

# ---------------- 4. tabella finale con le etichette dei campioni ----------------
awk -v OFS='\t' -v LAB="${L[*]}" '
  BEGIN { n = split(LAB, lab, " "); h = "gene" OFS "chr" OFS "start" OFS "end"; for (i = 1; i <= n; i++) h = h OFS lab[i]; print h }
  /^#/ || $1 == "Geneid" { next }
  { o = $1 OFS $2 OFS ($3 - 1) OFS $4; for (i = 7; i <= NF; i++) o = o OFS $i; print o }
' "$OUT/featureCounts_gene.txt" > "$OUT/conteggi_gene.tab"
echo "Fatto: $OUT/conteggi_gene.tab, $OUT/librerie.tsv, $OUT/geni_intersezione_MANE.bed"
