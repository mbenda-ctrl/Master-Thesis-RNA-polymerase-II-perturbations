#!/usr/bin/env bash
# PASSO 4 - QC dei BAM: flagstat + stima delle read multimappate per tutti i BAM delle cartelle indicate.
# Uso: BAM_DIRS="results/chipseq/bam" OUT=results/chipseq/bam_qc bash 04_bam_qc.sh
# Output: un .flagstat per BAM e un riepilogo unico bam_stats_summary.tsv
# Richiede samtools >= 1.13 (usa l'ambiente qc, non deeptools_env).
set -euo pipefail

read -r -a DIRS <<< "${BAM_DIRS:-results/chipseq/bam}"   # una o piu' cartelle di BAM, separate da spazio
OUT="${OUT:-results/chipseq/bam_qc}"
NP="${THREADS:-4}"
mkdir -p "$OUT"

samtools --version | head -1

SUMMARY="$OUT/bam_stats_summary.tsv"
printf "sample\ttotal\tmapped\tmapped_pct\tproperly_paired_pct\tduplicates\tsecondary\tprimary_mapped\tmapq_lt30\tmapq_lt30_pct\tmapq0\tmapq0_pct\tXS_tag\tXS_tag_pct\n" > "$SUMMARY"

for dir in "${DIRS[@]}"; do
  for bam in "$dir"/*.bam; do
    name=$(basename "$bam" .bam)
    echo "=== $name ==="

    # 1. flagstat completo (leggibile) + versione tsv per estrarre i valori
    samtools flagstat -@ "$NP" "$bam" > "$OUT/${name}.flagstat"
    samtools flagstat -@ "$NP" -O tsv "$bam" > "$OUT/${name}.flagstat.tsv"
    get() { awk -F'\t' -v k="$1" '$3==k {print $1; exit}' "$OUT/${name}.flagstat.tsv"; }

    total=$(get "total (QC-passed reads + QC-failed reads)")
    mapped=$(get "mapped")
    mapped_pct=$(get "mapped %")
    pp_pct=$(get "properly paired %")
    dups=$(get "duplicates")
    sec=$(get "secondary")

    # 2. multimapping, sulle sole read primarie mappate:
    #    MAPQ < 30  -> read con allineamento ambiguo (soglia usata per i file q30)
    #    MAPQ == 0  -> read con piu' posizioni equivalenti (multimapper "veri")
    #    tag XS     -> esiste un allineamento alternativo (bowtie2/bwa)
    read -r prim low z xs < <(
      samtools view -@ "$NP" -F 0x904 "$bam" | awk '
        { n++; if ($5 < 30) low++; if ($5 == 0) z++
          for (i = 12; i <= NF; i++) if ($i ~ /^XS:i:/) { xs++; break } }
        END { print n+0, low+0, z+0, xs+0 }')

    pct() { awk -v a="$1" -v b="$2" 'BEGIN{ if (b > 0) printf "%.2f", 100*a/b; else print "NA" }'; }

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$name" "$total" "$mapped" "$mapped_pct" "$pp_pct" "$dups" "$sec" \
      "$prim" "$low" "$(pct "$low" "$prim")" "$z" "$(pct "$z" "$prim")" "$xs" "$(pct "$xs" "$prim")" \
      >> "$SUMMARY"
  done
done

echo
column -t -s $'\t' "$SUMMARY"
echo "Riepilogo salvato in: $SUMMARY"

# report MultiQC dei flagstat
multiqc -f -n multiqc_bam_flagstat -o "$OUT" "$OUT"
