#!/usr/bin/env bash
# =============================================================================================
# METAGENE 4 - Metagene differenziali log2FC (trattato / controllo) con 04a_differential_metagene.py.
# Scale-regions come i plot finali: -3 kb | corpo scalato a 35 kb | +3 kb, bin 100 bp.
# Un set di geni differenziali per contrasto (Storm: tutti i geni). Per ogni contrasto una matrice con
# 4 campioni (2 trattati + 2 controlli), poi profilo + heatmap + barra delle lunghezze + scatter.
# NB: qui si usano i bigWig IP (non IP - Input): niente valori negativi -> log2FC sempre definito.
#
# Uso:  conda activate deeptools_env
#       bash 04_run_differential_metagene.sh 2>&1 | tee differential_metagene.log
# =============================================================================================
set -euo pipefail

SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"        # cartella di questo script

# ---------------- impostazioni ----------------
BW="${BW:-results/chipseq/bigwig_cpm}"          # bigWig IP (CPM)
GENES_ALL="${GENES_ALL:-results/chipseq/plot_finali/genes_MANE_min5kb.bed}"   # 12097 geni >= 5 kb (01_run_metagene_profiles.sh): usato per Storm
SETDIR="${SETDIR:-results/chipseq/diffbind_gene_sets}"              # cartella con i BED dei geni differenziali (genes_diff_*.bed)
MINLEN=5000                 # come i plot finali: geni < 5 kb esclusi (0 = tienili tutti)
OUT="${OUT:-results/chipseq/differential_bin100}"
PY="$SCRIPTS/04a_differential_metagene.py"
NP=3
BIN=100; FLANK=3000; BODY=35000    # bin, fianchi e corpo scalato (bp)
PSEUDO=0.1                  # pseudoconto nel log2((A + ps) / (B + ps))
ZMAX=1.5                    # scala colore heatmap fissa (-1.5 .. 1.5)
YMIN=-3; YMAX=3             # asse y del profilo log2FC, fisso per tutti i contrasti

# campione -> bigWig IP
declare -A F=(
  [Amr_r1]=S61454_IP_Amr_1_S33.bw     [Amr_r2]=S61466_IP_Amr_2_S2.bw
  [Fast_r1]=S61455_IP_E11_1_S50.bw      [Fast_r2]=S61467_IP_E22_2_S10.bw
  [Slow_r1]=S61456_IP_C412_1_S49.bw   [Slow_r2]=S61468_IP_C422_2_S47.bw
  [UT_r1]=S61451_IP_WT_1_S77.bw       [UT_r2]=S61465_IP_Nacl_2_S12.bw
  [Storm_r1]=S61452_IP_Storm_1_S58.bw [Storm_r2]=S61464_IP_Storm_2_S51.bw
  # dataset pubblico, da integrare:
  # [UTpub_r1]=HEK293T_UT_PolII_R1_CPM.bw  [UTpub_r2]=HEK293T_UT_PolII_R2_CPM.bw
  # [KCl_r1]=HEK293T_KCl_PolII_R1_CPM.bw   [KCl_r2]=HEK293T_KCl_PolII_R2_CPM.bw
)
# contrasto -> set di geni (BED in $SETDIR, creati da 01_chipseq/13_diffbind_gene_sets.py; piu' BED separati da spazio = uniti; ALL = tutti i geni >= 5 kb)
declare -A SET=(
  [Fast_vs_Amr_UP]="genes_diff_FastvsAmr_UP.bed"
  [Fast_vs_Amr_DOWN]="genes_diff_FastvsAmr_DOWN.bed"
  [Slow_vs_Amr]="genes_diff_SlowvsAmr.bed"
  [Storm_vs_UT]="ALL"
  # [KCl_vs_UTpub]="genes_diff_KCLvsUT.bed"   # dataset pubblico
)
# contrasto -> "trattato controllo"
declare -A C=([Fast_vs_Amr_UP]="Fast Amr" [Fast_vs_Amr_DOWN]="Fast Amr" [Slow_vs_Amr]="Slow Amr" [Storm_vs_UT]="Storm UT" [KCl_vs_UTpub]="KCl UTpub")
ORDER="Fast_vs_Amr_UP Fast_vs_Amr_DOWN Slow_vs_Amr Storm_vs_UT"   # + KCl_vs_UTpub (dataset pubblico)

# ---------------- controlli ----------------
for f in "$GENES_ALL" "$PY"; do [[ -s "$f" ]] || { echo "ERRORE: manca $f"; exit 1; }; done
for s in "${!F[@]}"; do [[ -s "$BW/${F[$s]}" ]] || { echo "ERRORE: manca $BW/${F[$s]}"; exit 1; }; done
for c in $ORDER; do
  for b in ${SET[$c]}; do [[ "$b" == ALL || -s "$SETDIR/$b" ]] || { echo "ERRORE: manca $SETDIR/$b"; exit 1; }; done
done
mkdir -p "$OUT/set_geni" "$OUT/matrici" "$OUT/plot" "$OUT/tmp"
export TMPDIR="$OUT/tmp"

for c in $ORDER; do
  read -r T K <<< "${C[$c]}"                                     # T = trattato, K = controllo
  # set di geni del contrasto: unione dei BED, solo geni >= MINLEN, senza doppioni
  G="$OUT/set_geni/geni_${c}.bed"
  if [[ "${SET[$c]}" == ALL ]]; then cp "$GENES_ALL" "$G"
  else
    for b in ${SET[$c]}; do cat "$SETDIR/$b"; done | awk -v m="$MINLEN" 'NF>=6 && $3-$2>=m' | sort -u > "$G"
  fi
  echo "$c: $(wc -l < "$G") geni (${SET[$c]})"
  M="$OUT/matrici/matrix_BODY_bin${BIN}_${c}_geniDiff.gz"
  if [[ ! -s "$M" ]]; then
    echo "=== computeMatrix $c ==="
    # --sortRegions keep: le righe della matrice restano nell'ordine del BED (serve per la barra delle lunghezze)
    computeMatrix scale-regions \
      -S "$BW/${F[${T}_r1]}" "$BW/${F[${T}_r2]}" "$BW/${F[${K}_r1]}" "$BW/${F[${K}_r2]}" \
      --samplesLabel "${T}_r1" "${T}_r2" "${K}_r1" "${K}_r2" \
      -R "$G" -b "$FLANK" -a "$FLANK" -m "$BODY" --binSize "$BIN" \
      --sortRegions keep --missingDataAsZero -p "$NP" -o "$M"
  fi
  # BED con le stesse righe della matrice, nello stesso ordine (computeMatrix puo' scartare regioni)
  BEDM="$OUT/matrici/regioni_${c}.bed"
  zcat "$M" | tail -n +2 | cut -f1-6 > "$BEDM"

  echo "=== differential metagene $c ==="
  python "$PY" "$M" --groupA "${T}_r1" "${T}_r2" --groupB "${K}_r1" "${K}_r2" \
    --nameA "${T/UTpub/UT*}" --nameB "${K/UTpub/UT*}" --pseudocount "$PSEUDO" \
    --zmax "$ZMAX" --ymin "$YMIN" --ymax "$YMAX" --bed "$BEDM" \
    --out "$OUT/plot/diff_${c}_geniDiff_bin${BIN}"
done
echo "Fatto. Output in: $OUT"
