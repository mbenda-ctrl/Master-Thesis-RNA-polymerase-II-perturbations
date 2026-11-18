#!/usr/bin/env bash
# =============================================================================================
# METAGENE 1 - Plot definitivi dei metagene Pol II.
# 3 contrasti x 3 regioni (TSS, TES, gene completo) x 2 scale (CPM, logCPM) = 18 plot.
#   - bigWig IP - Input (01_chipseq/06_bigwig_IPminusInput.sh), fianchi +-3 kb
#   - geni MANE >= 5 kb; bin 10 bp per TSS e TES (reference-point), 50 bp per il gene completo
#     (scale-regions: corpo di ogni gene scalato a 35 kb)
#   - per ogni regione la stessa scala Y in tutti i contrasti
# Crea anche genes_MANE_min5kb.bed, usato dai passi 5 e 8.
#
# Uso:  conda activate deeptools_env
#       bash 01_run_metagene_profiles.sh 2>&1 | tee plot_finali.log
#       (prova veloce su un solo campione: CONTRASTS=Amr_test bash 01_run_metagene_profiles.sh)
# =============================================================================================
set -euo pipefail

SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"        # cartella di questo script

# ---------------- impostazioni ----------------
BWDIR="${BWDIR:-results/chipseq/bigwig_IPminusInput}"
BED="${BED:-results/chipseq/consensus/genes_pc_union_MANE_n2_min1kb.bed}"   # geni MANE (tutti)
OUT="${OUT:-results/chipseq/plot_finali}"
PY="$SCRIPTS/01a_plot_metagene_profile.py"
NP=3                 # processori
BIN=10               # bin per TSS e TES
FLANK=3000           # +-3 kb per tutti i plot
BODY=35000           # unita' codificante scalata a 35 kb
BODY_BIN=50          # bin per l'unita' codificante
TITLE_BODY="Transcript scaled"

# asse Y (minimo sempre 0): TSS e unita' codificante 0-10, TES 0-1; versione log 0-3.5 / 0-1
TSS_YMAX=10;     TSS_YMAX_LOG=3.5
BODY_YMAX=10;    BODY_YMAX_LOG=3.5
TES_YMAX=1;      TES_YMAX_LOG=1

# versione "log": data = log2(CPM+1) su asse lineare (etichetta logCPM)
#                 axis = CPM su asse logaritmico (picco con lo stesso valore del lineare)
LOG_STYLE=data
MAKE_LOG=yes         # yes = genera anche le versioni log

# contrasti: nome -> etichette dei campioni (= nomi dei bigWig <etichetta>_IPminusInput.bw)
declare -A SAMPLES=(
  [AmrFastSlow]="Amr_r1 Amr_r2 Fast_r1 Fast_r2 Slow_r1 Slow_r2"
  [UTStorm]="UT_r1 UT_r2 Storm_r1 Storm_r2"
  # [KCl_UTpub]="UTpub_r1 UTpub_r2 KCl_r1 KCl_r2"   # dataset pubblico, da integrare
  [Amr_test]="Amr_r1 Amr_r2"
)
# di default i 3 contrasti; si puo' cambiare dall'esterno con la variabile CONTRASTS
CONTRASTS="${CONTRASTS:-AmrFastSlow UTStorm}"   # + KCl_UTpub (dataset pubblico)

# ---------------- preparazione ----------------
mkdir -p "$OUT/matrici" "$OUT/plot" "$OUT/tmp"
export TMPDIR="$OUT/tmp"
[[ -f "$PY" ]] || { echo "ERRORE: manca $PY"; exit 1; }
for c in $CONTRASTS; do                                          # tutti i bigWig devono esistere
  for s in ${SAMPLES[$c]}; do
    [[ -s "$BWDIR/${s}_IPminusInput.bw" ]] || { echo "ERRORE: manca $BWDIR/${s}_IPminusInput.bw"; exit 1; }
  done
done
GENES="$OUT/genes_MANE_min5kb.bed"
awk '$3-$2 >= 5000' "$BED" > "$GENES"                            # solo geni lunghi almeno 5 kb


# ---------------- 1. matrici ----------------
for c in $CONTRASTS; do
  S=(); L=()
  # etichetta UTpub mostrata come UT* nei plot
  for s in ${SAMPLES[$c]}; do S+=("$BWDIR/${s}_IPminusInput.bw"); L+=("${s/UTpub/UT*}"); done
  COMMON=(-S "${S[@]}" --samplesLabel "${L[@]}" --missingDataAsZero -p "$NP")   # opzioni comuni
  for R in TSS TES BODY; do
    M="$OUT/matrici/matrix_${R}_${c}.gz"
    [[ -s "$M" ]] && continue                                    # matrice gia' fatta
    echo "=== $c: $R ==="
    if [[ $R == BODY ]]; then
      # gene completo: TSS -> TES scalato a BODY bp, piu' FLANK a monte e a valle
      computeMatrix scale-regions "${COMMON[@]}" -R "$GENES" --binSize "$BODY_BIN" \
        -b "$FLANK" -a "$FLANK" -m "$BODY" -o "$M"
    else
      # TSS o TES: finestra di +-FLANK attorno al punto di riferimento
      computeMatrix reference-point --referencePoint "$R" "${COMMON[@]}" -R "$GENES" --binSize "$BIN" \
        -b "$FLANK" -a "$FLANK" -o "$M"
    fi
  done
done

# ---------------- 2. plot (scala Y comune per regione) ----------------
for R in TSS TES BODY; do
  ALL=(); for c in $CONTRASTS; do ALL+=("$OUT/matrici/matrix_${R}_${c}.gz"); done   # matrici della regione
  case $R in                                                     # titolo e massimo Y della regione
    TSS)  T="TSS";         Y=$TSS_YMAX;  YL=$TSS_YMAX_LOG ;;
    TES)  T="TES";         Y=$TES_YMAX;  YL=$TES_YMAX_LOG ;;
    BODY) T="$TITLE_BODY"; Y=$BODY_YMAX; YL=$BODY_YMAX_LOG ;;
  esac
  for c in $CONTRASTS; do
    M="$OUT/matrici/matrix_${R}_${c}.gz"
    python "$PY" "$M" --out "$OUT/plot/${c}_${R}_CPM"    --title "$T" --ymax "$Y"  --scale-from "${ALL[@]}"
    [[ $MAKE_LOG == yes ]] || continue
    if [[ $LOG_STYLE == axis ]]; then
      python "$PY" "$M" --out "$OUT/plot/${c}_${R}_logscale" --title "$T" --ymax "$Y" --scale-from "${ALL[@]}" --log-axis
    else
      python "$PY" "$M" --out "$OUT/plot/${c}_${R}_logCPM" --title "$T" --ymax "$YL" --scale-from "${ALL[@]}" --log
    fi
  done
done

echo "Fatto: $(ls "$OUT"/plot/*.png | wc -l) plot in $OUT/plot"
