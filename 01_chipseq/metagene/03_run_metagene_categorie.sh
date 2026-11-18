#!/usr/bin/env bash
# =============================================================================================
# METAGENE 3 - Metagene Pol II per categorie RNA-seq (tutti i confronti).
#   UP        padj < 0.05 e log2FC >  0.58
#   DOWN      padj < 0.05 e log2FC < -0.58
#   NO CHANGE -0.3 < log2FC < 0.3 e baseMean >= 50
# Per ogni confronto:
#   1) BED delle 3 categorie (02_make_bed_categorie.py: intersezione con i geni MANE >= 5 kb)
#   2) matrici computeMatrix TSS / TES / gene completo sui bigWig IP - Input (01_chipseq/06_bigwig_IPminusInput.sh), stessi parametri dei plot finali
#   3) plot: combinato (6 profili, 05a), profilo differenziale per categoria (05b), boxplot per regione (05c)
# Il controllo WT dell'RNA-seq e' lo stesso UT del ChIP: nei nomi e nei titoli e' sempre UT.
# Le matrici gia' presenti vengono saltate (per rifarle, cancellale).
#
# Uso:  conda activate deeptools_env
#       bash 03_run_metagene_categorie.sh 2>&1 | tee metagene_categorie.log
# =============================================================================================
set -euo pipefail          # si ferma al primo errore, su variabili non definite e su errori nelle pipe

# cartella in cui si trova questo script: gli script Python dei passi 3 e 5a-c stanno accanto
SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------- soglie ----------------
LFC=0.58          # |log2FC| per UP / DOWN
PADJ=0.05         # padj per UP / DOWN
NC=0.3            # NO CHANGE: |log2FC| < NC
MINBM=50          # NO CHANGE: baseMean minimo

# ---------------- percorsi ----------------
BWDIR="${BWDIR:-results/chipseq/bigwig_IPminusInput}"                    # bigWig IP - Input: <cond>_<r1|r2>_IPminusInput.bw (01_chipseq/06_bigwig_IPminusInput.sh)
GENES="${GENES:-results/chipseq/plot_finali/genes_MANE_min5kb.bed}"      # geni dei metagene (MANE >= 5 kb, creato da 01_run_metagene_profiles.sh)
RNA_DIR="${RNA_DIR:-results/rnaseq/RNAseq}"    # <confronto>_all_genes.csv (02_rnaseq/01_deseq2_rnaseq.R)
OUT="${OUT:-results/chipseq/metagene_categorie/lfc058_nc03}"           # cartella di uscita
PY_BED="$SCRIPTS/02_make_bed_categorie.py"
PY_COMB="$SCRIPTS/03a_plot_categorie_combinato.py"
PY_DIFF="$SCRIPTS/03b_plot_diff_categorie.py"
PY_BOX="$SCRIPTS/03c_plot_box_categorie.py"

# ---------------- parametri computeMatrix (come i plot finali) ----------------
NP=3              # processori
BIN=10            # bin per TSS e TES (bp)
FLANK=3000        # finestra a monte e a valle (bp)
BODY=35000        # lunghezza a cui viene scalato il corpo del gene (bp)
BODY_BIN=50       # bin per il gene completo (bp)

# ---------------- confronti ----------------
# nome ChIP -> "controllo_ChIP  trattato_ChIP  confronto_RNA (nome dei BED, gia' con UT)"
declare -A COMP
COMP[Fast_vs_Amr]="Amr Fast Fast_vs_Amr"
COMP[Slow_vs_Amr]="Amr Slow Slow_vs_Amr"
COMP[Storm_vs_UT]="UT Storm Storm_vs_UT"
COMP[KCl_vs_UTpub]="UTpub KCl NaCl_vs_UT"            # ChIP pubblico KCl / UT*, categorie RNA-seq NaCl vs UT
ORDER="Fast_vs_Amr Slow_vs_Amr Storm_vs_UT KCl_vs_UTpub"

# assi y fissi: TSS e gene completo stesso massimo (combinato); differenziale e boxplot +-3
declare -A YMAX=([TSS]=14 [BODY]=14 [TES]=1)
declare -A RLAB=([TSS]="TSS" [BODY]="Transcript scaled" [TES]="TES")
YLIM_DIFF=3

# ---------------- controlli ----------------
for f in "$GENES" "$PY_BED" "$PY_COMB" "$PY_DIFF" "$PY_BOX"; do
  [[ -s "$f" ]] || { echo "ERRORE: manca $f"; exit 1; }   # -s: il file esiste e non e' vuoto
done
mkdir -p "$OUT/bed_categorie" "$OUT/matrici" "$OUT/plot" "$OUT/tmp"
export TMPDIR="$OUT/tmp"                             # file temporanei di deepTools su E:

# ---------------- 1. BED delle categorie (tutti i confronti in una volta) ----------------
python "$PY_BED" --bed "$GENES" --rnaseq-dir "$RNA_DIR" --outdir "$OUT/bed_categorie" \
  --lfc "$LFC" --padj "$PADJ" --nc "$NC" --min-bm "$MINBM"

for C in $ORDER; do
  read -r CTRL TREAT RNA <<< "${COMP[$C]}"           # spacchetta "controllo trattato confronto_RNA"
  RNALAB=${RNA//_vs_/ vs }                            # Fast_vs_Amr -> "Fast vs Amr" per i titoli
  TITLE="geni RNA-seq $RNALAB: UP/DOWN |log2FC| > $LFC, NC |log2FC| < $NC"
  # i 3 BED nell'ordine UP, NOCHANGE, DOWN (gli script di plot si aspettano quest'ordine)
  R=("$OUT/bed_categorie/${RNA}_UP.bed" "$OUT/bed_categorie/${RNA}_NOCHANGE.bed" "$OUT/bed_categorie/${RNA}_DOWN.bed")
  echo "=== $C: UP $(wc -l < "${R[0]}") | NO CHANGE $(wc -l < "${R[1]}") | DOWN $(wc -l < "${R[2]}") ==="

  # bigWig: controllo r1, r2, trattato r1, r2 (le etichette <cond>_r1 servono agli script di plot)
  S=(); L=()
  for cond in "$CTRL" "$TREAT"; do for r in r1 r2; do
    f="$BWDIR/${cond}_${r}_IPminusInput.bw"; [[ -s "$f" ]] || { echo "ERRORE: manca $f"; exit 1; }
    S+=("$f"); L+=("${cond}_${r}")
  done; done

  # ---------------- 2. matrici ----------------
  # TSS: reference-point, +-3 kb attorno al TSS, bin 10 bp; --missingDataAsZero: bin senza dati = 0
  M="$OUT/matrici/matrix_TSS_${C}.gz"
  [[ -s "$M" ]] || computeMatrix reference-point --referencePoint TSS -S "${S[@]}" --samplesLabel "${L[@]}" \
    -R "${R[@]}" -b "$FLANK" -a "$FLANK" --binSize "$BIN" --missingDataAsZero -p "$NP" -o "$M"
  # TES: reference-point, +-3 kb attorno al TES
  M="$OUT/matrici/matrix_TES_${C}.gz"
  [[ -s "$M" ]] || computeMatrix reference-point --referencePoint TES -S "${S[@]}" --samplesLabel "${L[@]}" \
    -R "${R[@]}" -b "$FLANK" -a "$FLANK" --binSize "$BIN" --missingDataAsZero -p "$NP" -o "$M"
  # gene completo: scale-regions, ogni gene (TSS -> TES) scalato a 35 kb + 3 kb di fianchi, bin 50 bp
  M="$OUT/matrici/matrix_BODY_${C}.gz"
  [[ -s "$M" ]] || computeMatrix scale-regions -S "${S[@]}" --samplesLabel "${L[@]}" \
    -R "${R[@]}" -b "$FLANK" -a "$FLANK" -m "$BODY" --binSize "$BODY_BIN" --missingDataAsZero -p "$NP" -o "$M"

  # ---------------- 3. plot ----------------
  for Rg in TSS BODY TES; do
    M="$OUT/matrici/matrix_${Rg}_${C}.gz"
    # 6 profili: controllo e trattato x UP / NO CHANGE / DOWN
    python "$PY_COMB" "$M" --ctrl "$CTRL" --treat "$TREAT" --ymax "${YMAX[$Rg]}" \
      --out "$OUT/plot/${C}_${Rg}_combinato" --title "${RLAB[$Rg]} - $TITLE"
    # log2FC trattato / controllo per categoria, con banda IC 95%
    python "$PY_DIFF" "$M" --ctrl "$CTRL" --treat "$TREAT" --ylim "$YLIM_DIFF" \
      --out "$OUT/plot/${C}_${Rg}_differenziale" --title "${RLAB[$Rg]} - $TITLE"
  done
  # boxplot del log2FC per gene in promotore / gene body / terminazione (usa le 3 matrici)
  python "$PY_BOX" --tss "$OUT/matrici/matrix_TSS_${C}.gz" --body "$OUT/matrici/matrix_BODY_${C}.gz" \
    --tes "$OUT/matrici/matrix_TES_${C}.gz" --ctrl "$CTRL" --treat "$TREAT" --ylim "$YLIM_DIFF" \
    --out "$OUT/plot/${C}_boxplot_regioni" --title "Pol II per regione - $TITLE"
done
echo "Fatto. Output in: $OUT/plot"
