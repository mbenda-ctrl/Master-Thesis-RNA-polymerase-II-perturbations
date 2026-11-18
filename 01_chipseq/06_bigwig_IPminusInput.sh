#!/usr/bin/env bash
# =============================================================================================
# PASSO 0 - bigWig IP - Input (CPM) per tutti i campioni ChIP Pol II.
#
# Per ogni campione: bigwigCompare --operation subtract (IP - Input) a bin di 10 bp.
# I bigWig di partenza sono gia' normalizzati in CPM (bamCoverage --normalizeUsing CPM, --extendReads 200),
# quindi la sottrazione e' fatta su segnali confrontabili. I valori negativi (IP sotto l'input) restano nel
# bigWig; gli script a valle li portano a 0 quando serve.
#   - salta i file gia' presenti (si puo' interrompere e rilanciare)
#   - cerca l'input di ogni campione nel nome dei file (chiave in INPUT_KEY)
#   - se manca un IP o il suo input, lo segnala e passa al campione successivo
# Output: $OUT/<campione>_IPminusInput.bw  (usati dai passi 5, 6 e 7)
#
# Uso:  conda activate deeptools_env && bash 00_bigwig_IPminusInput.sh
# =============================================================================================
set -uo pipefail          # -u: errore su variabili non definite; pipefail: una pipe fallisce se fallisce un comando
                          # (niente -e: un campione con errore non deve fermare gli altri)

# ---------------- impostazioni ----------------
BW="${BW:-results/chipseq/bigwig_cpm}"              # bigWig IP e Input (CPM, frammenti estesi a 200 bp)
OUT="${OUT:-results/chipseq/bigwig_IPminusInput}"          # cartella di uscita
NP=3                                      # processori
BIN=10                                    # bin (bp) del bigWig in uscita

# etichetta del campione -> file IP
declare -A IP=(
  [Amr_r1]="S61454_IP_Amr_1_S33.bw"     [Amr_r2]="S61466_IP_Amr_2_S2.bw"
  [Fast_r1]="S61455_IP_E11_1_S50.bw"      [Fast_r2]="S61467_IP_E22_2_S10.bw"
  [Slow_r1]="S61456_IP_C412_1_S49.bw"   [Slow_r2]="S61468_IP_C422_2_S47.bw"
  [UT_r1]="S61451_IP_WT_1_S77.bw"       [UT_r2]="S61465_IP_Nacl_2_S12.bw"
  [Storm_r1]="S61452_IP_Storm_1_S58.bw" [Storm_r2]="S61464_IP_Storm_2_S51.bw"
  [UTpub_r1]="HEK293T_UT_PolII_R1_CPM.bw"  [UTpub_r2]="HEK293T_UT_PolII_R2_CPM.bw"     # dati pubblici (UT*)
  [KCl_r1]="HEK293T_KCl_PolII_R1_CPM.bw"   [KCl_r2]="HEK293T_KCl_PolII_R2_CPM.bw"
)
# etichetta -> parte del nome con cui riconoscere il bigWig dell'input corrispondente
declare -A INPUT_KEY=(
  [Amr_r1]="Input_Amr_1"    [Amr_r2]="Input_Amr_2"
  [Fast_r1]="Input_E11"     [Fast_r2]="Input_E22"
  [Slow_r1]="Input_C412"    [Slow_r2]="Input_C422"
  [UT_r1]="Input_WT_1"      [UT_r2]="Input_Nacl_2"
  [Storm_r1]="Input_Storm_1" [Storm_r2]="Input_Storm_2"
  [UTpub_r1]="UT_Input_R1"  [UTpub_r2]="UT_Input_R2"
  [KCl_r1]="KCl_Input_R1"   [KCl_r2]="KCl_Input_R2"
)
SAMPLES="Amr_r1 Amr_r2 Fast_r1 Fast_r2 Slow_r1 Slow_r2 UT_r1 UT_r2 Storm_r1 Storm_r2 UTpub_r1 UTpub_r2 KCl_r1 KCl_r2"

# ---------------- controlli ----------------
mkdir -p "$OUT/tmp" 2>/dev/null
# prova di scrittura: se il disco E: non e' montato ci si ferma subito
touch "$OUT/.test_scrittura" 2>/dev/null && rm -f "$OUT/.test_scrittura" || {
  echo "ERRORE: non riesco a scrivere in $OUT. Controlla il percorso di OUT"; exit 1; }
export TMPDIR="$OUT/tmp"                   # temporanei su E:, non dentro Ubuntu
command -v bigwigCompare >/dev/null || { echo "ERRORE: bigwigCompare non trovato (conda activate deeptools_env)"; exit 1; }

LOG="$OUT/log_$(date +%Y%m%d_%H%M).txt"   # log con data e ora
declare -a FATTI SALTATI MANCANTI ERRORI  # liste per il riepilogo finale

# ---------------- calcolo ----------------
for s in $SAMPLES; do
  o="$OUT/${s}_IPminusInput.bw"
  if [[ -s "$o" ]]; then SALTATI+=("$s"); echo "[$s] gia' presente, salto"; continue; fi   # gia' fatto

  ipf="$BW/${IP[$s]}"
  if [[ ! -f "$ipf" ]]; then MANCANTI+=("$s (IP: ${IP[$s]})"); echo "[$s] IP non trovato, salto"; continue; fi

  # input: deve esserci esattamente un bigWig che contiene la chiave nel nome
  mapfile -t hits < <(ls "$BW" | grep -i -- "${INPUT_KEY[$s]}" | grep -i '\.bw$')
  if [[ ${#hits[@]} -ne 1 ]]; then
    MANCANTI+=("$s (input '${INPUT_KEY[$s]}': trovati ${#hits[@]} file ${hits[*]:-})")
    echo "[$s] input non univoco o mancante, salto"; continue
  fi
  inf="$BW/${hits[0]}"

  echo "[$s] $(date +%H:%M)  IP=${IP[$s]}  Input=${hits[0]}" | tee -a "$LOG"
  # scrive prima su .part e rinomina solo se il comando e' finito bene (niente file a meta')
  if bigwigCompare -b1 "$ipf" -b2 "$inf" --operation subtract --binSize "$BIN" -p "$NP" -o "$o.part" >>"$LOG" 2>&1; then
    mv "$o.part" "$o"; FATTI+=("$s")
  else
    rm -f "$o.part"; ERRORI+=("$s"); echo "[$s] ERRORE, vedi $LOG"
  fi
done

# ---------------- riepilogo ----------------
echo
echo "Creati:   ${FATTI[*]:-nessuno}"
echo "Saltati (gia' presenti): ${SALTATI[*]:-nessuno}"
echo "Errori:   ${ERRORI[*]:-nessuno}"
if [[ ${#MANCANTI[@]} -gt 0 ]]; then
  echo "File mancanti o da sistemare in INPUT_KEY:"; printf '  - %s\n' "${MANCANTI[@]}"
fi
echo "Output in: $OUT"
