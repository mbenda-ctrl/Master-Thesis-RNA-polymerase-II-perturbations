#!/usr/bin/env python3
"""
METAGENE 2 - Intersezione dei geni selezionati dall'RNA-seq con i geni usati per il ChIP (BED dei metagene).

Categorie RNA-seq (soglie definitive):
  UP        padj < 0.05  e  log2FC >  0.58      (fold change > 1.5)
  DOWN      padj < 0.05  e  log2FC < -0.58      (fold change < 1/1.5)
  NOCHANGE  -0.3 < log2FC < 0.3  e  baseMean >= 50
            (baseMean >= 50: i geni poco espressi hanno log2FC vicino a 0 solo perche' DESeq2 non ha
             abbastanza conte per stimare una variazione, non perche' siano davvero invariati)
  tutti gli altri geni non entrano in nessun BED.

log2FC (MLE), padj e baseMean sono quelli di DESeq2 sull'RNA-seq, letti dalle tabelle
<confronto>_all_genes.csv di 02_rnaseq/01_deseq2_rnaseq.R (colonna symbol per il simbolo del gene).
Per NaCl vs UT la soglia di padj e' 0.10 (--padj-nacl), come nelle liste DE.
Il collegamento con il BED dei metagene (geni MANE >= 5 kb, simbolo nella 4a colonna) avviene per simbolo.
Se un simbolo compare piu' volte si tiene il gene con baseMean piu' alto.

Nomi: nei file di input il controllo dell'RNA-seq si chiama WT; nei BED in uscita diventa UT
(Storm_vs_WT -> Storm_vs_UT, NaCl_vs_WT -> NaCl_vs_UT), come nel ChIP.

Uso:
  python 02_make_bed_categorie.py --bed genes_MANE_min5kb.bed --rnaseq-dir RNAseq --outdir bed_categorie
Output: <confronto>_<UP|NOCHANGE|DOWN>.bed  +  intersezione_summary.tsv
"""
import argparse, os
import pandas as pd

# ------------------------------------------------------------------------------------------
# parametri (i default sono le soglie decise; si possono cambiare da riga di comando)
# ------------------------------------------------------------------------------------------
ap = argparse.ArgumentParser()
ap.add_argument("--bed", required=True, help="BED dei geni dei metagene (simbolo nella 4a colonna)")
ap.add_argument("--rnaseq-dir", required=True, help="cartella con <confronto>_all_genes.csv di 02_rnaseq/01_deseq2_rnaseq.R")
ap.add_argument("--outdir", required=True)
ap.add_argument("--lfc", type=float, default=0.58, help="soglia |log2FC| per UP e DOWN")
ap.add_argument("--padj", type=float, default=0.05, help="soglia padj per UP e DOWN")
ap.add_argument("--padj-nacl", type=float, default=0.10, help="soglia padj per NaCl vs UT")
ap.add_argument("--nc", type=float, default=0.3, help="NOCHANGE: |log2FC| < nc")
ap.add_argument("--min-bm", type=float, default=50, help="NOCHANGE: baseMean minimo")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)                         # crea la cartella di uscita

# confronto nei file di input -> nome usato in uscita (WT -> UT)
COMPS = {"Fast_vs_Amr": "Fast_vs_Amr", "Slow_vs_Amr": "Slow_vs_Amr",
         "Storm_vs_UT": "Storm_vs_UT", "NaCl_vs_UT": "NaCl_vs_UT"}
CATS = ["UP", "NOCHANGE", "DOWN"]                            # ordine fisso, lo stesso dei gruppi nelle matrici

# BED dei metagene: colonne chr, start, end, simbolo, score, strand (nessuna intestazione)
bed = pd.read_csv(a.bed, sep="\t", header=None)
print(f"BED metagene: {len(bed)} geni")
print(f"Soglie: UP/DOWN padj < {a.padj:g} e |log2FC| > {a.lfc:g};  NOCHANGE |log2FC| < {a.nc:g} e baseMean >= {a.min_bm:g}\n")

rows = []                                                    # righe del riepilogo
for c_in, c_out in COMPS.items():
    # 1. tabella DESeq2 RNA-seq del confronto
    I = pd.read_csv(os.path.join(a.rnaseq_dir, f"{c_in}_all_genes.csv"))
    I = I.rename(columns={"symbol": "gene_name", "log2FC_MLE": "RNAseq_log2FoldChange",
                          "padj": "RNAseq_padj", "baseMean": "RNAseq_baseMean"})
    I = I.dropna(subset=["gene_name", "RNAseq_log2FoldChange"])          # geni senza simbolo o non testati
    # 2. un solo gene per simbolo: il piu' espresso
    I = I.sort_values("RNAseq_baseMean", ascending=False).drop_duplicates("gene_name")
    padj_cut = a.padj_nacl if c_in.startswith("NaCl") else a.padj

    # 3. categorie con le soglie nuove (padj NaN = gene filtrato da DESeq2 -> non significativo)
    l, p, bm = I.RNAseq_log2FoldChange, I.RNAseq_padj.fillna(1), I.RNAseq_baseMean
    cat = pd.Series("other", index=I.index)                  # di partenza nessuna categoria
    cat[(p < padj_cut) & (l >= a.lfc)] = "UP"
    cat[(p < padj_cut) & (l <= -a.lfc)] = "DOWN"
    cat[(l.abs() < a.nc) & (bm >= a.min_bm)] = "NOCHANGE"
    gene2cat = dict(zip(I.gene_name, cat))                   # simbolo -> categoria

    # 4. intersezione con il BED dei metagene e scrittura di un BED per categoria
    row = {"confronto": c_out, "geni_BED_con_RNAseq": int(bed[3].isin(gene2cat).sum())}
    for k in CATS:
        sel = bed[bed[3].map(gene2cat) == k]                 # righe del BED la cui categoria e' k
        sel.to_csv(os.path.join(a.outdir, f"{c_out}_{k}.bed"), sep="\t", header=False, index=False)
        row[f"{k}_nel_BED"] = len(sel)                  # geni che finiscono nel metagene
        row[f"{k}_RNAseq_tot"] = int((cat == k).sum())   # geni della categoria in tutto l'RNA-seq
    rows.append(row)

# 5. riepilogo: quanti geni per categoria, prima e dopo l'intersezione con il BED (geni >= 5 kb)
S = pd.DataFrame(rows)
S.to_csv(os.path.join(a.outdir, "intersezione_summary.tsv"), sep="\t", index=False)
print(S.to_string(index=False))
