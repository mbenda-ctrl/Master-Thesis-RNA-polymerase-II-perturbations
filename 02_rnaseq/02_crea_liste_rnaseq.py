#!/usr/bin/env python3
"""
Liste dei geni differenziali RNA-seq (UP / DOWN) dalle tabelle complete di 01_deseq2_rnaseq.R.
Criteri: padj (BH) < 0.05 (NaCl vs UT: < 0.10) e |log2FC MLE| >= 0.58 (fold change 1.5x).
Uso: python 02_crea_liste_rnaseq.py <cartella_RNAseq> <cartella_liste>
"""
import sys, os, pandas as pd
IN, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
CONF = {"Slow_vs_Amr": ("Slow_vs_Amr", 0.05), "Fast_vs_Amr": ("Fast_vs_Amr", 0.05),
        "Storm_vs_UT": ("STORM_vs_WT", 0.05), "NaCl_vs_UT": ("NACL_vs_WT", 0.10)}
for tab, (nome, padj) in CONF.items():
    A = pd.read_csv(f"{IN}/{tab}_all_genes.csv")
    A = A.set_index("gene_id")[["baseMean", "log2FC_MLE", "lfcSE", "stat", "pvalue", "padj"]].rename(columns={"log2FC_MLE": "log2FoldChange"})
    A.index.name = ""
    sig = A[A.padj < padj]
    for d, m in (("UP", sig.log2FoldChange >= 0.58), ("DOWN", sig.log2FoldChange <= -0.58)):
        sig[m].sort_values("padj").to_csv(f"{OUT}/{nome}_{d}.csv")
    print(nome, "UP", int((sig.log2FoldChange >= 0.58).sum()), "DOWN", int((sig.log2FoldChange <= -0.58).sum()))
