#!/usr/bin/env python3
"""
Volcano plot da un report DiffBind (dba.report esportato in CSV/TSV): log2FC (colonna Fold)
sull'asse x, -log10(FDR) sull'asse y, colorando le regioni significative up/down.

Uso:
  python diffbind_volcano.py report.tsv --nameA Fast --nameB Amr --out volcano_fast_amr \
      [--fdr 0.05] [--label-top 10]
"""
import argparse
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("report")
ap.add_argument("--nameA", required=True, help="nome del gruppo con Fold>0 (es. Fast)")
ap.add_argument("--nameB", required=True, help="nome del gruppo di riferimento, Fold<0 (es. Amr)")
ap.add_argument("--out", required=True)
ap.add_argument("--fdr", type=float, default=0.05)
ap.add_argument("--min-abs-fold", type=float, default=0.0, help="soglia su |log2FC| (es. 1.5): sopra = colorato, sotto = grigio")
ap.add_argument("--label-top", type=int, default=0, help="etichetta le N regioni piu' significative/estreme (0=nessuna)")
ap.add_argument("--label-by", choices=["fdr", "fold"], default="fdr", help="criterio di selezione dei top hit da etichettare: 'fdr' (FDR piu' basso) o 'fold' (|log2FC| piu' alto)")
ap.add_argument("--bed", default=None, help="BED con nome gene in colonna 4 (es. genes_pc_union_MANE_n2_min1kb.bed): se fornito, le etichette dei top hit mostrano il nome del gene invece delle coordinate")
ap.add_argument("--ymax", type=float, default=None, help="limite superiore asse Y (-log10 FDR), per uniformare la scala tra plot diversi")
ap.add_argument("--xmax", type=float, default=None, help="limite assoluto asse X (|log2FC|), applicato simmetricamente")
a = ap.parse_args()

sep = "," if a.report.endswith(".csv") else "\t"
D = pd.read_csv(a.report, sep=sep)
D.columns = [c.strip('"') for c in D.columns]
D["negLog10FDR"] = -np.log10(D["FDR"].clip(lower=1e-300))

sig_up = (D.FDR < a.fdr) & (D.Fold > a.min_abs_fold)
sig_down = (D.FDR < a.fdr) & (D.Fold < -a.min_abs_fold)
notsig = ~(sig_up | sig_down)

fig, ax = plt.subplots(figsize=(8, 6))
if notsig.sum() > 0:
    ax.scatter(D.Fold[notsig], D.negLog10FDR[notsig], s=8, color="lightgrey", alpha=0.6, linewidths=0, label=f"non signif. (n={notsig.sum()})")
ax.scatter(D.Fold[sig_down], D.negLog10FDR[sig_down], s=10, color="#3A6EA5", alpha=0.7, linewidths=0, label=f"{a.nameB} > {a.nameA} (n={sig_down.sum()})")
ax.scatter(D.Fold[sig_up], D.negLog10FDR[sig_up], s=10, color="#C1666B", alpha=0.7, linewidths=0, label=f"{a.nameA} > {a.nameB} (n={sig_up.sum()})")

ax.axhline(-np.log10(a.fdr), color="grey", lw=0.8, ls="--")
ax.axvline(0, color="grey", lw=0.8, ls="--")
if a.min_abs_fold > 0:
    for s in (-a.min_abs_fold, a.min_abs_fold):
        ax.axvline(s, color="grey", lw=0.8, ls=":")
ax.set_xlabel(f"log2FC ({a.nameA} / {a.nameB})")
if a.ymax is not None:
    ax.set_ylim(top=a.ymax)
if a.xmax is not None:
    ax.set_xlim(-a.xmax, a.xmax)
ax.set_ylabel("-log10(FDR)")
ax.set_title(f"Volcano plot: {a.nameA} vs {a.nameB}", fontsize=12)
ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
ax.spines[["top","right"]].set_visible(False)

if a.label_top > 0:
    if a.label_by == "fold":
        top = D.reindex(D.Fold.abs().sort_values(ascending=False).index[:a.label_top])
    else:
        top = D.nsmallest(a.label_top, "FDR")

    gene_of = {}
    if a.bed:
        bed = pd.read_csv(a.bed, sep="\t", header=None, usecols=[0,1,2,3,5], names=["chrom","start","end","name","strand"])
        # finestra TU: corpo genico + 3kb a monte del TSS (i picchi cadono spesso nel promotore,
        # non dentro il gene stesso -- stessa logica gia' usata altrove nel progetto per assegnare i picchi ai geni)
        bed["tu_start"] = np.where(bed.strand == "+", bed.start - 3000, bed.start)
        bed["tu_end"] = np.where(bed.strand == "+", bed.end, bed.end + 3000)
        bed["tu_start"] = bed["tu_start"].clip(lower=0)
        for chrom, g in bed.groupby("chrom"):
            rows = D.index[(D.seqnames == chrom)]
            if len(rows) == 0: continue
            sub = D.loc[rows]
            for idx, row in zip(rows, sub.itertuples()):
                ov = g[(g.tu_end > row.start) & (g.tu_start < row.end)]
                if len(ov) > 0:
                    gene_of[idx] = "/".join(ov.name.unique()[:2])

    texts = []
    for idx, row in top.iterrows():
        label = gene_of.get(idx, f"{row.seqnames}:{row.start}")
        t = ax.annotate(label, (row.Fold, row.negLog10FDR),
                        fontsize=6.5, ha="center", style="italic" if idx in gene_of else "normal")
        texts.append(t)

    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))
    except ImportError:
        pass  # se adjustText non e' installato, le etichette restano nella posizione di default (possibile sovrapposizione)

fig.tight_layout()
fig.savefig(f"{a.out}.png", dpi=200, bbox_inches="tight")
print(f"{a.nameA} vs {a.nameB}: {len(D)} regioni testate | signif. FDR<{a.fdr}: {sig_up.sum()+sig_down.sum()} ({sig_up.sum()} su, {sig_down.sum()} giu)")
print(f"salvato: {a.out}.png")
