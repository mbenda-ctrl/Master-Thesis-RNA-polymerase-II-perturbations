#!/usr/bin/env python3
"""
Boxplot del log2FC di Pol II (trattato / controllo) per gene, in 3 regioni, separato per categoria RNA-seq.

Regioni (ognuna presa dalla matrice adatta):
  promotore     = TSS +- 300 bp                                  (matrice TSS, reference-point)
  gene body     = corpo scalato, escluso il primo e l'ultimo 10%  (matrice BODY, scale-regions)
  terminazione  = da TES a TES + 3 kb                            (matrice TES, reference-point)

Passaggi per ogni gene e regione:
  1. segnale di ogni replicato = media dei bin della finestra (IP - Input, CPM)
  2. media dei 2 replicati della condizione; valori negativi portati a 0
  3. log2FC = log2((trattato + ps) / (controllo + ps)), ps = 0.1
Test: Mann-Whitney (due code) di UP contro NO CHANGE e di DOWN contro NO CHANGE in ogni regione;
il p-value e' scritto sopra il box. Le mediane, n e p finiscono anche nel file .tsv.

Uso: python 03c_plot_box_categorie.py --tss m_TSS.gz --body m_BODY.gz --tes m_TES.gz --ctrl Amr --treat Fast \
        --out nome --title "..." [--ylim 3]
"""
import argparse, gzip, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# ------------------------------------------------------------------------------------------
# parametri
# ------------------------------------------------------------------------------------------
ap = argparse.ArgumentParser()
ap.add_argument("--tss", required=True); ap.add_argument("--body", required=True); ap.add_argument("--tes", required=True)
ap.add_argument("--ctrl", required=True); ap.add_argument("--treat", required=True)
ap.add_argument("--out", required=True); ap.add_argument("--title", required=True)
ap.add_argument("--ylim", type=float, default=3, help="asse y fisso: da -ylim a +ylim")
ap.add_argument("--pseudo", type=float, default=0.1)
a = ap.parse_args()

CATS = [("UP", "#d62728"), ("NO CHANGE", "#7f7f7f"), ("DOWN", "#1f77b4")]   # stesso ordine dei gruppi nella matrice
LAB = lambda s: s.replace("UTpub", "UT*")


def region_lfc(path, which):
    """log2FC per gene nella finestra della regione; restituisce una lista [UP, NOCHANGE, DOWN] di array"""
    with gzip.open(path, "rt") as fh:
        h = json.loads(fh.readline()[1:])
    X = pd.read_csv(path, sep="\t", skiprows=1, header=None, compression="gzip").iloc[:, 6:].to_numpy(float)
    lab, sb, gb = h["sample_labels"], h["sample_boundaries"], h["group_boundaries"]
    bs = h["bin size"][0]; up = h["upstream"][0] // bs; body = h["body"][0] // bs
    # finestra in bin
    if which == "promotore":
        w = slice(up - 300 // bs, up + 300 // bs)                       # TSS +- 300 bp
    elif which == "genebody":
        w = slice(up + int(0.1 * body), up + int(0.9 * body))          # 10-90% del corpo scalato
    else:
        w = slice(up, sb[1] - sb[0])                                    # TES -> fine della finestra (+3 kb)

    def cond(c):
        # passaggi 1-2: media dei bin della finestra per replicato, poi media dei replicati, negativi -> 0
        idx = [i for i, s in enumerate(lab) if s.rpartition("_r")[0] == c]
        return np.clip(np.mean([X[:, sb[i]:sb[i + 1]][:, w].mean(axis=1) for i in idx], axis=0), 0, None)

    L = np.log2((cond(a.treat) + a.pseudo) / (cond(a.ctrl) + a.pseudo))   # passaggio 3
    return [L[gb[g]:gb[g + 1]] for g in range(3)]


REG = [("promotore\n(TSS ± 300 bp)", region_lfc(a.tss, "promotore")),
       ("gene body\n(10-90%)", region_lfc(a.body, "genebody")),
       ("terminazione\n(TES → +3 kb)", region_lfc(a.tes, "terminazione"))]

# ------------------------------------------------------------------------------------------
# figura: per ogni regione 3 box (UP, NO CHANGE, DOWN) + test contro NO CHANGE
# ------------------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
rows = []
for r, (rname, groups) in enumerate(REG):
    base = r * 4                                         # 3 box + uno spazio per regione
    for g, (cat, col) in enumerate(CATS):
        if len(groups[g]) == 0:
            rows.append({"regione": rname.split("\n")[0], "categoria": cat, "n": 0, "mediana_log2FC": np.nan})
            continue
        bp = ax.boxplot(groups[g], positions=[base + g], widths=0.8, patch_artist=True, showfliers=False,
                        medianprops=dict(color="black"))
        bp["boxes"][0].set(facecolor=col, alpha=0.6)
        rows.append({"regione": rname.split("\n")[0], "categoria": cat, "n": len(groups[g]),
                     "mediana_log2FC": np.median(groups[g])})
    for g in (0, 2):                                     # UP e DOWN contro NO CHANGE
        if len(groups[g]) == 0 or len(groups[1]) == 0:
            continue
        p = mannwhitneyu(groups[g], groups[1]).pvalue
        rows[-3 + g]["p_vs_NC"] = p
        ax.text(base + g, a.ylim * 0.9, f"p={p:.1e}", ha="center", fontsize=7, rotation=90, va="top")
ax.axhline(0, color="black", lw=0.8, ls="--")
ax.set_xticks([r * 4 + 1 for r in range(3)]); ax.set_xticklabels([r[0] for r in REG])
ax.set_ylim(-a.ylim, a.ylim)
ax.set_ylabel(f"log2FC Pol II ({LAB(a.treat)} / {LAB(a.ctrl)}) per gene")
ax.set_title(a.title, fontsize=10)
ax.legend([plt.Rectangle((0, 0), 1, 1, fc=c, alpha=0.6) for _, c in CATS],
          [f"{c} (n={len(REG[0][1][i])})" for i, (c, _) in enumerate(CATS)], frameon=False, loc="lower right", fontsize=8)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{a.out}.{ext}", dpi=300)
plt.close(fig)
pd.DataFrame(rows).to_csv(f"{a.out}.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).to_string(index=False))
print("salvato:", a.out + ".png")
