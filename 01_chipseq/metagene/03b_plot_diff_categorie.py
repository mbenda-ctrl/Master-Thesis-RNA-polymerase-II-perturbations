#!/usr/bin/env python3
"""
Profilo DIFFERENZIALE di Pol II per categoria RNA-seq: una curva per UP, NO CHANGE e DOWN.

Passaggi (per ogni gene e per ogni bin della matrice di computeMatrix):
  1. segnale = IP - Input (bigWig gia' sottratti); valori negativi portati a 0 (sotto l'input = nessun segnale)
  2. media dei 2 replicati del trattato e media dei 2 replicati del controllo
  3. log2FC = log2((trattato + ps) / (controllo + ps)), ps = pseudoconto (default 0.1 CPM) per evitare log(0)
  4. media mobile su --smooth bin (default 5) lungo il profilo di ogni gene, per ridurre il rumore
  5. curva = media del log2FC sui geni della categoria (media dei log2, non log2 della media)
  6. banda = intervallo di confidenza 95% bootstrap: si ricampionano i geni con reinserimento --nboot volte
     e si prendono il 2.5 e il 97.5 percentile delle medie
La curva NO CHANGE fa da riferimento: lo spostamento globale del segnale tra le condizioni
(es. diversa efficienza dell'IP, senza spike-in) compare anche li'; conta la distanza di UP e DOWN da NO CHANGE.

Uso: python 03b_plot_diff_categorie.py matrice.gz --ctrl Amr --treat Fast --out nome --title "..." [--ylim 3]
"""
import argparse, gzip, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------------------
# parametri
# ------------------------------------------------------------------------------------------
ap = argparse.ArgumentParser()
ap.add_argument("matrix", help="matrice computeMatrix con 3 gruppi di regioni: UP, NOCHANGE, DOWN (in quest'ordine)")
ap.add_argument("--ctrl", required=True, help="prefisso dei campioni di controllo (es. Amr -> Amr_r1, Amr_r2)")
ap.add_argument("--treat", required=True, help="prefisso dei campioni trattati")
ap.add_argument("--out", required=True, help="nome dei file in uscita, senza estensione")
ap.add_argument("--title", required=True)
ap.add_argument("--ylim", type=float, default=3, help="asse y fisso e simmetrico: da -ylim a +ylim")
ap.add_argument("--pseudo", type=float, default=0.1)
ap.add_argument("--smooth", type=int, default=5, help="media mobile in bin (1 = nessuna)")
ap.add_argument("--nboot", type=int, default=500)
a = ap.parse_args()

CATS = [("UP", "#d62728"), ("NO CHANGE", "#7f7f7f"), ("DOWN", "#1f77b4")]   # stesso ordine dei gruppi nella matrice
LAB = lambda s: s.replace("UTpub", "UT*")                                   # nome per le etichette

# ------------------------------------------------------------------------------------------
# lettura della matrice: prima riga = intestazione JSON (campioni, confini, parametri), poi una riga per gene
# (6 colonne BED + i bin di tutti i campioni affiancati)
# ------------------------------------------------------------------------------------------
with gzip.open(a.matrix, "rt") as fh:
    h = json.loads(fh.readline()[1:])
X = pd.read_csv(a.matrix, sep="\t", skiprows=1, header=None, compression="gzip").iloc[:, 6:].to_numpy(float)
labels = h["sample_labels"]         # es. Amr_r1, Amr_r2, Fast_r1, Fast_r2
sb = h["sample_boundaries"]         # colonne dove inizia/finisce ogni campione
gb = h["group_boundaries"]          # righe dove inizia/finisce ogni gruppo (UP, NOCHANGE, DOWN)
nb = sb[1] - sb[0]                  # numero di bin per campione


def cond_mean(cond):
    """media dei replicati di una condizione (bin per bin, gene per gene), con IP - Input < 0 portato a 0"""
    idx = [i for i, s in enumerate(labels) if s.rpartition("_r")[0] == cond]
    if not idx:
        raise SystemExit(f"nessun campione {cond}_r* in {labels}")
    return np.clip(np.mean([X[:, sb[i]:sb[i + 1]] for i in idx], axis=0), 0, None)


# passaggi 1-3: log2FC per gene e per bin
L = np.log2((cond_mean(a.treat) + a.pseudo) / (cond_mean(a.ctrl) + a.pseudo))
# passaggio 4: media mobile lungo il profilo di ogni gene
if a.smooth > 1:
    k = np.ones(a.smooth) / a.smooth
    L = np.apply_along_axis(lambda r: np.convolve(r, k, mode="same"), 1, L)

# asse x: posizioni di TSS / TES (o del punto di riferimento) in bin
bs = h["bin size"][0]; up = h["upstream"][0] // bs; dn = h["downstream"][0] // bs; body = h["body"][0] // bs
if body == 0:      # matrice reference-point (TSS o TES)
    ticks, tlab, vl = [0, up, nb - 1], [f"-{up*bs/1000:g} kb", h["ref point"][0], f"+{dn*bs/1000:g} kb"], [up]
else:              # matrice scale-regions (gene completo)
    tes = up + body
    ticks, tlab, vl = [0, up, tes, nb - 1], [f"-{up*bs/1000:g} kb", "TSS", "TES", f"+{dn*bs/1000:g} kb"], [up, tes]

# ------------------------------------------------------------------------------------------
# passaggi 5-6: curva media e banda bootstrap per ogni categoria
# ------------------------------------------------------------------------------------------
rng = np.random.default_rng(1)       # seme fisso: stesse bande a ogni esecuzione
fig, ax = plt.subplots(figsize=(6.5, 5))
x = np.arange(nb)
for g, (cat, col) in enumerate(CATS):
    M = L[gb[g]:gb[g + 1]]           # geni della categoria
    if len(M) == 0:
        continue
    if len(M) < 5:
        print(f"  ATTENZIONE: {cat} ha solo {len(M)} geni")
    mean = M.mean(axis=0)
    boot = np.array([M[rng.integers(0, len(M), len(M))].mean(axis=0) for _ in range(a.nboot)])
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
    ax.fill_between(x, lo, hi, color=col, alpha=0.2, lw=0)
    ax.plot(x, mean, color=col, lw=1.8, label=f"{cat} (n={len(M)})")
ax.axhline(0, color="black", lw=0.8, ls="--")
for v in vl:
    ax.axvline(v, color="grey", lw=0.8, ls=":")
ax.set_xticks(ticks); ax.set_xticklabels(tlab); ax.set_xlim(0, nb - 1)
ax.set_ylim(-a.ylim, a.ylim)
ax.set_ylabel(f"log2FC Pol II ({LAB(a.treat)} / {LAB(a.ctrl)})")
ax.set_title(a.title, fontsize=10)
ax.legend(frameon=False, loc="upper right")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{a.out}.{ext}", dpi=300)
plt.close(fig)
print("salvato:", a.out + ".png")
