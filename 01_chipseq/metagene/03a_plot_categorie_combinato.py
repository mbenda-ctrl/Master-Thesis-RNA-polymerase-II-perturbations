#!/usr/bin/env python3
"""
Profilo di Pol II di controllo e trattato con le 3 categorie RNA-seq nello stesso plot: 6 profili.

Passaggi:
  1. segnale = IP - Input (bigWig gia' sottratti), in CPM
  2. per ogni gene e bin: media dei 2 replicati della condizione
  3. profilo = media sui geni della categoria (bin per bin)
  4. con --log: log2(profilo + 1), cioe' prima la media e poi il log (valori negativi portati a 0 prima del log)
Colore = condizione (come nei plot finali); sfumatura = categoria: UP scuro, NO CHANGE medio, DOWN chiaro.
Asse y fisso (--ymin / --ymax): se un profilo supera ymax lo script lo segnala.

Uso: python 03a_plot_categorie_combinato.py matrice.gz --ctrl Amr --treat Fast --out nome --title "..." --ymax 14 [--log]
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
ap.add_argument("matrix", help="matrice computeMatrix con 3 gruppi: UP, NOCHANGE, DOWN (in quest'ordine)")
ap.add_argument("--ctrl", required=True)
ap.add_argument("--treat", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--title", required=True)
ap.add_argument("--ymin", type=float, default=0)
ap.add_argument("--ymax", type=float, required=True)
ap.add_argument("--log", action="store_true")
ap.add_argument("--ylabel", default=None)
a = ap.parse_args()

# condizione -> (UP scuro, NO CHANGE medio, DOWN chiaro); il medio e' il colore dei plot finali
SHADES = {
    "Amr":   ("#145214", "#2ca02c", "#98df8a"),
    "Fast":  ("#7f0000", "#d62728", "#ff9896"),
    "Slow":  ("#3f2266", "#9467bd", "#c5b0d5"),
    "UT":    ("#0b3d6b", "#1f77b4", "#aec7e8"),
    "Storm": ("#994700", "#ff7f0e", "#ffbb78"),
    "UTpub": ("#333333", "#7f7f7f", "#c7c7c7"),
    "KCl":   ("#0a6b75", "#17becf", "#9edae5"),
}
CATS = ["UP", "NO CHANGE", "DOWN"]

# ------------------------------------------------------------------------------------------
# lettura della matrice (intestazione JSON + una riga per gene)
# ------------------------------------------------------------------------------------------
with gzip.open(a.matrix, "rt") as fh:
    h = json.loads(fh.readline()[1:])
X = pd.read_csv(a.matrix, sep="\t", skiprows=1, header=None, compression="gzip").iloc[:, 6:].to_numpy(float)
labels, sb, gb = h["sample_labels"], h["sample_boundaries"], h["group_boundaries"]
if len(gb) - 1 != 3:
    raise SystemExit(f"attesi 3 gruppi (UP, NOCHANGE, DOWN), trovati {len(gb) - 1}")
nb = sb[1] - sb[0]


def cond_mean(cond):
    """passaggio 2: media dei replicati della condizione, gene per gene e bin per bin"""
    idx = [i for i, s in enumerate(labels) if s.rpartition("_r")[0] == cond]
    if not idx:
        raise SystemExit(f"nessun campione {cond}_r* in {labels}")
    return np.mean([X[:, sb[i]:sb[i + 1]] for i in idx], axis=0)


# asse x: posizioni di TSS / TES in bin
bs = h["bin size"][0]; up = h["upstream"][0] // bs; dn = h["downstream"][0] // bs; body = h["body"][0] // bs
if body == 0:
    ticks, tlab, vl = [0, up, nb - 1], [f"-{up*bs/1000:g} kb", h["ref point"][0], f"+{dn*bs/1000:g} kb"], [up]
else:
    tes = up + body
    ticks, tlab, vl = [0, up, tes, nb - 1], [f"-{up*bs/1000:g} kb", "TSS", "TES", f"+{dn*bs/1000:g} kb"], [up, tes]

# ------------------------------------------------------------------------------------------
# passaggi 3-4: un profilo per condizione e categoria
# ------------------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.5, 5))
for cond in (a.ctrl, a.treat):
    G = cond_mean(cond)
    for g, cat in enumerate(CATS):
        rows = G[gb[g]:gb[g + 1]]
        if len(rows) == 0:
            continue
        prof = np.nanmean(rows, axis=0)
        if a.log:
            prof = np.log2(np.clip(prof, 0, None) + 1)
        if prof.max() > a.ymax:
            print(f"  ATTENZIONE: {cond} {cat} supera ymax ({prof.max():.2f} > {a.ymax:g})")
        ax.plot(np.arange(nb), prof, color=SHADES.get(cond, ("black",) * 3)[g], lw=1.8,
                label=f"{cond.replace('UTpub', 'UT*')} {cat} (n={len(rows)})")
for v in vl:
    ax.axvline(v, color="grey", lw=0.8, ls=":")
ax.set_xticks(ticks); ax.set_xticklabels(tlab); ax.set_xlim(0, nb - 1)
ax.set_ylim(a.ymin, a.ymax)
ax.set_ylabel(a.ylabel or ("log2(CPM + 1)" if a.log else "CPM"))
ax.set_title(a.title, fontsize=10)
ax.legend(frameon=False, loc="upper right", fontsize=8)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{a.out}.{ext}", dpi=300)
plt.close(fig)
print("salvato:", a.out + ".png")
