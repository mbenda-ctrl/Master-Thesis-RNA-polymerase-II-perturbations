#!/usr/bin/env python3
"""
METAGENE 1a - Profilo medio da una matrice computeMatrix, versione definitiva per la tesi
(sostituisce plotProfile di deepTools per avere colori, assi e scala log sotto controllo).

- un colore per condizione, replicato 1 scuro / replicato 2 chiaro
- asse Y sempre LINEARE
- con --log il logaritmo e' applicato al PROFILO MEDIO (prima la media, poi il log):
      log2(max(media, 0) + 1)
  come nei plot in scala log precedenti; con --log-mode genes il log e' applicato
  gene per gene prima della media
- etichette dell'asse X lette dall'header della matrice (reference-point o scale-regions)

Uso:
  python 01a_plot_metagene_profile.py matrice.gz --out nome --title TSS --ymin 0 --ymax 10 [--log]
"""
import argparse, gzip, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("matrix")
ap.add_argument("--out", required=True, help="nome file senza estensione")
ap.add_argument("--title", required=True)
ap.add_argument("--log", action="store_true", help="logCPM: log2(CPM+1)")
ap.add_argument("--log-axis", action="store_true",
                help="come i vecchi plot: profilo in CPM su asse Y logaritmico (tacche in CPM)")
ap.add_argument("--log-mode", choices=["profile", "genes"], default="profile",
                help="profile = media poi log (default); genes = log gene per gene poi media")
ap.add_argument("--ymin", type=float, default=0)
ap.add_argument("--ymax", default="auto", help="numero, oppure 'auto' = 10%% sopra il massimo delle matrici in --scale-from")
ap.add_argument("--scale-from", nargs="*", default=None,
                help="matrici da cui calcolare un ymax comune (default: solo questa)")
ap.add_argument("--ylabel", default=None, help="default: CPM o logCPM")
ap.add_argument("--start-label", default="TSS")
ap.add_argument("--end-label", default="TES")
a = ap.parse_args()

# ---------- colori: condizione -> (replicato 1 scuro, replicato 2 chiaro) ----------
PAIRS = {
    "Amr":   ("#2ca02c", "#98df8a"),   # verde
    "Fast":  ("#d62728", "#ff9896"),   # rosso
    "Slow":  ("#9467bd", "#c5b0d5"),   # viola
    "UT":    ("#1f77b4", "#aec7e8"),   # blu
    "Storm": ("#ff7f0e", "#ffbb78"),   # arancio
    "UT*":   ("#7f7f7f", "#c7c7c7"),   # grigio
    "KCl":   ("#17becf", "#9edae5"),   # ciano
}


def color_of(name):
    """colore di un campione: Amr_r1 -> condizione Amr, replicato 1 -> tono scuro"""
    cond, _, rep = name.rpartition("_r")
    pair = PAIRS.get(cond, ("black", "grey"))
    return pair[0] if rep == "1" else pair[1]


def profiles(path):
    """profili medi (lineari o logCPM) di tutti i campioni di una matrice."""
    with gzip.open(path, "rt") as fh:
        hh = json.loads(fh.readline()[1:])                # intestazione JSON della matrice
    M = pd.read_csv(path, sep="\t", skiprows=1, header=None, compression="gzip").iloc[:, 6:].to_numpy(float)   # solo i bin
    if a.log and a.log_mode == "genes":
        M = np.log2(np.clip(M, 0, None) + 1)
    out = []
    for i in range(len(hh["sample_labels"])):
        p = np.nanmean(M[:, hh["sample_boundaries"][i]:hh["sample_boundaries"][i + 1]], axis=0)   # media sui geni
        if a.log and a.log_mode == "profile":
            p = np.log2(np.clip(p, 0, None) + 1)
        out.append(p)
    return out


def nice_ceil(x):
    """arrotonda per eccesso a 1, 2, 2.5 o 5 x 10^n."""
    e = 10 ** np.floor(np.log10(x))
    for m in (1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10):
        if m * e >= x:
            return m * e


# ---------- matrice ----------
with gzip.open(a.matrix, "rt") as fh:
    h = json.loads(fh.readline()[1:])                    # prima riga: "@" + intestazione JSON
X = pd.read_csv(a.matrix, sep="\t", skiprows=1, header=None, compression="gzip").iloc[:, 6:].to_numpy(float)
if a.log and a.log_mode == "genes":
    X = np.log2(np.clip(X, 0, None) + 1)
labels, sb = h["sample_labels"], h["sample_boundaries"]   # nomi dei campioni e colonne dove iniziano

# misure della matrice espresse in numero di bin
bs = h["bin size"][0]                                # bin (bp)
up = h["upstream"][0] // bs                          # bin a monte
dn = h["downstream"][0] // bs                        # bin a valle
body = h["body"][0] // bs                            # bin del corpo scalato (0 = reference-point)
# parti non scalate opzionali (--unscaled5prime / --unscaled3prime di computeMatrix; 0 se assenti)
u5 = h.get("unscaled 5 prime", [0])[0] // bs
u3 = h.get("unscaled 3 prime", [0])[0] // bs
nb = sb[1] - sb[0]                                   # bin totali per campione

if body == 0:                                   # reference-point
    ref = h["ref point"][0] or "ref"
    ticks = [0, up, nb - 1]
    tlab = [f"-{up*bs/1000:g} kb", ref, f"+{dn*bs/1000:g} kb"]
    vlines = [up]
else:                                           # scale-regions
    tes = up + u5 + body + u3
    ticks = [0, up, tes, nb - 1]
    tlab = [f"-{up*bs/1000:g} kb", a.start_label, a.end_label, f"+{dn*bs/1000:g} kb"]
    vlines = [up, tes]

# ---------- figura ----------
fig, ax = plt.subplots(figsize=(6.5, 5))
for i, s in enumerate(labels):
    prof = np.nanmean(X[:, sb[i]:sb[i + 1]], axis=0)   # profilo medio del campione i
    if a.log and a.log_mode == "profile":
        prof = np.log2(np.clip(prof, 0, None) + 1)
    ax.plot(np.arange(nb), prof, color=color_of(s), lw=1.8, label=s)   # una linea per campione
for v in vlines:
    ax.axvline(v, color="grey", lw=0.8, ls=":")
ax.set_xticks(ticks)
ax.set_xticklabels(tlab)
ax.set_xlim(0, nb - 1)
if a.ymax == "auto":
    # ymax comune: picco piu' alto tra tutte le matrici indicate, +10%, arrotondato
    peak = max(p.max() for m in (a.scale_from or [a.matrix]) for p in profiles(m))
    ymax = nice_ceil(peak * 1.1)
    print(f"  ymax automatico = {ymax:g} (picco {peak:.2f})")
else:
    ymax = float(a.ymax)
if a.log_axis:
    from matplotlib.ticker import ScalarFormatter
    ax.set_yscale("log")
    ax.set_ylim(None, ymax)                     # il minimo deve essere > 0: automatico
    ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.set_ylabel(a.ylabel or "CPM (scala log)")
else:
    ax.set_ylim(a.ymin, ymax)
    ax.set_ylabel(a.ylabel or ("logCPM" if a.log else "CPM"))
ax.set_title(a.title)
ax.legend(frameon=False, loc="upper right")          # legenda senza riquadro
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
for ext in ("png", "pdf"):                           # PNG e PDF
    fig.savefig(f"{a.out}.{ext}", dpi=300)
print("salvato:", a.out + ".png")
