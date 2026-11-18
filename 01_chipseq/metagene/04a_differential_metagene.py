#!/usr/bin/env python3
"""
METAGENE 4a - Costruisce un profilo e una heatmap di log2 fold-change tra due gruppi di campioni,
a partire da una matrice di computeMatrix che contiene SOLO i campioni dei due gruppi
(subset con computeMatrixOperations subset, o ricalcolata con -S solo su quei bigwig).

Uso:
  python 04a_differential_metagene.py matrix.gz --groupA IP_E11_1 IP_E22_2 --groupB IP_Amr_1 IP_Amr_2 \
      --nameA Fast --nameB Amr --out prefix [--pseudocount 0.1]

Produce:
  <prefix>_log2fc_profile.png   profilo medio log2(A/B) lungo la regione
  <prefix>_log2fc_heatmap.png   heatmap log2(A/B) per regione (righe = regioni, ordinate per media)
  <prefix>_log2fc.tsv           matrice numerica dei valori log2FC per regione x bin
  <prefix>_length_vs_log2fc_scatter.png   (con --bed) lunghezza del gene vs log2FC medio

Calcolo: per ogni gene e bin, media dei replicati di A e di B, poi log2((A + ps) / (B + ps)).
Il profilo e' la media dei log2FC sui geni (media dei log, non log della media).
"""
import argparse, gzip, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

ap = argparse.ArgumentParser()
ap.add_argument("matrix"); ap.add_argument("--groupA", nargs="+", required=True); ap.add_argument("--groupB", nargs="+", required=True)
ap.add_argument("--nameA", default="A"); ap.add_argument("--nameB", default="B")
ap.add_argument("--out", required=True); ap.add_argument("--pseudocount", type=float, default=0.1)
ap.add_argument("--ymin", type=float, default=None, help="limite inferiore asse y del profilo (per confrontare piu' plot sulla stessa scala)")
ap.add_argument("--ymax", type=float, default=None, help="limite superiore asse y del profilo")
ap.add_argument("--zmax", type=float, default=None, help="valore massimo fisso per la scala colore della heatmap (default: 95° percentile dei dati, calcolato automaticamente)")
ap.add_argument("--bed", default=None, help="stesso file BED passato a computeMatrix -R: se fornito, aggiunge una barra con la lunghezza genica accanto alla heatmap")
ap.add_argument("--length-bins", default="1,2,3,5,10,20,50,100,200",
                help="soglie (in kb) per gli intervalli discreti della barra lunghezza, separate da virgola (default: 1,2,3,5,10,20,50,100,200 kb)")
a = ap.parse_args()

# intestazione JSON della matrice (prima riga, preceduta da "@")
with gzip.open(a.matrix, "rt") as fh: h = json.loads(fh.readline().lstrip("@"))
# valori: una riga per gene, dopo le 6 colonne BED i bin di tutti i campioni affiancati
X = pd.read_csv(a.matrix, sep="\t", skiprows=1, header=None, compression="gzip").iloc[:, 6:].to_numpy(float)
labels = h["sample_labels"]; sb = h["sample_boundaries"]     # nomi dei campioni e colonne di inizio/fine
# ogni campione richiesto deve essere nella matrice
for s in a.groupA + a.groupB: assert s in labels, f"campione {s} non nella matrice: {labels}"

def group_mean(names):
    idx = [labels.index(n) for n in names]
    arrs = [X[:, sb[i]:sb[i+1]] for i in idx]
    return np.nanmean(np.stack(arrs, axis=0), axis=0)   # media tra repliche, per regione x bin

Am = group_mean(a.groupA); Bm = group_mean(a.groupB)        # media dei replicati, geni x bin
L2FC = np.log2((Am + a.pseudocount) / (Bm + a.pseudocount))   # log2FC per gene e per bin

# misure in numero di bin; body = 0 vuol dire matrice reference-point
bs = h["bin size"][0]; up = h["upstream"][0] // bs; body = h["body"][0] // bs; down = h["downstream"][0] // bs
isref = body == 0

# --- profilo medio ---
prof = np.nanmean(L2FC, axis=0)                              # media sui geni, bin per bin
x = np.arange(len(prof))
fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(x, prof, color="black", lw=1.8)
ax.axhline(0, color="grey", lw=0.8, ls="--")
if isref:
    refpoint = h.get("ref point", ["center"])
    refpoint = refpoint[0] if isinstance(refpoint, list) else refpoint
    refpoint = {"TSS": "TSS", "TES": "TES", "center": "center"}.get(refpoint, refpoint or "center")
    ax.axvline(up, color="grey", lw=0.8)
    ax.set_xticks([0, up, len(prof)-1]); ax.set_xticklabels([f"-{up*bs/1000:g} kb", refpoint, f"+{down*bs/1000:g} kb"])
else:
    ax.axvline(up, color="grey", lw=0.8, ls=":"); ax.axvline(up+body, color="grey", lw=0.8, ls=":")
    # "start"/"end" del corpo scalato = TSS e TES di ogni regione (gene), resi espliciti
    ax.set_xticks([0, up, up+body, len(prof)-1]); ax.set_xticklabels([f"-{up*bs/1000:g} kb", "TSS", "TES", f"+{down*bs/1000:g} kb"])
refpoint_title = refpoint if isref else "TSS\u2192TES"
if a.ymin is not None or a.ymax is not None: ax.set_ylim(a.ymin, a.ymax)   # asse y fisso se richiesto
ax.set_ylabel(f"log2FC ({a.nameA} / {a.nameB})"); ax.set_title(f"Differential profile: {a.nameA} vs {a.nameB} ({refpoint_title}, n={len(L2FC)} regioni)")
ax.spines[["top","right"]].set_visible(False); fig.tight_layout(); fig.savefig(f"{a.out}_log2fc_profile.png", dpi=200)

# --- heatmap, ordinata per media log2FC della regione (decrescente) ---
order = np.argsort(-np.nanmean(L2FC, axis=1))
vmax = a.zmax if a.zmax is not None else np.nanpercentile(np.abs(L2FC), 95)   # scala colore simmetrica
n_genes = len(L2FC)

if a.bed:
    bed = pd.read_csv(a.bed, sep="\t", header=None, usecols=[1,2], names=["start","end"])
    lengths = (bed["end"] - bed["start"]).to_numpy()             # lunghezza di ogni gene (bp)
    assert len(lengths) == n_genes, f"il BED ha {len(lengths)} righe ma la matrice ne ha {n_genes}: devono corrispondere 1:1 (stesso file passato a computeMatrix -R, con --sortRegions keep)"

    # GridSpec esplicito: ogni asse-dato ha la SUA riga di altezza garantita identica (sharey),
    # le colorbar vivono in assi propri e non "rubano" mai spazio agli assi con i dati (bug della v4).
    fig2 = plt.figure(figsize=(11.5, 7.5))
    gs = fig2.add_gridspec(2, 4, width_ratios=[6, 0.25, 0.9, 1.6], height_ratios=[20, 1], hspace=0.35, wspace=0.15)
    ax2 = fig2.add_subplot(gs[0, 0])
    cax_heat = fig2.add_subplot(gs[1, 0])
    ax3 = fig2.add_subplot(gs[0, 2], sharey=ax2)
    cax_len = fig2.add_subplot(gs[1, 2])
    ax4 = fig2.add_subplot(gs[0, 3], sharey=ax2)
else:
    fig2, ax2 = plt.subplots(figsize=(6, 7))

im = ax2.imshow(L2FC[order], aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)   # righe ordinate
if isref:
    ax2.axvline(up, color="black", lw=0.8)
    ax2.set_xticks([0, up, L2FC.shape[1]-1]); ax2.set_xticklabels([f"-{up*bs/1000:g} kb", refpoint, f"+{down*bs/1000:g} kb"])
else:
    ax2.axvline(up, color="black", lw=0.8, ls=":"); ax2.axvline(up+body, color="black", lw=0.8, ls=":")
    ax2.set_xticks([0, up, up+body, L2FC.shape[1]-1]); ax2.set_xticklabels([f"-{up*bs/1000:g} kb", "TSS", "TES", f"+{down*bs/1000:g} kb"])
ax2.set_ylabel("regioni (geni)"); ax2.set_xlabel(f"posizione lungo il metagene ({refpoint_title})"); ax2.set_title(f"log2FC {a.nameA}/{a.nameB}")
ax2.set_ylim(n_genes - 0.5, -0.5)  # garantisce esplicitamente lo stesso range su tutti gli assi condivisi

if a.bed:
    fig2.colorbar(im, cax=cax_heat, orientation="horizontal", label="log2FC")

    len_kb = lengths[order] / 1000                               # lunghezze nello stesso ordine della heatmap
    len_kb_safe = np.clip(len_kb, 0.1, None)  # LogNorm non accetta zero

    # --- pannello 1: striscia colorata, scala LOGARITMICA (i geni coprono ordini di grandezza diversi) ---
    norm_log = mcolors.LogNorm(vmin=max(len_kb_safe.min(), 0.5), vmax=len_kb_safe.max())
    im3 = ax3.imshow(len_kb_safe.reshape(-1, 1), aspect="auto", cmap="viridis", norm=norm_log,
                     extent=[0, 1, n_genes - 0.5, -0.5], interpolation="nearest")
    ax3.set_xticks([]); ax3.set_title("lunghezza\ngeni", fontsize=9)
    ax3.spines[["top","right","left","bottom"]].set_visible(False)
    ax3.tick_params(left=False, labelleft=False, bottom=False)
    cbar3 = fig2.colorbar(im3, cax=cax_len, orientation="horizontal", label="lunghezza (kb, scala log)")
    cbar3.set_ticks([1, 10, 100, 1000]); cbar3.set_ticklabels(["1", "10", "100", "1000"])
    cbar3.ax.tick_params(labelsize=8)

    # --- pannello 2: mediana mobile di log10(lunghezza) lungo l'ordinamento -- il vero strumento
    #     diagnostico: mostra un trend anche dove la striscia colorata soffre di aliasing ---
    win = min(200, max(10, n_genes // 20))
    roll = pd.Series(np.log10(len_kb_safe)).rolling(win, center=True, min_periods=max(5, win//4)).median()
    ax4.plot(roll.values, np.arange(n_genes), color="#2E5C3E", lw=1.5)
    ax4.set_xlabel(f"mediana mobile\nlog10(kb), finestra={win}"); ax4.set_title("trend\nlunghezza", fontsize=9)
    ax4.spines[["top","right"]].set_visible(False)
    ax4.tick_params(left=False, labelleft=False)

    # --- correlazione Spearman: lunghezza (log10) vs log2FC medio del gene -- il numero da citare ---
    from scipy.stats import spearmanr
    mean_l2fc_per_gene = np.nanmean(L2FC, axis=1)                # un valore per gene (media sui bin)
    # solo geni con valori finiti e lunghezza > 0
    valid = np.isfinite(mean_l2fc_per_gene) & np.isfinite(lengths) & (lengths > 0)
    rho, pval = spearmanr(np.log10(lengths[valid]), mean_l2fc_per_gene[valid])
    print(f"\ncorrelazione Spearman log10(lunghezza) vs log2FC medio per gene: rho={rho:.3f} (p={pval:.2e}, n={valid.sum()})")
    if abs(rho) > 0.2:
        print("ATTENZIONE: |rho| > 0.2 -- una parte non trascurabile del pattern nella heatmap potrebbe riflettere")
        print("la compressione differenziale dei geni lunghi vs corti nello scaling TSS->TES, non solo biologia.")

    # scatter separato (lunghezza log vs log2FC per gene), file a parte -- non fa parte della griglia sopra
    fig5, ax5 = plt.subplots(figsize=(5.5, 5))
    ax5.scatter(lengths[valid]/1000, mean_l2fc_per_gene[valid], s=4, alpha=0.25, color="#3A6EA5", linewidths=0)
    ax5.set_xscale("log")
    ax5.set_xlabel("lunghezza gene (kb, scala log)"); ax5.set_ylabel(f"log2FC medio ({a.nameA}/{a.nameB})")
    ax5.set_title(f"Spearman rho={rho:.3f} (p={pval:.1e}, n={valid.sum()})", fontsize=10)
    ax5.axhline(0, color="grey", lw=0.7, ls="--")
    ax5.spines[["top","right"]].set_visible(False)
    fig5.tight_layout(); fig5.savefig(f"{a.out}_length_vs_log2fc_scatter.png", dpi=200)
else:
    fig2.colorbar(im, ax=ax2, label="log2FC", shrink=.7)

fig2.savefig(f"{a.out}_log2fc_heatmap.png", dpi=200)         # heatmap (con o senza pannelli lunghezza)

pd.DataFrame(L2FC).to_csv(f"{a.out}_log2fc.tsv", sep="\t", header=False, index=False)   # valori numerici (ordine originale)
print(f"regioni: {n_genes} | log2FC medio globale: {np.nanmean(L2FC):.3f} | mediana: {np.nanmedian(L2FC):.3f}")
outputs = [f"{a.out}_log2fc_profile.png", f"{a.out}_log2fc_heatmap.png", f"{a.out}_log2fc.tsv"]
if a.bed: outputs.append(f"{a.out}_length_vs_log2fc_scatter.png")
print("salvati:", *outputs)
