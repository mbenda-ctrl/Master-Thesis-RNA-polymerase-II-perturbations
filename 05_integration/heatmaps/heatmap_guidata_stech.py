#!/usr/bin/env python3
"""
Heatmap GUIDATA DALLA STECHIOMETRIA m6A: tutti i geni con stechiometria variata
(variazione_stechiometria = aumenta/diminuisce: |log2FC_stech| > 0.58 e FDR < 0.05), senza altri filtri.

Colonne (log2FC, stessa scala blu-rosso):
  m6A stechiometria | Pol II (DiffBind) | RNA-seq | RNA nascent
  - Pol II: media pesata dei Fold DiffBind delle regioni nel gene +-1 kb (solo geni del BED);
            a colori se significativo (Slow FDR<0.05; KCl FDR<0.01 e |Fold|>1.5; Fast/Storm tutti i siti),
            grigio se non significativo, bianco se il gene non è nel BED o non ha regioni.
  - RNA-seq: liste UP/DOWN a colori, altrimenti log2FC MLE in grigio; bianco se non testato.
  - nascent: log2FC MLE, a colori se status UP/DOWN, grigio altrimenti; bianco se non testato.
Clustering: Ward sulle 4 colonne (vuoti = 0, valori tagliati a +-3), k gruppi.

Uso: python heatmap_guidata_stech.py <cartella_PolII_tutti> <cartella_PolII_sig> <outdir> [k]
"""
import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

POL_ALL, POL_SIG, OUT = sys.argv[1:4]
K = int(sys.argv[4]) if len(sys.argv) > 4 else 4
W_FIG = float(sys.argv[5]) if len(sys.argv) > 5 else None      # larghezza figura (pollici), opzionale
H_FIG = float(sys.argv[6]) if len(sys.argv) > 6 else None      # altezza figura (pollici), opzionale
SOLO = sys.argv[7].split(",") if len(sys.argv) > 7 else None   # solo questi confronti
LIM, CLIP = 2, 3
os.makedirs(OUT, exist_ok=True)
XL = "../m6a/geni_variazione_siti_e_stechiometria_log2FC.xlsx"
# confronto m6A: (heatmap ChIP, liste RNA-seq, tabelle RNA-seq/nascent, Pol II sig da file (None = tutti), titolo)
CONF = {"Slow_vs_Amr": ("Slow_vs_Amr", "Slow_vs_Amr", "Slow_vs_Amr", True, "Slow vs Amr", "Pol II Slow vs Amr (FDR < 0.05)"),
        "Fast_vs_Amr": ("Fast_vs_Amr", "Fast_vs_Amr", "Fast_vs_Amr", False, "Fast vs Amr", "Pol II Fast vs Amr (tutti i siti)"),
        "Storm_vs_WT": ("Storm_vs_UT", "STORM_vs_WT", "Storm_vs_UT", False, "Storm vs UT", "Pol II Storm vs UT (tutti i siti)"),
        "NaCl_vs_WT":  ("KCl_vs_UTpub", "NACL_vs_WT", "NaCl_vs_UT", True, "NaCl vs UT", "Pol II KCl vs UT* (FDR < 0.01, |Fold| > 1.5)")}
CMAP_LFC = LinearSegmentedColormap.from_list("l", ["#2166ac", "#9cc3dc", "#f5efe0", "#ef9a7c", "#b2182b"])
CMAP_GRIGIO = LinearSegmentedColormap.from_list("g", ["#f2f2f2", "#9e9e9e", "#303030"])
nv = lambda s: pd.Series(s).astype(str).str.replace(r"\.\d+$", "", regex=True).values
riep = []
for sh, (hc, lst, tab, polsig, titolo, pol_lab) in CONF.items():
    if SOLO and sh not in SOLO: continue
    M = pd.read_excel(XL, sheet_name=sh); M.index = nv(M.gene_id); M = M[~M.index.duplicated()]
    M = M[M.variazione_stechiometria.isin(["aumenta", "diminuisce"])]
    D = pd.DataFrame(index=M.index); D["gene"] = M.gene_name.values
    D["m6A stechiometria"] = M.log2FC_stech
    # Pol II (per simbolo dal file ChIP; tutti i geni del BED con regioni)
    P = pd.read_csv(f"{POL_ALL}/heatmap_chip_{hc}.tsv", sep="\t", index_col=0)
    D["Pol II"] = P["log2FC Pol II"].reindex(D.gene).values
    if polsig:
        ps = set(pd.read_csv(f"{POL_SIG}/heatmap_chip_{hc}.tsv", sep="\t", index_col=0).index)
        D["sig Pol II"] = D.gene.isin(ps).values
    else:
        D["sig Pol II"] = D["Pol II"].notna()
    # RNA-seq
    L = pd.concat([pd.read_csv(f"../rnaseq_lists/{lst}_{d}.csv", index_col=0) for d in ("UP", "DOWN")]); L.index = nv(L.index)
    RA = pd.read_csv(f"../RNAseq/{tab}_all_genes.csv"); RA.index = nv(RA.gene_id)
    D["RNA-seq"] = L.log2FoldChange.reindex(D.index).combine_first(RA.log2FC_MLE.reindex(D.index))
    D["sig RNA-seq"] = D.index.isin(L.index)
    # nascent
    NA_ = pd.read_csv(f"../nascent/{tab}_all_genes.csv"); NA_.index = nv(NA_.gene_id); NA_ = NA_[~NA_.index.duplicated()]
    D["RNA nascent"] = NA_.log2FC_MLE.reindex(D.index)
    D["sig RNA nascent"] = NA_.status.reindex(D.index).isin(["UP", "DOWN"]).values
    D["sig m6A stechiometria"] = True
    cols = ["m6A stechiometria", "Pol II", "RNA-seq", "RNA nascent"]
    # clustering
    k = min(K, len(D))
    X = D[cols].fillna(0).values.clip(-CLIP, CLIP)
    Z = linkage(X, method="ward"); D["gruppo"] = fcluster(Z, k, criterion="maxclust")
    D = D.iloc[dendrogram(Z, no_plot=True)["leaves"]]
    D["gruppo"] = D.gruppo.map(dict(zip(pd.unique(D.gruppo), range(1, k + 1))))
    D.to_csv(f"{OUT}/heatmap_stech_{sh}.tsv", sep="\t", index_label="gene_id")
    # figura
    n = len(D); nomi = n <= 250
    h = H_FIG or (max(6, min(0.11 * n + 2.5, 30)) if nomi else 10)
    fig = plt.figure(figsize=(W_FIG or (8.8 + (1.6 if nomi else 0)), h))
    wr = [1.1, 0.2, 1, 0.3, 1, 1, 1] + ([1.6] if nomi else [])
    gs = fig.add_gridspec(2, len(wr), width_ratios=wr, height_ratios=[22, 0.35], wspace=0.06, hspace=0.6 / h * 10 * 0.12)
    axd = fig.add_subplot(gs[0, 0])
    dendrogram(Z, orientation="left", ax=axd, no_labels=True, color_threshold=0,
               above_threshold_color="#444444", link_color_func=lambda _: "#444444")
    axd.set_ylim(10 * n, 0); axd.axis("off")
    axs = fig.add_subplot(gs[0, 1])
    axs.imshow((D.gruppo.values % 2)[:, None], aspect="auto", interpolation="nearest",
               cmap=LinearSegmentedColormap.from_list("gr", ["#bdbdbd", "#7f7f7f"]), vmin=0, vmax=1)
    axs.set_xticks([]); axs.tick_params(left=False, labelleft=False)
    bounds = np.flatnonzero(np.diff(D.gruppo.values)) + 0.5
    for b0, b1 in zip(np.r_[-0.5, bounds], np.r_[bounds, n - 0.5]):
        axs.text(0, (b0 + b1) / 2, str(D.gruppo.values[int(b0 + 0.5)]), ha="center", va="center", fontsize=10, fontweight="bold", color="white")
    for j, col in enumerate(cols):
        ax = fig.add_subplot(gs[0, 2 + j + (1 if j else 0)], sharey=axs)
        v = D[col].values.astype(float); sg = D[f"sig {col}"].values.astype(bool)
        im = ax.imshow(np.ma.masked_invalid(np.where(sg, v, np.nan))[:, None], aspect="auto", cmap=CMAP_LFC, vmin=-LIM, vmax=LIM, interpolation="nearest")
        img = ax.imshow(np.ma.masked_invalid(np.where(~sg, v, np.nan))[:, None], aspect="auto", cmap=CMAP_GRIGIO, vmin=-LIM, vmax=LIM, interpolation="nearest")
        ax.set_xticks([0]); ax.set_xticklabels([col.replace("m6A ", "m6A\n")], fontsize=9); ax.xaxis.tick_top()
        ax.tick_params(left=False, labelleft=False, bottom=False, top=False)
        for b in bounds: ax.axhline(b, color="k", lw=0.9)
        if j == 2: ax.set_title(f"log2FC ({titolo})", fontsize=10, pad=22)
        if nomi and j == len(cols) - 1:
            ax.yaxis.tick_right(); ax.tick_params(right=False, labelright=True)
            ax.set_yticks(range(n)); ax.set_yticklabels(D.gene.values, fontsize=max(2.5, min(7, 0.85 * 72 * (h - 2.2) / n)))
    cb = fig.colorbar(im, cax=fig.add_subplot(gs[1, 2:7]), orientation="horizontal", extend="both")
    cb.set_label(f"log2FC - colori saturati a ±{LIM}; grigio = non significativo; bianco = dato assente "
                 f"(Pol II: gene fuori dal BED o senza regioni)\n{pol_lab}", fontsize=7); cb.ax.tick_params(labelsize=7)
    fig.subplots_adjust(top=1 - 1.3 / h, bottom=0.9 / h)
    fig.suptitle(f"{titolo} - tutti i geni con stechiometria m6A variata (n = {n})", fontsize=12, y=1 - 0.35 / h)
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/heatmap_stech_{sh}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    g = D.groupby("gruppo")
    S = pd.DataFrame({"n": g.size(), "stech": g["m6A stechiometria"].median(), "PolII": g["Pol II"].median(),
                      "RNAseq": g["RNA-seq"].median(), "nascent": g["RNA nascent"].median(),
                      "RNAseq_sig": g["sig RNA-seq"].sum(), "nascent_sig": g["sig RNA nascent"].sum(),
                      "PolII_presente": g["Pol II"].count(), "PolII_sig": g["sig Pol II"].sum()}).round(2)
    print(f"\n=== {sh}: {n} geni | nel BED con Pol II {int(D['Pol II'].notna().sum())} | RNA-seq sig {int(D['sig RNA-seq'].sum())}"
          f" | nascent sig {int(D['sig RNA nascent'].sum())} | RNA-seq o nascent sig {int((D['sig RNA-seq'] | D['sig RNA nascent']).sum())}")
    print(S.to_string())
    riep.append({"confronto": sh, "stech_variata": n, "con_PolII": int(D["Pol II"].notna().sum()), "PolII_sig": int(D["sig Pol II"].sum()),
                 "RNAseq_sig": int(D["sig RNA-seq"].sum()), "nascent_sig": int(D["sig RNA nascent"].sum()),
                 "RNAseq_o_nascent_sig": int((D["sig RNA-seq"] | D["sig RNA nascent"]).sum())})
pd.DataFrame(riep).to_csv(f"{OUT}/riepilogo.tsv", sep="\t", index=False)
