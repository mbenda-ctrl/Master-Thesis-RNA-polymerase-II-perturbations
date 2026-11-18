#!/usr/bin/env python3
"""
Heatmap integrata GUIDATA DA RNA-seq / RNA nascent (6 colonne: 3 baseline | 3 log2FC), un confronto per figura.

v2 - differenze rispetto alla versione originale:
  * RNA-seq: geni significativi = liste UP/DOWN (meno restrittive: nessun filtro >= 10 conte per campione;
    padj < 0.05, NaCl < 0.1, |log2FC| >= 0.58). log2FC RNA-seq: valore della lista; per i geni non in lista
    (significativi solo nel nascent) log2FC MLE dalla tabella DESeq2 completa (--rna-all-dir).
  * nascent: significativi = colonna "status" (UP/DOWN) dei risultati definitivi (size factor globali sulle
    reads modificate, >= 5 reads in ogni campione del confronto); log2FC = log2FC_MLE.
  * baseline: invariate rispetto all'originale.
  * file DiffBind .csv (KCl) letti con la virgola; confronti senza conteggi ChIP: baseline Pol II vuota.

Geni: unione dei geni SIGNIFICATIVI in RNA-seq o in RNA nascent, tra i geni del BED intersecato (collega x MANE):
  RNA-seq : conte grezze >= --min-rna (10) in tutti e 4 i campioni + padj < --padj e |log2FC| > 0.58 (DESeq2)
  nascent : read modificate >= --min-new (5) in tutti e 4 i campioni + padj < --padj e |log2FC| > 0.58 (DESeq2)
  (--padj-kcl per il solo confronto KCl / NaCl, es. 0.1)
Baseline (media dei 2 controlli):
  Pol II  = log2 CPM (IP - Input), CPM = (IP - Input x read IP / read Input) / read IP x 1e6
  RNA-seq = log2(conte normalizzate DESeq2 + 1)   (size factor median-of-ratios sui campioni dell'esperimento)
  nascent = log2(read modificate normalizzate DESeq2 + 1)   (size factor sui 4 campioni del confronto)
log2FC:
  RNA-seq = liste UP/DOWN (o DESeq2 completo per i geni fuori lista); nascent = log2FC_MLE DESeq2
  Pol II, confronti con DiffBind (--diffbind: Slow vs Amr, KCl vs UT*) = Fold DiffBind (log2, gia' con shrinkage):
      media pesata (peso 2^Conc) dei siti il cui centro cade nel gene +- --fianco (1 kb).
      Gene SIGNIFICATIVO nel ChIP se ha almeno un sito con FDR < 0.05 (v2: senza filtro Conc; --copertura-chip per riattivarlo):
      colori rosso/blu; NON significativo: SCALA DI GRIGI (barra dei colori a parte). Nessun sito: bianco.
  Pol II, altri confronti (Fast vs Amr, Storm vs UT) = log2FC nominale sul CPM (IP - Input):
      media dei log2(CPM trattato / CPM controllo) dei replicati; solo se read nette >= --min-chip (20) in tutti e 4,
      altrimenti bianco.
Clustering: Ward sui 3 log2FC (valori mancanti = 0 solo per le distanze, limite +- --clip-cluster), --k gruppi.
La barra laterale mostra solo i gruppi (niente rosso / blu).

Uso:
  python heatmap_guidata_rna.py --counts-dir input \
      --diffbind Slow_vs_Amr ../chip/Slow_vs_Amr_DiffBind_FULL.tsv --diffbind KCl_vs_UTpub ../chip/dba_kcl_report.csv \
      --rna-lists-dir ../rnaseq_lists --rna-all-dir ../RNAseq --nascent-dir ../nascent \
      --rna-counts ../../expressed_genes_counts_matrix.tsv --nascent-csv-dir ../../nascent/input \
      --outdir heatmap_v2 --confronti Slow_vs_Amr KCl_vs_UTpub --k 4 --clip-cluster 3
"""
import argparse, glob, os, re                                             # argomenti, file, percorsi, espressioni regolari
import numpy as np                                                        # calcolo numerico
import pandas as pd                                                       # tabelle
import matplotlib                                                         # grafici
matplotlib.use("Agg")                                                     # figure solo su file
import matplotlib.pyplot as plt                                           # interfaccia per disegnare
from matplotlib.colors import LinearSegmentedColormap                     # scale di colore
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster         # clustering gerarchico

ap = argparse.ArgumentParser()                                            # lettore degli argomenti
ap.add_argument("--counts-dir", required=True, help="conteggi_gene.tab e librerie.tsv (nostri ChIP)")
ap.add_argument("--counts-dir-kcl", default=None, help="conteggi_gene.tab e librerie.tsv (ChIP KCl / UT*)")
ap.add_argument("--diffbind", nargs=2, action="append", default=[], metavar=("CONFRONTO", "TSV"),
                help="confronti con il Pol II da DiffBind (es. Slow_vs_Amr file.tsv)")
ap.add_argument("--rna-lists-dir", required=True)                         # liste RNA-seq <nome>_UP.csv / _DOWN.csv
ap.add_argument("--rna-all-dir", required=True)                           # tabelle DESeq2 RNA-seq complete (<confronto>_all_genes.csv)
ap.add_argument("--nascent-dir", required=True)                           # risultati nascent (<confronto>_all_genes.csv)
ap.add_argument("--rna-counts", required=True)                            # conte grezze RNA-seq
ap.add_argument("--nascent-csv-dir", required=True)                       # CSV RNAkinet
ap.add_argument("--outdir", required=True)
ap.add_argument("--confronti", nargs="+", default=None, help="solo questi confronti")
ap.add_argument("--min-rna", type=float, default=10); ap.add_argument("--min-new", type=float, default=5)
ap.add_argument("--min-chip", type=float, default=20)                     # read nette minime (ChIP nominale e DiffBind)
ap.add_argument("--padj", type=float, default=0.05); ap.add_argument("--lfc", type=float, default=0.58)
ap.add_argument("--padj-kcl", type=float, default=None, help="padj per il solo confronto KCl / NaCl (default = --padj)")
ap.add_argument("--fdr", type=float, default=0.05, help="FDR DiffBind per il ChIP significativo")
ap.add_argument("--fianco", type=int, default=1000, help="bp aggiunte ai due lati del gene per assegnare i siti DiffBind")
ap.add_argument("--k", type=int, default=4); ap.add_argument("--clip-cluster", type=float, default=3)
ap.add_argument("--lim", type=float, default=2, help="i colori dei log2FC saturano a +-lim")
ap.add_argument("--tieni-senza-nascent", action="store_true",
                help="tiene i geni senza log2FC nascent (cella vuota; 0 nel clustering)")
ap.add_argument("--suffisso", default="", help="aggiunto ai nomi dei file di uscita")
ap.add_argument("--copertura-chip", action="store_true",
                help="richiede anche Conc >= log2(--min-chip) in entrambi i gruppi per il Pol II significativo (come l'originale)")
ap.add_argument("--solo-stech-variata", action="store_true", help="con --solo-m6a: tiene solo i geni con stechiometria variata")
ap.add_argument("--solo-m6a", action="store_true", help="tiene solo i geni con dato m6A e fa il clustering sulle colonne m6A")
ap.add_argument("--m6a-colonne", default="siti,stech",
                help="colonne m6A da mostrare: siti, stech o siti,stech (stechiometria a colori solo se variata, grigia se no)")
ap.add_argument("--m6a", default=None,
                help="Excel con un foglio per confronto (log2FC_n_siti, log2FC_stech): aggiunge 2 colonne m6A "
                     "sulla stessa scala dei log2FC (non usate nel clustering)")
ap.add_argument("--soglie-chip", nargs=3, action="append", default=[], metavar=("CONFRONTO", "FDR", "FOLD"),
                help="soglie ChIP specifiche per confronto: regione significativa se FDR < FDR e |Fold| > FOLD "
                     "(senza il criterio Conc); es. --soglie-chip KCl_vs_UTpub 0.01 1.5")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

# confronto -> (ChIP tratt, ChIP ctrl, colonne Conc DiffBind (tratt, ctrl), integrazione, RNA ctrl, RNA tratt,
#               campioni esperimento RNA, nascent ctrl, nascent tratt, titolo)
INFO = {
    "Fast_vs_Amr": ("Fast", "Amr", ("Fast", "Amr"), "Fast_vs_Amr", "Amr", "E11", ["Amr", "C412", "E11"], "AMR", "E11", "Fast vs Amr"),
    "Slow_vs_Amr": ("Slow", "Amr", ("Slow", "Amr"), "Slow_vs_Amr", "Amr", "C412", ["Amr", "C412", "E11"], "AMR", "C412", "Slow vs Amr"),
    "Storm_vs_UT": ("Storm", "UT", ("Storm", "UT"), "Storm_vs_WT", "WT", "STORM", ["WT", "STORM", "NACL"], "WT", "STORM", "Storm vs UT"),
    "KCl_vs_UTpub": ("KCl", "UTpub", ("KCl", "UT"), "NaCl_vs_WT", "WT", "NACL", ["WT", "STORM", "NACL"], "WT", "NACL",
                     "ChIP KCl vs UT*  |  RNA NaCl vs UT"),
}
DIFFBIND = dict(a.diffbind)                                               # confronti con DiffBind
# confronto -> (prefisso liste RNA-seq, tabella RNA-seq completa, tabella nascent)
SRC = {"Fast_vs_Amr": ("Fast_vs_Amr", "Fast_vs_Amr", "Fast_vs_Amr"),
       "Slow_vs_Amr": ("Slow_vs_Amr", "Slow_vs_Amr", "Slow_vs_Amr"),
       "Storm_vs_UT": ("STORM_vs_WT", "Storm_vs_UT", "Storm_vs_UT"),
       "KCl_vs_UTpub": ("NACL_vs_WT", "NaCl_vs_UT", "NaCl_vs_UT")}

# ---------------- ChIP: conte per gene (nostri + pubblici) e librerie ----------------
CH = pd.read_csv(f"{a.counts_dir}/conteggi_gene.tab", sep="\t", index_col=0)
COORD = CH[["chr", "start", "end"]].copy()                                # coordinate dei geni (per DiffBind)
COORD["chr"] = "chr" + COORD.chr.astype(str).str.replace("^chr", "", regex=True)
CHIP = CH.drop(columns=["chr", "start", "end"])
LIB = pd.read_csv(f"{a.counts_dir}/librerie.tsv", sep="\t", index_col=0)["read_mappate"]
if a.counts_dir_kcl:                                                      # conteggi pubblici KCl / UT*
    K = pd.read_csv(f"{a.counts_dir_kcl}/conteggi_gene.tab", sep="\t", index_col=0).drop(columns=["chr", "start", "end"])
    CHIP = CHIP.join(K[~K.index.duplicated()], how="left")
    LIB = pd.concat([LIB, pd.read_csv(f"{a.counts_dir_kcl}/librerie.tsv", sep="\t", index_col=0)["read_mappate"]])
CHIP = CHIP[~CHIP.index.duplicated()]
U = CHIP.index                                                            # universo: geni del BED intersecato
COORD = COORD[~COORD.index.duplicated()]

# ---------------- RNA-seq e nascent: conte grezze (ENSG) ----------------
raw = pd.read_csv(a.rna_counts, sep="\t", header=None, index_col=0, dtype=str)
hdr = raw.loc["gene_id"].tolist(); raw = raw.drop(index="gene_id"); raw.columns = [x.replace("HEK_", "") for x in hdr]
RNA = raw.astype(float)
RNA.index = RNA.index.str.replace(r"\.\d+$", "", regex=True); RNA = RNA.groupby(level=0).sum()
MOD = {}; NAME = {}
files = sorted(glob.glob(os.path.join(a.nascent_csv_dir, "*rnakinet_unique_clean.csv")))
for f in files:
    s = re.search(r"HEK_(\w+?_R\d)", os.path.basename(f)).group(1); d = pd.read_csv(f)
    d["gid"] = d.gene_id.str.replace(r"\.\d+$", "", regex=True)
    MOD[s] = d.groupby("gid")["n_reads_modificate_5eu"].sum()
    NAME.update(dict(zip(d.gid, d.gene_name)))
MOD = pd.DataFrame(MOD).fillna(0)


def size_factors(C):
    """Size factor DESeq2 (median-of-ratios)."""
    P = C[(C > 0).all(axis=1)]; lg = np.log(P)
    return np.exp(lg.sub(lg.mean(axis=1), axis=0).median(axis=0))


SOGLIE = {k: (float(x), float(y)) for k, x, y in a.soglie_chip}


def diffbind_gene(f, t, c, comp=None):
    """Per gene: Fold DiffBind (media pesata dei siti) e significativita' (FDR + copertura Conc >= log2(min-chip))."""
    S = pd.read_csv(f, sep="," if f.endswith(".csv") else "\t"); S.columns = [x.replace("p-value", "p.value") for x in S.columns]
    S["centro"] = (S.start + S.end) // 2
    if comp in SOGLIE:                                                    # soglie specifiche del confronto
        fdr_c, fold_c = SOGLIE[comp]
        S["sig"] = (S.FDR < fdr_c) & (S.Fold.abs() > fold_c)
    elif a.copertura_chip:                                                # FDR + copertura Conc >= log2(min-chip)
        S["sig"] = (S.FDR < a.fdr) & (S[f"Conc_{t}"] >= np.log2(a.min_chip)) & (S[f"Conc_{c}"] >= np.log2(a.min_chip))
    else:                                                                 # default v2: solo FDR
        S["sig"] = S.FDR < a.fdr
    S["w"] = 2 ** S.Conc                                                  # peso: read normalizzate del sito
    out = {}
    for ch, s in S.groupby("seqnames"):
        g = COORD[COORD.chr == ch]
        if g.empty: continue
        cc = s.centro.values[:, None]
        dentro = (cc >= g.start.values - a.fianco) & (cc <= g.end.values + a.fianco)
        for i, j in zip(*np.nonzero(dentro)):
            out.setdefault(g.index[j], []).append((s.w.values[i], s.Fold.values[i], s.sig.values[i]))
    lfc = pd.Series({k: sum(w * f_ for w, f_, _ in v) / sum(w for w, _, _ in v) for k, v in out.items()}, dtype=float)
    sig = pd.Series({k: any(x for _, _, x in v) for k, v in out.items()}, dtype=bool)
    return lfc, sig


CMAP_BASE = LinearSegmentedColormap.from_list("b", ["#f5efe0", "#e9a07a", "#b2182b"])          # baseline
CMAP_LFC = LinearSegmentedColormap.from_list("l", ["#2166ac", "#9cc3dc", "#f5efe0", "#ef9a7c", "#b2182b"])   # log2FC
CMAP_GRIGIO = LinearSegmentedColormap.from_list("g", ["#f2f2f2", "#9e9e9e", "#303030"])        # ChIP non significativo

riepilogo = []
for comp, (t, c, (ct, cc), integ, rc, rt, rexp, nc, nt, titolo) in INFO.items():
    if a.confronti and comp not in a.confronti: continue
    print(f"\n=== {comp} ===")
    padj = 0.1 if comp == "KCl_vs_UTpub" else 0.05                     # soglie usate nelle liste (solo per il riepilogo)
    lst, rall, nas = SRC[comp]
    nv = lambda x: x.str.replace(r"\.\d+$", "", regex=True)
    up = pd.read_csv(f"{a.rna_lists_dir}/{lst}_UP.csv", index_col=0); dn = pd.read_csv(f"{a.rna_lists_dir}/{lst}_DOWN.csv", index_col=0)
    L = pd.concat([up, dn]); L.index = nv(L.index.to_series())               # liste RNA-seq (significativi)
    RA = pd.read_csv(f"{a.rna_all_dir}/{rall}_all_genes.csv"); RA.index = nv(RA.gene_id)
    NA_ = pd.read_csv(f"{a.nascent_dir}/{nas}_all_genes.csv"); NA_.index = nv(NA_.gene_id)
    I = pd.DataFrame(index=RA.index.union(NA_.index).union(L.index))
    I["gene_id_nover"] = I.index
    I["gene_name"] = [NAME.get(g, s_) for g, s_ in zip(I.index, RA.symbol.reindex(I.index))]
    I["RNAseq_baseMean"] = L.baseMean.combine_first(RA.baseMean).reindex(I.index)
    I["RNAseq_log2FoldChange"] = L.log2FoldChange.combine_first(RA.log2FC_MLE).reindex(I.index)
    I["RNAseq_sig"] = I.index.isin(L.index)
    I["nuovoRNA_log2FoldChange"] = NA_.log2FC_MLE.reindex(I.index)
    I["nuovoRNA_sig"] = NA_.status.reindex(I.index).isin(["UP", "DOWN"]).values
    I = I.dropna(subset=["gene_name"]).sort_values("RNAseq_baseMean", ascending=False)
    I = I.drop_duplicates("gene_name").set_index("gene_name").reindex(U)
    ens = I.gene_id_nover                                                 # ENSG di ogni simbolo
    # --- significativi: liste RNA-seq (nessun filtro di copertura) e status del nascent ---
    sig_r = I.RNAseq_sig.fillna(False).astype(bool)
    sig_n = I.nuovoRNA_sig.fillna(False).astype(bool)
    geni = U[(sig_r | sig_n).values]
    D = pd.DataFrame(index=geni)
    D["log2FC RNA-seq"] = I.RNAseq_log2FoldChange.reindex(geni)
    D["log2FC RNA nascent"] = I.nuovoRNA_log2FoldChange.reindex(geni)
    n0 = len(D)
    if a.tieni_senza_nascent:                                             # basta il log2FC RNA-seq
        D = D.dropna(subset=["log2FC RNA-seq"])
    else:                                                                 # servono entrambi i log2FC DESeq2
        D = D.dropna()
    print(f"  significativi (padj < {padj:g}): RNA-seq {int(sig_r.sum())}, nascent {int(sig_n.sum())}, unione {n0}"
          f" -> in heatmap: {len(D)} (senza log2FC nascent: {int(D['log2FC RNA nascent'].isna().sum())})")
    # --- baseline RNA-seq e nascent: conte normalizzate DESeq2 ---
    CR = RNA[[x for x in RNA.columns if x.rsplit("_", 1)[0] in rexp]]   # esperimento RNA-seq
    NR = CR / size_factors(CR)
    D["baseline RNA-seq"] = np.log2(NR.reindex(ens.reindex(D.index).values)[[f"{rc}_r1", f"{rc}_r2"]].mean(axis=1).values + 1)
    CN = MOD[[f"{nc}_R1", f"{nc}_R2", f"{nt}_R1", f"{nt}_R2"]]
    NN = CN / size_factors(CN)
    D["baseline RNA nascent"] = np.log2(NN.reindex(ens.reindex(D.index).values)[[f"{nc}_R1", f"{nc}_R2"]].mean(axis=1).values + 1)
    # --- Pol II: CPM (IP - Input) ---
    reps = [r for r in (1, 2) if f"{c}_r{r}" in CHIP.columns and f"{t}_r{r}" in CHIP.columns]
    s4 = [f"{c}_r{r}" for r in reps] + [f"{t}_r{r}" for r in reps]; nr = len(reps)
    if nr:
        NET = pd.DataFrame({s: CHIP[s].astype(float) - CHIP[f"{s}_input"].astype(float) * LIB[s] / LIB[f"{s}_input"] for s in s4})
        CPM = NET / LIB[s4].values * 1e6
        base_pol = CPM.loc[D.index, s4[:nr]].mean(axis=1)
        D["baseline Pol II"] = np.log2(base_pol.where(base_pol > 0))      # IP <= Input: baseline non definita
    else:
        print(f"  ATTENZIONE: conteggi ChIP di {c}/{t} non disponibili -> baseline Pol II vuota")
        D["baseline Pol II"] = np.nan
    sig_chip = None
    if comp in DIFFBIND:                                                  # Pol II da DiffBind
        lfc_db, sig_db = diffbind_gene(DIFFBIND[comp], ct, cc, comp)
        D["log2FC Pol II"] = lfc_db.reindex(D.index)
        sig_chip = sig_db.reindex(D.index).fillna(False).astype(bool)
        print(f"  Pol II DiffBind: {int(D['log2FC Pol II'].notna().sum())} geni con siti, {int(sig_chip.sum())} significativi")
    else:                                                                 # log2FC nominale sul CPM
        ok = (NET.loc[D.index, s4] >= a.min_chip).all(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            lf = sum(np.log2(CPM.loc[D.index, s4[nr + i]] / CPM.loc[D.index, s4[i]]) for i in range(nr)) / nr
        D["log2FC Pol II"] = lf.where(ok)
        print(f"  Pol II nominale: {int(ok.sum())} geni con read nette >= {a.min_chip:g} in tutti i campioni")
    # --- m6A (opzionale): log2FC numero di siti e stechiometria ---
    M6A_SHEET = {"Fast_vs_Amr": "Fast_vs_Amr", "Slow_vs_Amr": "Slow_vs_Amr",
                 "Storm_vs_UT": "Storm_vs_WT", "KCl_vs_UTpub": "NaCl_vs_WT"}
    extra_cols = []
    if a.m6a:
        Mx = pd.read_excel(a.m6a, sheet_name=M6A_SHEET[comp])
        Mx.index = Mx.gene_id.str.replace(r"\.\d+$", "", regex=True)
        Mx = Mx[~Mx.index.duplicated()]
        e = ens.reindex(D.index).values
        D["m6A n. siti"] = Mx.log2FC_n_siti.reindex(e).values
        D["m6A stechiometria"] = Mx.log2FC_stech.reindex(e).values
        D["sig m6A stechiometria"] = Mx.variazione_stechiometria.reindex(e).isin(["aumenta", "diminuisce"]).values
        scelta = a.m6a_colonne.split(",")
        extra_cols = [c for k, c in (("siti", "m6A n. siti"), ("stech", "m6A stechiometria")) if k in scelta]
        print(f"  m6A: {int(D['m6A n. siti'].notna().sum())} geni della heatmap con variazione m6A")
    # --- clustering (sui 3 log2FC Pol II / RNA-seq / nascent, oppure sulle colonne m6A con --solo-m6a) ---
    lfc_cols = ["log2FC Pol II", "log2FC RNA-seq", "log2FC RNA nascent"]
    if a.solo_m6a:
        D = D[D["m6A stechiometria"].notna() & (D["sig m6A stechiometria"] | (not a.solo_stech_variata))].copy()
        print(f"  solo geni con m6A: {len(D)}")
    clu_cols = extra_cols if a.solo_m6a else lfc_cols
    X = D[clu_cols].fillna(0).values.clip(-a.clip_cluster, a.clip_cluster)
    Z = linkage(X, method="ward"); D["gruppo"] = fcluster(Z, a.k, criterion="maxclust")
    ordine = dendrogram(Z, no_plot=True)["leaves"]
    D = D.iloc[ordine]; D["gruppo"] = D.gruppo.map(dict(zip(pd.unique(D.gruppo), range(1, a.k + 1))))
    if sig_chip is not None: sig_chip = sig_chip.reindex(D.index)
    D.to_csv(os.path.join(a.outdir, f"heatmap_rna_{comp}{a.suffisso}.tsv"), sep="\t", index_label="gene")
    # --- figura ---
    n = len(D)
    ne = len(extra_cols)
    wr = [1.2, 0.22, 1, 1, 1, 0.35, 1, 1, 1] + ([0.35] + [1] * ne if ne else [])
    fig = plt.figure(figsize=(10.5 + 2.3 * ne, 10))
    gs = fig.add_gridspec(2, len(wr), width_ratios=wr, height_ratios=[22, 0.5], wspace=0.06, hspace=0.12)
    axd = fig.add_subplot(gs[0, 0])                                       # dendrogramma
    dendrogram(Z, orientation="left", ax=axd, no_labels=True, color_threshold=0,
               above_threshold_color="#444444", link_color_func=lambda k: "#444444")
    axd.set_ylim(10 * n, 0); axd.axis("off")
    axs = fig.add_subplot(gs[0, 1])                                       # barra dei gruppi (grigi alternati, senza rosso/blu)
    axs.imshow((D.gruppo.values % 2)[:, None], aspect="auto", interpolation="nearest",
               cmap=LinearSegmentedColormap.from_list("gr", ["#bdbdbd", "#7f7f7f"]), vmin=0, vmax=1)
    axs.set_xticks([]); axs.tick_params(left=False, labelleft=False)
    bounds = np.flatnonzero(np.diff(D.gruppo.values)) + 0.5
    for b0, b1 in zip(np.r_[-0.5, bounds], np.r_[bounds, n - 0.5]):
        axs.text(0, (b0 + b1) / 2, str(D.gruppo.values[int(b0 + 0.5)]), ha="center", va="center", fontsize=10, fontweight="bold", color="white")
    for j, col in enumerate(["Pol II", "RNA-seq", "RNA nascent"]):        # baseline
        ax = fig.add_subplot(gs[0, 2 + j], sharey=axs)
        v = D[f"baseline {col}"].values.astype(float)
        lo, hi = np.nanpercentile(v, [2, 98]) if np.isfinite(v).any() else (0, 1)
        im = ax.imshow(np.ma.masked_invalid(v)[:, None], aspect="auto", cmap=CMAP_BASE, vmin=lo, vmax=hi, interpolation="nearest")
        ax.set_xticks([0]); ax.set_xticklabels([col], fontsize=9); ax.xaxis.tick_top()
        ax.tick_params(left=False, labelleft=False, bottom=False, top=False)
        for b in bounds: ax.axhline(b, color="k", lw=0.9)
        cb = fig.colorbar(im, cax=fig.add_subplot(gs[1, 2 + j]), orientation="horizontal"); cb.ax.tick_params(labelsize=6)
        cb.set_label("log2 CPM\n(IP - Input)" if col == "Pol II" else "log2 conte\nnorm. DESeq2", fontsize=7)   # etichette su 2 righe
        if j == 1: ax.set_title(f"baseline ({'UT*' if c == 'UTpub' else c}{', RNA UT' if c == 'UTpub' else ''})", fontsize=10, pad=20)
    for j, col in enumerate(lfc_cols):                                    # log2FC
        ax = fig.add_subplot(gs[0, 6 + j], sharey=axs)
        v = D[col].values.astype(float)
        if col == "log2FC Pol II" and sig_chip is not None:               # DiffBind: significativi a colori, gli altri in grigio
            vs = np.where(sig_chip.values, v, np.nan); vn = np.where(~sig_chip.values, v, np.nan)
            im = ax.imshow(np.ma.masked_invalid(vs)[:, None], aspect="auto", cmap=CMAP_LFC, vmin=-a.lim, vmax=a.lim, interpolation="nearest")
            img = ax.imshow(np.ma.masked_invalid(vn)[:, None], aspect="auto", cmap=CMAP_GRIGIO, vmin=-a.lim, vmax=a.lim, interpolation="nearest")
            if np.isfinite(vn).any():                                     # barra del grigio solo se serve
                cbg = fig.colorbar(img, cax=fig.add_axes([0.93, 0.32, 0.012, 0.25]), extend="both")
                cbg.set_label("Pol II DiffBind non significativo (log2FC)", fontsize=7); cbg.ax.tick_params(labelsize=6)
        else:
            im = ax.imshow(np.ma.masked_invalid(v)[:, None], aspect="auto", cmap=CMAP_LFC, vmin=-a.lim, vmax=a.lim, interpolation="nearest")
        ax.set_xticks([0]); ax.set_xticklabels([col.replace("log2FC ", "")], fontsize=9); ax.xaxis.tick_top()
        ax.tick_params(left=False, labelleft=False, bottom=False, top=False)
        for b in bounds: ax.axhline(b, color="k", lw=0.9)
        if j == 1:
            fc = f"{t} / {'UT*' if c == 'UTpub' else c}" + (", RNA NaCl / UT" if c == "UTpub" else "")
            ax.set_title(f"log2FC ({fc})", fontsize=10, pad=20)
    for j, col in enumerate(extra_cols):                                  # m6A: stessa scala dei log2FC
        ax = fig.add_subplot(gs[0, 10 + j], sharey=axs)
        v = D[col].values.astype(float)
        if col == "m6A stechiometria":                                    # variata a colori, non variata in grigio
            sg = D["sig m6A stechiometria"].values.astype(bool)
            vs = np.where(sg, v, np.nan); vn = np.where(~sg, v, np.nan)
            im = ax.imshow(np.ma.masked_invalid(vs)[:, None], aspect="auto", cmap=CMAP_LFC, vmin=-a.lim, vmax=a.lim, interpolation="nearest")
            ax.imshow(np.ma.masked_invalid(vn)[:, None], aspect="auto", cmap=CMAP_GRIGIO, vmin=-a.lim, vmax=a.lim, interpolation="nearest")
        else:
            im = ax.imshow(np.ma.masked_invalid(v)[:, None], aspect="auto", cmap=CMAP_LFC, vmin=-a.lim, vmax=a.lim, interpolation="nearest")
        ax.set_xticks([0]); ax.set_xticklabels([col.replace("m6A ", "")], fontsize=9); ax.xaxis.tick_top()
        ax.tick_params(left=False, labelleft=False, bottom=False, top=False)
        for b in bounds: ax.axhline(b, color="k", lw=0.9)
        if j == 0:
            ax.set_title("m6A (log2FC)", fontsize=10, pad=20, loc="left")
    cb = fig.colorbar(im, cax=fig.add_subplot(gs[1, 6:(9 if not ne else 10 + ne)]), orientation="horizontal", extend="both")
    cb.set_label(("log2FC (Pol II DiffBind, RNA DESeq2" if sig_chip is not None else "log2FC (Pol II nominale su CPM, RNA DESeq2")
                 + ((", m6A " + " e ".join(c.replace("m6A ", "") for c in extra_cols) + ")") if ne else ")") + f" - colori saturati a ±{a.lim:g}; bianco = dato assente" + ("\nstechiometria m6A in grigio = non variata (FDR >= 0.05 o |log2FC| < 0.58)" if "m6A stechiometria" in extra_cols else ""),
                 fontsize=8); cb.ax.tick_params(labelsize=7)
    fig.suptitle(titolo + (" - geni DE RNA-seq / nascent con variazione m6A (gruppi sulle colonne m6A)" if a.solo_m6a
                           else " - geni DE RNA-seq / nascent"), fontsize=12, y=1.0)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(a.outdir, f"heatmap_rna_{comp}{a.suffisso}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    M = D.groupby("gruppo")[lfc_cols + extra_cols].median().round(2).assign(n=D.groupby("gruppo").size())
    print(M.to_string())
    riepilogo.append({"confronto": comp, "padj": padj, "sig_RNAseq": int(sig_r.sum()), "sig_nascent": int(sig_n.sum()),
                      "unione": n0, "in_heatmap": len(D), "senza_nascent": int(D["log2FC RNA nascent"].isna().sum()), "PolII_con_valore": int(D["log2FC Pol II"].notna().sum()),
                      "PolII_significativi": int(sig_chip.sum()) if sig_chip is not None else np.nan})
pd.DataFrame(riepilogo).to_csv(os.path.join(a.outdir, f"riepilogo_heatmap_rna{a.suffisso}.tsv"), sep="\t", index=False)
print("\n" + pd.DataFrame(riepilogo).to_string(index=False))
