#!/usr/bin/env python3
"""
PASSO 13 - Geni con occupazione differenziale di Pol II (DiffBind) per i metagene differenziali.

Un gene del set consensus (passo 8) e' "differenziale" se almeno un sito DiffBind significativo
(FDR < --fdr e |Fold| > --min-fold) si sovrappone alla sua unita' trascrizionale
(TSS - 3 kb ... TES, rispettando il filamento).
Per ogni confronto scrive il BED dei geni con siti in aumento (UP), in calo (DOWN) e tutti insieme,
con i nomi attesi da metagene/04_run_differential_metagene.sh: genes_diff_<Trattato>vs<Controllo>[_UP|_DOWN].bed

Uso: python 13_diffbind_gene_sets.py --bed genes_pc_union_MANE_n2_min1kb.bed \
        --diffbind-dir results/chipseq/diffbind --outdir results/chipseq/diffbind_gene_sets \
        [--contrasts Slow_vs_Amr Fast_vs_Amr] [--fdr 0.05 --min-fold 1 --upstream 3000]
"""
import argparse, os
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--bed", required=True, help="BED6 dei geni consensus (nome in colonna 4, filamento in colonna 6)")
ap.add_argument("--diffbind-dir", required=True, help="cartella con <confronto>_DiffBind_all_sites.tsv (passo 10)")
ap.add_argument("--outdir", required=True)
ap.add_argument("--contrasts", nargs="+", default=["Slow_vs_Amr", "Fast_vs_Amr"])
ap.add_argument("--fdr", type=float, default=0.05)
ap.add_argument("--min-fold", type=float, default=1.0, help="|Fold| (log2) minimo")
ap.add_argument("--upstream", type=int, default=3000, help="bp a monte del TSS incluse nell'unita' trascrizionale")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

G = pd.read_csv(a.bed, sep="\t", header=None).iloc[:, :6]
G.columns = ["chrom", "start", "end", "name", "score", "strand"]
G["tu_start"] = np.where(G.strand == "+", G.start - a.upstream, G.start).clip(min=0)
G["tu_end"] = np.where(G.strand == "+", G.end, G.end + a.upstream)

def genes_hit(regions):
    """nomi dei geni la cui unita' trascrizionale si sovrappone ad almeno una regione"""
    hit = set()
    for chrom, g in G.groupby("chrom"):
        r = regions[regions.seqnames == chrom]
        if r.empty:
            continue
        gs = g.tu_start.to_numpy(); ge = g.tu_end.to_numpy(); names = g.name.to_numpy()
        for s, e in zip(r.start.to_numpy(), r.end.to_numpy()):
            hit.update(names[(ge > s) & (gs < e)])
    return hit

rows = []
for c in a.contrasts:
    f = os.path.join(a.diffbind_dir, f"{c}_DiffBind_all_sites.tsv")
    if not os.path.exists(f):
        print(f"ATTENZIONE: manca {f}, salto"); continue
    D = pd.read_csv(f, sep="\t")
    D["seqnames"] = D.seqnames.astype(str)
    sig = D[(D.FDR < a.fdr) & (D.Fold.abs() > a.min_fold)]
    trt, ctl = c.split("_vs_")
    tag = f"genes_diff_{trt}vs{ctl}"
    sets = {"": genes_hit(sig), "_UP": genes_hit(sig[sig.Fold > 0]), "_DOWN": genes_hit(sig[sig.Fold < 0])}
    for suf, names in sets.items():
        G[G.name.isin(names)][["chrom", "start", "end", "name", "score", "strand"]].to_csv(
            os.path.join(a.outdir, f"{tag}{suf}.bed"), sep="\t", header=False, index=False)
    rows.append(dict(confronto=c, siti_testati=len(D), siti_significativi=len(sig),
                     geni=len(sets[""]), geni_UP=len(sets["_UP"]), geni_DOWN=len(sets["_DOWN"])))
S = pd.DataFrame(rows)
S.to_csv(os.path.join(a.outdir, "riepilogo_geni_diffbind.tsv"), sep="\t", index=False)
print(S.to_string(index=False))
