#!/usr/bin/env python3
"""
PASSO 8a - Geni GENCODE (o altro GTF) con almeno un picco dentro l'unita trascrizionale.

Unita trascrizionale = corpo del gene (feature 'gene' del GTF) piu, opzionalmente, una finestra a monte
del TSS (--up) e a valle del TES (--down), tenendo conto dello strand. Un picco "cade" nel gene se si
sovrappone anche di 1 bp alla finestra. Nessuna assegnazione al gene piu vicino: un picco puo contare
per piu geni (es. promoter bidirezionali, geni annidati).

Input geni: GTF/GTF.gz (usa le righe 'gene') oppure TSV senza header con colonne
            chrom, start0, end, strand, gene_id, gene_name, gene_type
Input picchi: uno o piu file, BED/narrowPeak (start 0-based) oppure TSV annotati con header
              seqnames/start/end (start 1-based, come quelli di ChIPseeker).

Uso:
  python peaks_to_genes.py --genes gencode.vNN.annotation.gtf.gz --peaks picchi/*.tsv --out risultati \
         [--up 3000 --down 300] [--biotype protein_coding,lncRNA] [--min-len 0] \
         [--core-exclude Nacl_1,C412_1,C422_2]
"""
import argparse, gzip, os, re, sys
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--genes", required=True); ap.add_argument("--peaks", nargs="+", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--up", type=int, default=3000, help="bp a monte del TSS inclusi nell'unita (default 3000)")
ap.add_argument("--down", type=int, default=300, help="bp a valle del TES inclusi (default 300)")
ap.add_argument("--biotype", default=None, help="lista separata da virgole, es. protein_coding,lncRNA")
ap.add_argument("--min-len", type=int, default=0)
ap.add_argument("--min-samples", type=int, default=2, help="set consensus: geni con picchi in almeno N campioni (default 2)")
ap.add_argument("--core-exclude", default="", help="sottostringhe dei campioni da escludere per l'intersezione 'core'")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)

# ---------------- geni ----------------
def read_genes(path):
    if re.search(r"\.gtf(\.gz)?$", path):
        op = gzip.open if path.endswith(".gz") else open
        rows = []
        get = lambda s, k: (m.group(1) if (m := re.search(rf'{k} "([^"]+)"', s)) else None)
        with op(path, "rt") as fh:
            for l in fh:
                if l[0] == "#": continue
                f = l.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "gene": continue
                rows.append([f[0], int(f[3]) - 1, int(f[4]), f[6], get(f[8], "gene_id"),
                             get(f[8], "gene_name"), get(f[8], "gene_type") or get(f[8], "gene_biotype")])
        return pd.DataFrame(rows, columns=["chrom", "start0", "end", "strand", "gene_id", "gene_name", "gene_type"])
    return pd.read_csv(path, sep="\t", header=None,
                       names=["chrom", "start0", "end", "strand", "gene_id", "gene_name", "gene_type"])

G = read_genes(a.genes)
n_all = len(G)
if a.biotype: G = G[G.gene_type.isin(a.biotype.split(","))]
G = G[(G.end - G.start0) >= a.min_len].drop_duplicates("gene_id").reset_index(drop=True)
plus = G.strand == "+"
G["win_start"] = np.where(plus, G.start0 - a.up, G.start0 - a.down).clip(min=0)
G["win_end"] = np.where(plus, G.end + a.down, G.end + a.up)
G["TSS"] = np.where(plus, G.start0 + 1, G.end); G["TES"] = np.where(plus, G.end, G.start0 + 1)

# ---------------- picchi ----------------
def read_peaks(path):
    with open(path) as fh: first = fh.readline()
    if first.startswith("seqnames"):
        d = pd.read_csv(path, sep="\t", usecols=["seqnames", "start", "end"])
        return pd.DataFrame({"chrom": d.seqnames, "start": d.start - 1, "end": d.end})
    d = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"], comment="#")
    return d

def sname(p):
    b = os.path.basename(p); b = re.sub(r"^\d{8,}_", "", b)
    return re.sub(r"(_narrow_peaks_annotation\.tsv|_peaks\.narrowPeak|\.narrowPeak|\.bed|\.tsv)$", "", b)

samples, counts, summ = [], {}, []
for p in sorted(a.peaks):
    s = sname(p); P = read_peaks(p); samples.append(s)
    c = np.zeros(len(G), dtype=int); in_win = 0
    for chrom, gi in G.groupby("chrom").groups.items():
        Pc = P[P.chrom == chrom]
        if Pc.empty: continue
        ps, pe = np.sort(Pc.start.values), np.sort(Pc.end.values)
        ws, we = G.loc[gi, "win_start"].values, G.loc[gi, "win_end"].values
        # n. picchi con start < we e end > ws  (vale per intervalli qualsiasi, array ordinati separatamente)
        c[G.index.get_indexer(gi)] = np.searchsorted(ps, we, "left") - np.searchsorted(pe, ws, "right")
    counts[s] = c
    # picchi che cadono in almeno una finestra (unione delle finestre)
    for chrom, Pc in P.groupby("chrom"):
        W = G[G.chrom == chrom].sort_values("win_start")
        if W.empty: continue
        ms, me = [], []
        for st, en in zip(W.win_start.values, W.win_end.values):
            if ms and st <= me[-1]: me[-1] = max(me[-1], en)
            else: ms.append(st); me.append(en)
        ms, me = np.array(ms), np.array(me)
        idx = np.searchsorted(ms, Pc.end.values, "left") - 1
        ok = idx >= 0
        in_win += int((me[idx[ok]] > Pc.start.values[ok]).sum())
    summ.append(dict(sample=s, n_peaks=len(P), peaks_in_gene_windows=in_win,
                     pct_in_gene_windows=round(100 * in_win / len(P), 1),
                     genes_with_peak=int((c > 0).sum())))

C = pd.DataFrame(counts)
pres = (C > 0).astype(int)
excl = [x for x in a.core_exclude.split(",") if x]
core = [s for s in samples if not any(x in s for x in excl)]
T = pd.concat([G[["gene_id", "gene_name", "gene_type", "chrom", "start0", "end", "strand", "TSS", "TES", "win_start", "win_end"]], C], axis=1)
T["n_samples"] = pres.sum(axis=1).values
if excl: T["n_core"] = pres[core].sum(axis=1).values
T.to_csv(f"{a.out}/gene_peak_counts.tsv", sep="\t", index=False)
pd.DataFrame(summ).to_csv(f"{a.out}/summary_per_sample.tsv", sep="\t", index=False)

def save(name, mask):
    D = T[mask].sort_values(["chrom", "start0"])
    D.to_csv(f"{a.out}/gene_set_{name}.tsv", sep="\t", index=False)
    pd.DataFrame({"c": D.chrom, "s": D.start0, "e": D.end, "n": D.gene_name.fillna(D.gene_id), "sc": 0, "st": D.strand}
                 ).to_csv(f"{a.out}/gene_set_{name}.bed", sep="\t", header=False, index=False)
    print(f"{name:22s} {mask.sum():6d} geni")

print(f"geni nel GTF: {n_all} | dopo filtri: {len(G)} | campioni: {len(samples)}")
save("union", T.n_samples >= 1)
save(f"union_n{a.min_samples}", T.n_samples >= a.min_samples)        # set consensus (input del passo 08b)
save("intersection", T.n_samples == len(samples))
if excl: save("intersection_core", T.n_core == len(core))
print(pd.DataFrame(summ).to_string(index=False))
