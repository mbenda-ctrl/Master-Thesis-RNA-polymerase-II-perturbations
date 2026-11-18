#!/usr/bin/env python3
"""
PASSO 12a - Geni protein-coding GENCODE senza sovrapposizioni con altri geni protein-coding.

Serve per assegnare in modo univoco le regioni ChIP ai geni: un gene viene escluso se il suo intervallo
(gene intero, righe 'gene' del GTF) si sovrappone anche di 1 bp a quello di un altro gene protein-coding,
su qualsiasi filamento. Esclusi chrY e chrM.
L'intersezione per nome con il set consensus del passo 8 (geni con picchi) si fa nel passo 12b.

Uso: python 12a_nonoverlapping_genes.py gencode.v48.basic.annotation.gtf.gz genes_pc_nonoverlap.bed
Output: BED4 (cromosoma, inizio 0-based, fine, nome del gene), ordinato per coordinate
"""
import gzip, re, sys
import pandas as pd

gtf, out = sys.argv[1], sys.argv[2]
op = gzip.open if gtf.endswith(".gz") else open
get = lambda s, k: (m.group(1) if (m := re.search(rf'{k} "([^"]+)"', s)) else None)

rows = []
with op(gtf, "rt") as fh:
    for l in fh:
        if l.startswith("#"):
            continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene" or get(f[8], "gene_type") != "protein_coding":
            continue
        if f[0] in ("chrY", "chrM", "Y", "M", "MT"):
            continue
        rows.append((f[0], int(f[3]) - 1, int(f[4]), get(f[8], "gene_name") or get(f[8], "gene_id")))

G = pd.DataFrame(rows, columns=["chrom", "start", "end", "name"]).sort_values(["chrom", "start", "end"]).reset_index(drop=True)

# un gene si sovrappone a un altro se inizia prima della fine massima dei geni precedenti (stesso cromosoma)
# oppure se il gene successivo inizia prima della sua fine
keep = pd.Series(True, index=G.index)
for _, g in G.groupby("chrom", sort=False):
    idx = g.index.to_numpy(); s = g.start.to_numpy(); e = g.end.to_numpy()
    prev_max_end = pd.Series(e).cummax().shift(1, fill_value=-1).to_numpy()
    next_start = pd.Series(s).shift(-1, fill_value=10**12).to_numpy()
    overlap = (s < prev_max_end) | (next_start < e)
    keep[idx[overlap]] = False

# nomi duplicati (rari): si tiene la prima copia
out_df = G[keep].drop_duplicates("name")
out_df.to_csv(out, sep="\t", header=False, index=False)
print(f"geni protein-coding: {len(G)} | senza sovrapposizioni: {len(out_df)} -> {out}")
