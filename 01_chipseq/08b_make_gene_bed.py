#!/usr/bin/env python3
"""
PASSO 8b - Costruisce un BED6 (e una tabella con TSS/TES) di UN trascritto per gene, a partire da un GTF
e da una tabella di geni (gene_set_*.tsv del passo 08a: colonna gene_id o ENSEMBL).

Trascritto per gene: MANE Select (tag "MANE_Select" del GTF GENCODE basic). I geni senza MANE Select
sono esclusi; con --allow-fallback si usa invece Ensembl_canonical e poi il trascritto piu lungo.
Funziona con GTF GENCODE, Ensembl e MANE (Ensembl). Legge anche .gz.

Uso:
  python 08b_make_gene_bed.py gencode.v48.basic.annotation.gtf.gz gene_set_union_n2.tsv \
         consensus/genes_pc_union_MANE_n2_min1kb --min-len 1000 --biotype protein_coding
                          [--biotype protein_coding] [--add-chr]
"""
import argparse, gzip, re, sys
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("gtf"); ap.add_argument("genes"); ap.add_argument("prefix")
ap.add_argument("--min-len", type=int, default=0, help="scarta trascritti piu corti (bp)")
ap.add_argument("--biotype", default=None, help="es. protein_coding (usa gene_type/gene_biotype se presente)")
ap.add_argument("--allow-fallback", action="store_true", help="tiene anche i geni senza MANE Select (Ensembl_canonical, poi il piu lungo)")
ap.add_argument("--add-chr", action="store_true", help="aggiunge 'chr' ai cromosomi (GTF Ensembl: 1,2,...,MT)")
a = ap.parse_args()

opener = gzip.open if a.gtf.endswith(".gz") else open
attr = lambda s, k: (m.group(1) if (m := re.search(rf'{k} "([^"]+)"', s)) else None)

tx = {}   # gene_id (senza versione) -> lista trascritti
with opener(a.gtf, "rt") as fh:
    for l in fh:
        if l[0] == "#": continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "transcript": continue
        gid = attr(f[8], "gene_id")
        if gid is None: continue
        gid = gid.split(".")[0]
        tags = set(re.findall(r'tag "([^"]+)"', f[8]))
        chrom = f[0]
        if a.add_chr:
            chrom = "chrM" if chrom == "MT" else ("chr" + chrom if not chrom.startswith("chr") else chrom)
        tx.setdefault(gid, []).append(dict(
            chrom=chrom, start=int(f[3]), end=int(f[4]), strand=f[6],
            tx=attr(f[8], "transcript_id"), name=attr(f[8], "gene_name"),
            biotype=attr(f[8], "gene_type") or attr(f[8], "gene_biotype"), tags=tags))

def pick(lst):
    if not a.allow_fallback:
        c = [t for t in lst if "MANE_Select" in t["tags"]]
        return (max(c, key=lambda t: t["end"] - t["start"]), "MANE_Select") if c else (None, "no_MANE")
    for tag in ("MANE_Select", "Ensembl_canonical"):
        c = [t for t in lst if tag in t["tags"]]
        if c: return max(c, key=lambda t: t["end"] - t["start"]), tag
    return max(lst, key=lambda t: t["end"] - t["start"]), "longest"

G = pd.read_csv(a.genes, sep="\t")
out, n = [], dict(no_ensembl=0, not_in_gtf=0, no_MANE=0, biotype=0, too_short=0, ok=0)
for _, r in G.iterrows():
    e = r.get("ENSEMBL", r.get("gene_id"))                 # gene_set_*.tsv del passo 08a: colonna gene_id
    if pd.isna(e): n["no_ensembl"] += 1; continue
    lst = tx.get(str(e).split(".")[0])
    if not lst: n["not_in_gtf"] += 1; continue
    t, why = pick(lst)
    if t is None: n["no_MANE"] += 1; continue
    if a.biotype and t["biotype"] and t["biotype"] != a.biotype: n["biotype"] += 1; continue
    L = t["end"] - t["start"] + 1
    if L < a.min_len: n["too_short"] += 1; continue
    tss, tes = (t["start"], t["end"]) if t["strand"] == "+" else (t["end"], t["start"])
    out.append(dict(chrom=t["chrom"], start0=t["start"] - 1, end=t["end"],
                    name=t["name"] or r.get("symbol", r.get("gene_name")), score=0, strand=t["strand"],
                    EntrezID=r.get("EntrezID"), ENSEMBL=e, transcript=t["tx"], selected_by=why,
                    TSS=tss, TES=tes, length=L, n_samples=r.get("n_samples")))
    n["ok"] += 1

D = pd.DataFrame(out).sort_values(["chrom", "start0"])
D[["chrom", "start0", "end", "name", "score", "strand"]].to_csv(a.prefix + ".bed", sep="\t", header=False, index=False)
D.to_csv(a.prefix + ".tsv", sep="\t", index=False)
print(f"geni in input: {len(G)} | nel BED: {n['ok']}")
print(f"scartati -> senza ID Ensembl: {n['no_ensembl']}, non nel GTF: {n['not_in_gtf']}, senza MANE Select: {n['no_MANE']}, "
      f"biotype: {n['biotype']}, troppo corti: {n['too_short']}")
print("trascritto scelto da:", D.selected_by.value_counts().to_dict())
