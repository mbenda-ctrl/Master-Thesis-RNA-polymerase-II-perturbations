#!/usr/bin/env python3
"""
PASSO 2 - Filtro dei siti modkit per campione: copertura N_valid_cov >= 98 e fraction_modified >= 20 %.

Legge i file per campione <...>HEK_<CONDIZIONE>_<r|R><n>_all_mods_filtCov_annotated.tsv (output modkit
annotato), traduce i codici di modifica (a = m6A, m = m5C, 17596 = inosina, 17802 = pseudouridina),
aggiunge regione genica e trascritto dalla tabella di lookup e scrive un file per campione.

Uso: python 02_filtra_siti.py --input-dir data/modkit --lookup lookup_regione_siti.csv.gz \
                              --outdir results/modkit/siti_filtrati_N98_pct20 [--min-cov 98 --min-pct 20]
"""
import argparse, pandas as pd, glob, re, os
ap = argparse.ArgumentParser()
ap.add_argument("--input-dir", default="data/modkit")
ap.add_argument("--lookup", required=True, help="tabella chrom/start/end/gene_id -> region, transcript_used, gene_symbol")
ap.add_argument("--outdir", default="results/modkit/siti_filtrati_N98_pct20")
ap.add_argument("--min-cov", type=int, default=98)
ap.add_argument("--min-pct", type=float, default=20)
a = ap.parse_args()
out = a.outdir; os.makedirs(out, exist_ok=True)
lk=pd.read_csv(a.lookup)
lk=lk.drop_duplicates(["chrom","start","end","gene_id"])
names={'a':'m6A','m':'m5C','17596':'Inosina','17802':'Pseudouridina'}
rows=[]
for f in sorted(glob.glob(os.path.join(a.input_dir, "*_all_mods_filtCov_annotated.tsv"))):
    m=re.search(r"HEK_(\w+?)_[rR](\d)_all_mods",f); cond,rep=m.group(1).upper() if m.group(1).upper() in("STORM","WT","NACL") else m.group(1),"R"+m.group(2)
    d=pd.read_csv(f,sep="\t",dtype={'mod_code':str})
    k=d[(d.N_valid_cov>=a.min_cov)&(d.fraction_modified_pct>=a.min_pct)].copy()
    k["mod_name"]=k.mod_code.map(names)
    k=k.merge(lk[["chrom","start","end","gene_id","region","transcript_used","gene_symbol"]],on=["chrom","start","end","gene_id"],how="left")
    k.insert(0,"replicate",rep); k.insert(0,"condition",cond)
    k.to_csv(f"{out}/siti_{cond}_{rep}_N98_pct20.tsv.gz",sep="\t",index=False)
    c=k.mod_name.value_counts()
    rows.append(dict(condition=cond,replicate=rep,siti_grezzi=len(d),siti_filtrati=len(k),**{n:int(c.get(n,0)) for n in names.values()}))
r=pd.DataFrame(rows); r.to_csv(f"{out}/riepilogo_siti_filtrati.csv",index=False); print(r.to_string())
print("region NaN %:",round(k.region.isna().mean()*100,1))
