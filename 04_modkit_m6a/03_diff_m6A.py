#!/usr/bin/env python3
"""
PASSO 3 - Analisi differenziale m6A per gene (geni protein-coding), trattato vs controllo.

Per ogni gene presente in entrambe le condizioni (siti filtrati del passo 2):
  - stechiometria: reads modificate / copertura sommate sui siti del gene, test quasi-binomiale
    (GLM binomiale con dispersione stimata, minimo 3 osservazioni) oppure Fisher sui totali;
    FDR di Benjamini-Hochberg. Gene differenzialmente metilato (colonna diff, file
    diff_metilati_log2FC058_*): |log2FC stechiometria| > 0.58 e FDR < 0.05
    (stech_diff, criterio precedente: FDR < 0.05 e |delta %| >= 10)
  - numero di siti: siti_diff se la differenza e' di almeno 2 siti e i replicati sono coerenti
  - carico (somma delle frazioni modificate): carico_diff se |delta| >= 1 con replicati coerenti
Confronti: STORM vs WT, NACL vs WT, E11 (Fast) vs Amr, C412 (Slow) vs Amr.

Uso: python 03_diff_m6A.py --sites-dir results/modkit/siti_filtrati_N98_pct20 --outdir results/modkit/diff_m6A
"""
import argparse, os
import pandas as pd, numpy as np, warnings; warnings.filterwarnings("ignore")
import statsmodels.api as sm
from scipy.stats import fisher_exact, norm
from statsmodels.stats.multitest import multipletests
ap = argparse.ArgumentParser()
ap.add_argument("--sites-dir", default="results/modkit/siti_filtrati_N98_pct20")
ap.add_argument("--outdir", default="results/modkit/diff_m6A")
a = ap.parse_args()
D = a.sites_dir.rstrip("/") + "/"; O = a.outdir.rstrip("/") + "/"; os.makedirs(O, exist_ok=True)
agg={}
for c in ["WT","STORM","NACL","Amr","E11","C412"]:
    for r in ("R1","R2"):
        d=pd.read_csv(f"{D}siti_{c}_{r}_N98_pct20.tsv.gz",sep="\t",usecols=["mod_name","gene_id","gene_name","gene_type","N_mod","N_valid_cov","fraction_modified_pct"])
        d=d[(d.mod_name=="m6A")&(d.gene_type=="protein_coding")]
        agg[(c,r)]=d.groupby("gene_id").agg(n=("N_mod","size"),mod=("N_mod","sum"),cov=("N_valid_cov","sum"),carico=("fraction_modified_pct",lambda s:s.sum()/100),nome=("gene_name","first"))
summ=[]
for base,test in [("WT","STORM"),("WT","NACL"),("Amr","E11"),("Amr","C412")]:
    A=lambda c:set(agg[(c,"R1")].index)|set(agg[(c,"R2")].index)
    comuni=sorted(A(base)&A(test)); rows=[]
    for gid in comuni:
        obs=[]
        for c,x in ((base,0),(test,1)):
            for r in ("R1","R2"):
                a=agg[(c,r)]
                if gid in a.index: obs.append((x,a.at[gid,"mod"],a.at[gid,"cov"],a.at[gid,"n"],a.at[gid,"carico"],a.at[gid,"nome"]))
        x=np.array([o[0] for o in obs]); mod=np.array([o[1] for o in obs],float); cov=np.array([o[2] for o in obs],float)
        n=np.array([o[3] for o in obs],float); car=np.array([o[4] for o in obs],float)
        pb=mod[x==0].sum()/cov[x==0].sum()*100; pt=mod[x==1].sum()/cov[x==1].sum()*100
        p=np.nan; metodo="pooled_Fisher"
        if len(obs)>=3:
            try:
                m=sm.GLM(np.column_stack([mod,cov-mod]),sm.add_constant(x.astype(float)),family=sm.families.Binomial()).fit()
                disp=max(1.0,m.pearson_chi2/m.df_resid); p=2*norm.sf(abs(m.params[1]/(m.bse[1]*np.sqrt(disp)))); metodo="quasi_binomiale"
            except Exception: pass
        if np.isnan(p):
            tab=[[mod[x==0].sum(),cov[x==0].sum()-mod[x==0].sum()],[mod[x==1].sum(),cov[x==1].sum()-mod[x==1].sum()]]; p=fisher_exact(tab)[1]
        def coerente(v):
            b,t=v[x==0],v[x==1]
            return 1 if t.min()>b.max() else (-1 if t.max()<b.min() else 0)
        rows.append(dict(gene_id=gid,gene_name=obs[0][5],n_replicati=len(obs),
            siti_base=n[x==0].mean(),siti_test=n[x==1].mean(),delta_siti=n[x==1].mean()-n[x==0].mean(),siti_coerenti=coerente(n),
            pct_base=pb,pct_test=pt,delta_pct=pt-pb,metodo=metodo,p_stech=p,
            carico_base=car[x==0].mean(),carico_test=car[x==1].mean(),delta_carico=car[x==1].mean()-car[x==0].mean(),carico_coerente=coerente(car)))
    r=pd.DataFrame(rows); r["FDR_stech"]=multipletests(r.p_stech,method="fdr_bh")[1]
    r["stech_diff"]=np.where((r.FDR_stech<0.05)&(r.delta_pct.abs()>=10),np.sign(r.delta_pct),0).astype(int)
    r["siti_diff"]=np.where((r.delta_siti.abs()>=2)&(r.siti_coerenti!=0)&(np.sign(r.delta_siti)==r.siti_coerenti),np.sign(r.delta_siti),0).astype(int)
    r["carico_diff"]=np.where((r.delta_carico.abs()>=1)&(r.carico_coerente!=0)&(np.sign(r.delta_carico)==r.carico_coerente),np.sign(r.delta_carico),0).astype(int)
    r["categoria"]=np.select([(r.stech_diff!=0)&(r.siti_diff!=0)&(r.stech_diff==r.siti_diff),(r.stech_diff!=0)&(r.siti_diff!=0),(r.stech_diff!=0),(r.siti_diff!=0)],["stech+siti concordi","stech+siti discordi","solo stechiometria","solo n_siti"],"nessuna")
    r.sort_values("p_stech").to_csv(O+f"diff_completo_{test}_vs_{base}.csv",index=False)
    # classificazione finale dei geni differenzialmente metilati: |log2FC stechiometria| > 0.58 e FDR < 0.05
    r["log2FC_stech"]=np.log2(r.pct_test/r.pct_base)
    r["diff"]=np.where((r.log2FC_stech.abs()>0.58)&(r.FDR_stech<0.05),np.sign(r.log2FC_stech),0).astype(int)
    r.sort_values("p_stech").to_csv(O+f"diff_metilati_log2FC058_{test}_vs_{base}.csv",index=False)
    # persi/guadagnati
    rb=set(agg[(base,"R1")].index)&set(agg[(base,"R2")].index); rt=set(agg[(test,"R1")].index)&set(agg[(test,"R2")].index)
    s=dict(confronto=f"{test}_vs_{base}",geni=len(r),FDR05=int((r.FDR_stech<0.05).sum()),stech_su=int((r.stech_diff>0).sum()),stech_giu=int((r.stech_diff<0).sum()),siti_su=int((r.siti_diff>0).sum()),siti_giu=int((r.siti_diff<0).sum()),carico_su=int((r.carico_diff>0).sum()),carico_giu=int((r.carico_diff<0).sum()),persi=len(rb-A(test)),guadagnati=len(rt-A(base)))
    s.update(r.categoria.value_counts().to_dict()); summ.append(s)
S=pd.DataFrame(summ).fillna(0); S.to_csv(O+"riepilogo_diff_completo.csv",index=False); print(S.T.to_string())
