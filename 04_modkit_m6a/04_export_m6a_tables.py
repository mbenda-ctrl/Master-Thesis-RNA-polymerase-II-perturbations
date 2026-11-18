#!/usr/bin/env python3
"""
PASSO 4 - Tabelle m6A per l'integrazione (Excel letti da 05_integration).

Dai risultati del passo 3 (diff_metilati_log2FC058_<test>_vs_<base>.csv) e dai siti filtrati del passo 2:
  1. categoria di ogni gene: stech+siti concordi / discordi, solo stechiometria, solo n. siti, nessuna
  2. geni persi / guadagnati: m6A in entrambi i replicati di una condizione e in nessun replicato dell'altra
  3. geni_diff_metilati_e_perdita_siti.xlsx  -> fogli <Confronto>_diff_metilati, _meno_siti, _persi
     (letto da 05_integration/03_m6a_polII_venn.R)
  4. geni_variazione_siti_e_stechiometria_log2FC.xlsx -> un foglio per confronto con log2FC del numero
     di siti e della stechiometria (letto dalle heatmap con --m6a)

Uso: python 04_export_m6a_tables.py --sites-dir results/modkit/siti_filtrati_N98_pct20 \
                                    --diff-dir results/modkit/diff_m6A
"""
import argparse, os
import numpy as np, pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--sites-dir", default="results/modkit/siti_filtrati_N98_pct20")
ap.add_argument("--diff-dir", default="results/modkit/diff_m6A", help="output del passo 3 (qui si scrivono anche gli Excel)")
a = ap.parse_args()
O = a.diff_dir.rstrip("/") + "/"

# file del passo 3 -> nome del confronto usato a valle
COMP = [("STORM", "WT", "Storm_vs_WT"), ("NACL", "WT", "NaCl_vs_WT"), ("E11", "Amr", "Fast_vs_Amr"), ("C412", "Amr", "Slow_vs_Amr")]

def genes_m6a(c, r):
    f = os.path.join(a.sites_dir, f"siti_{c}_{r}_N98_pct20.tsv.gz")
    d = pd.read_csv(f, sep="\t", usecols=["mod_name", "gene_id", "gene_name", "gene_type"])
    d = d[(d.mod_name == "m6A") & (d.gene_type == "protein_coding")]
    return d.drop_duplicates("gene_id").set_index("gene_id").gene_name

cols = ["gene_id", "gene_name", "n_replicati", "siti_base", "siti_test", "delta_siti", "pct_base", "pct_test", "delta_pct",
        "log2FC_stech", "p_stech", "FDR_stech", "carico_base", "carico_test", "delta_carico", "categoria"]
cnt = []
w1 = pd.ExcelWriter(O + "geni_diff_metilati_e_perdita_siti.xlsx")
w2 = pd.ExcelWriter(O + "geni_variazione_siti_e_stechiometria_log2FC.xlsx")
for test, base, n in COMP:
    r = pd.read_csv(O + f"diff_metilati_log2FC058_{test}_vs_{base}.csv")
    # 1. categoria con la classificazione finale della stechiometria (colonna diff)
    s_ok, d = r.siti_diff != 0, r["diff"] != 0
    r["categoria"] = np.select([d & s_ok & (np.sign(r["diff"]) == r.siti_diff), d & s_ok, d, s_ok],
                               ["stech+siti concordi", "stech+siti discordi", "solo stechiometria", "solo n. siti"], "nessuna")
    r.to_csv(O + f"diff_metilati_log2FC058_{test}_vs_{base}.csv", index=False)

    # 2. geni persi / guadagnati
    gb = {k: genes_m6a(base, k) for k in ("R1", "R2")}; gt = {k: genes_m6a(test, k) for k in ("R1", "R2")}
    names = pd.concat(list(gb.values()) + list(gt.values())); names = names[~names.index.duplicated()]
    persi = sorted((set(gb["R1"].index) & set(gb["R2"].index)) - (set(gt["R1"].index) | set(gt["R2"].index)))
    guad = sorted((set(gt["R1"].index) & set(gt["R2"].index)) - (set(gb["R1"].index) | set(gb["R2"].index)))
    pg = pd.DataFrame({"gene_id": persi + guad, "gene_name": [names.get(g) for g in persi + guad],
                       "stato": [f"persi_in_{test}"] * len(persi) + [f"guadagnati_in_{test}"] * len(guad)})
    pg.to_csv(O + f"geni_persi_guadagnati_{test}_vs_{base}.csv", index=False)

    # 3. Excel per il Venn m6A x Pol II
    dm = r[r["diff"] != 0].copy()
    dm.insert(2, "direzione", dm["diff"].map({1: "stechiometria_su", -1: "stechiometria_giu"}))
    ms = r[r.siti_diff < 0].copy()
    dm[["gene_id", "gene_name", "direzione"] + cols[2:]].sort_values("FDR_stech").to_excel(w1, sheet_name=f"{n}_diff_metilati"[:31], index=False)
    ms[cols].sort_values("delta_siti").to_excel(w1, sheet_name=f"{n}_meno_siti"[:31], index=False)
    pg[pg.stato.str.startswith("persi")].to_excel(w1, sheet_name=f"{n}_persi"[:31], index=False)

    # 4. Excel per le heatmap (variazione del numero di siti e della stechiometria)
    r["log2FC_n_siti"] = np.log2(r.siti_test / r.siti_base)
    r["variazione_n_siti"] = r.siti_diff.map({1: "aumenta", -1: "diminuisce", 0: "invariato"})
    r["variazione_stechiometria"] = r["diff"].map({1: "aumenta", -1: "diminuisce", 0: "invariata"})
    s = r[(r.siti_diff != 0) | (r["diff"] != 0)]
    s[["gene_id", "gene_name", "variazione_n_siti", "siti_base", "siti_test", "delta_siti", "log2FC_n_siti",
       "variazione_stechiometria", "pct_base", "pct_test", "delta_pct", "log2FC_stech", "p_stech", "FDR_stech", "categoria"]
      ].sort_values(["variazione_n_siti", "variazione_stechiometria", "log2FC_stech"]).to_excel(w2, sheet_name=n, index=False)
    cnt.append(dict(confronto=n, diff_metilati=len(dm), meno_siti=len(ms), persi=len(persi), guadagnati=len(guad),
                    variazione_siti_o_stech=len(s)))
pd.DataFrame(cnt).to_excel(w1, sheet_name="LEGGIMI_conteggi", index=False)
w1.close(); w2.close()
print(pd.DataFrame(cnt).to_string(index=False))
