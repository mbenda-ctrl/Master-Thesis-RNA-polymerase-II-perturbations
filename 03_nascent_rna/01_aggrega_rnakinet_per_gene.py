#!/usr/bin/env python3
"""
Aggrega il file RNAkinet per-read (read_id, 5eu_mod_score, 5eu_modified_prediction,
gene_id, gene_name, sample) a livello di GENE e CAMPIONE.

Lavora a blocchi (chunk) per non caricare tutto il file in RAM in una volta —
necessario per file da ~1GB con decine di milioni di righe.

USO:
    python3 aggrega_rnakinet_per_gene.py /percorso/al/file_rnakinet.tsv /percorso/output_aggregato.csv

Se il separatore non e' il tab, cambia SEP sotto (es. ',' per CSV).
"""
import sys
import pandas as pd

SEP = '\t'          # cambia in ',' se il file e' separato da virgole
CHUNKSIZE = 2_000_000  # righe per blocco; riduci a 500_000 se il laptop fatica

def main(path_input, path_output):
    accumulatore = {}  # (sample, gene_id, gene_name) -> [n_reads, n_modificate, somma_score]

    n_blocchi = 0
    for chunk in pd.read_csv(path_input, sep=SEP, chunksize=CHUNKSIZE,
                              usecols=['sample', 'gene_id', 'gene_name',
                                       '5eu_mod_score', '5eu_modified_prediction']):
        n_blocchi += 1
        print(f'Blocco {n_blocchi} letto ({len(chunk):,} righe)...', flush=True)

        # assicura che la colonna booleana sia trattata correttamente anche se letta come stringa
        if chunk['5eu_modified_prediction'].dtype == object:
            chunk['5eu_modified_prediction'] = chunk['5eu_modified_prediction'].astype(str).str.strip().str.lower() == 'true'

        gruppo = chunk.groupby(['sample', 'gene_id', 'gene_name']).agg(
            n_reads=('5eu_mod_score', 'size'),
            n_modificate=('5eu_modified_prediction', 'sum'),
            somma_score=('5eu_mod_score', 'sum'),
        )

        for idx, riga in gruppo.iterrows():
            if idx not in accumulatore:
                accumulatore[idx] = [0, 0, 0.0]
            accumulatore[idx][0] += riga['n_reads']
            accumulatore[idx][1] += riga['n_modificate']
            accumulatore[idx][2] += riga['somma_score']

    print(f'Totale blocchi processati: {n_blocchi}. Costruisco la tabella finale...')

    righe_finali = []
    for (sample, gene_id, gene_name), (n_reads, n_mod, somma_score) in accumulatore.items():
        righe_finali.append({
            'sample': sample,
            'gene_id': gene_id,
            'gene_name': gene_name,
            'n_reads': n_reads,
            'n_reads_modificate_5eu': n_mod,
            'pct_modificate_5eu': round(100 * n_mod / n_reads, 3) if n_reads > 0 else None,
            'score_medio_5eu': round(somma_score / n_reads, 4) if n_reads > 0 else None,
        })

    df_finale = pd.DataFrame(righe_finali).sort_values(['sample', 'n_reads'], ascending=[True, False])
    df_finale.to_csv(path_output, index=False)
    print(f'Salvato: {path_output} ({len(df_finale):,} righe gene x campione)')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Uso: python3 aggrega_rnakinet_per_gene.py <file_input> <file_output.csv>')
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
