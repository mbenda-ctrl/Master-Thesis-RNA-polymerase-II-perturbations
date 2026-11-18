#!/usr/bin/env python3
"""
modkit_analysis.py
===================

Pipeline per l'analisi comparativa di siti di modifica dell'RNA
(m6A, m5C, Pseudouridina, Inosina) chiamati da MODKIT (formato bedMethyl
esteso, annotato con gene_id/gene_name/gene_type).

COSA FA:
1. Carica tutti i file *_annotated.tsv di una cartella, ne deduce
   condizione/replicato dal nome file e li concatena in un unico dataframe
   "long" con una colonna `mod_name` leggibile (m6A, m5C, Pseudouridina,
   Inosina) al posto del `mod_code` grezzo.
2. Calcola i totali di siti modificati per: modifica, campione, condizione
   (media tra replicati), gene_type.
3. Classifica ogni sito in base alla % di modifica (fraction_modified_pct)
   in bassa / media / alta (soglie configurabili).
4. (OPZIONALE, richiede un file di annotazione GTF/GFF3) Classifica ogni
   sito nella regione genica di appartenenza: 5'UTR, CDS/esone codificante,
   esone non codificante, introne, 3'UTR - usando il gene_id già presente
   nel file MODKIT.
5. Produce una tabella riassuntiva per gene ("stato di modifica del gene":
   quante e quali modifiche, in quale regione, con quale livello).
6. Salva tabelle di riepilogo (CSV) e alcuni grafici diagnostici (PNG).

USO:
    python modkit_analysis.py \
        --input-dir /path/to/tsv_files \
        --outdir /path/to/results \
        --gtf /path/to/annotation.gtf.gz   # opzionale ma consigliato

Se --gtf non viene fornito, tutti gli step che non richiedono l'annotazione
di regione vengono comunque eseguiti; la colonna `region` risultera' NA.
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# 1. COSTANTI
# --------------------------------------------------------------------------

# Mappa dei codici MODKIT/ChEBI verso un nome leggibile.
# a = m6A, m = m5C, 17596 = Inosina, 17802 = Pseudouridina
# (rif. documentazione modkit / SAM MM/ML tag spec)
MOD_MAP = {
    "a": "m6A",
    "m": "m5C",
    "17596": "Inosina",
    "17802": "Pseudouridina",
}

BEDMETHYL_COLS = [
    "chrom", "start", "end", "mod_code", "score", "strand",
    "thick_start", "thick_end", "color", "N_valid_cov",
    "fraction_modified_pct", "N_mod", "N_canonical", "N_other_mod",
    "N_delete", "N_fail", "N_diff", "N_nocall",
    "gene_id", "gene_name", "gene_type",
]

# Soglie di default per la classificazione del livello di modifica.
# (vedi MOD_LEVEL_BINS/MOD_LEVEL_LABELS piu' sotto per i quartili usati)


# --------------------------------------------------------------------------
# 2. PARSING DEL NOME FILE -> campione / condizione / replicato
# --------------------------------------------------------------------------

def parse_sample_info(filename: str):
    """
    Deduce condizione e replicato dal nome file, es:
    'HEK_Amr_r1_all_mods_filtCov_annotated.tsv' -> ('Amr', 'r1')
    'HEK_STORM_R1_all_mods_filtCov_annotated.tsv' -> ('STORM', 'R1')

    Adatta questa funzione se la nomenclatura dei tuoi file cambia.
    """
    stem = Path(filename).stem
    m = re.search(r"HEK_(?P<condition>[A-Za-z0-9]+)_(?P<rep>[Rr]\d+)", stem)
    if not m:
        # fallback: usa l'intero nome come condizione, replicato ignoto
        return stem, "NA"
    condition = m.group("condition")
    rep = m.group("rep").upper()
    return condition, rep


# --------------------------------------------------------------------------
# 3. CARICAMENTO DEI DATI
# --------------------------------------------------------------------------

def load_one_file(filepath: Path) -> pd.DataFrame:
    df = pd.read_csv(filepath, sep="\t")
    # alcuni file potrebbero non avere header identico -> normalizza
    df.columns = [c.strip() for c in df.columns]

    if "mod_code" not in df.columns:
        sys.exit(
            f"ERRORE: '{filepath.name}' non ha la colonna 'mod_code' -> non sembra un file "
            f"di output modkit valido (colonne trovate: {list(df.columns)[:8]}...). "
            f"Controlla che nella cartella di input non ci siano altri file .tsv "
            f"che non siano output modkit (es. tabelle di annotazione con un nome diverso "
            f"da quelli in ANNOTATION_TABLE_NAMES)."
        )

    condition, rep = parse_sample_info(filepath.name)
    df["sample_id"] = filepath.stem
    df["condition"] = condition
    df["replicate"] = rep
    df["mod_name"] = df["mod_code"].astype(str).map(MOD_MAP).fillna(df["mod_code"].astype(str))

    # riduzione memoria: le colonne testuali si ripetono su milioni di righe
    # ma hanno poche migliaia di valori distinti -> category. NOTA: gene_id e
    # gene_name sono ESCLUSI di proposito: hanno decine di migliaia di valori
    # distinti, e alcune operazioni pandas (es. .pivot(), che a differenza di
    # .groupby() NON ha il parametro observed=) possono "esplodere" in memoria
    # enumerando combinazioni di categorie mai osservate nei dati. Il
    # risparmio di RAM su queste due colonne non vale il rischio: restano
    # semplici stringhe.
    for col in ["chrom", "mod_code", "mod_name", "strand",
                "gene_type", "sample_id", "condition", "replicate"]:
        if col in df.columns:
            df[col] = df[col].astype("category")
    if "fraction_modified_pct" in df.columns:
        df["fraction_modified_pct"] = df["fraction_modified_pct"].astype("float32")
    for col in ["start", "end", "N_valid_cov", "N_mod", "N_canonical", "N_other_mod",
                "N_delete", "N_fail", "N_diff", "N_nocall", "score"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("float32").astype("Int32")

    return df


# Nomi delle tabelle di annotazione leggere: se si trovano nella STESSA
# cartella dei file modkit (uso comune/comodo), vanno escluse dal caricamento
# -- altrimenti load_all() le scambierebbe per output modkit e fallirebbe
# (mancano le colonne come mod_code).
ANNOTATION_TABLE_NAMES = {"mane_transcript_spans.tsv", "mane_region_features.tsv", "all_gene_spans.tsv"}


def load_all(input_dir: Path) -> pd.DataFrame:
    files = sorted(
        f for f in input_dir.glob("*.tsv")
        if f.name not in ANNOTATION_TABLE_NAMES
    )
    if not files:
        sys.exit(f"Nessun file .tsv trovato in {input_dir}")
    frames = []
    for f in files:
        print(f"  caricamento {f.name} ...")
        frames.append(load_one_file(f))
    df = pd.concat(frames, ignore_index=True)
    del frames
    # ricompatta le category dopo la concat (ogni file aveva le proprie);
    # gene_id/gene_name NON sono tra queste, vedi nota in load_one_file
    for col in ["chrom", "mod_code", "mod_name", "strand",
                "gene_type", "sample_id", "condition", "replicate"]:
        if col in df.columns:
            df[col] = df[col].astype("category")
    print(f"Totale righe caricate: {len(df):,} da {len(files)} file "
          f"(~{df.memory_usage(deep=True).sum() / 1e9:.2f} GB in RAM)")
    return df


# --------------------------------------------------------------------------
# 3bis. FILTRO DI QUALITA' DEI SITI
# --------------------------------------------------------------------------

def filter_modkit_sites(df: pd.DataFrame, min_coverage: int = 10,
                         m6a_fail_fraction: float = 0.25,
                         adaptive_percentile: float = 90.0) -> pd.DataFrame:
    """
    Filtra i siti MODKIT in base alla qualita' della chiamata.

    - N_valid_cov >= min_coverage (default 10): copertura minima, uguale
      per tutte le modifiche (in pratica non e' mai il fattore limitante:
      la coverage mediana e' 52-65x per tutte e 4 le modifiche).

    - fail_fraction (= N_fail / (N_valid_cov + N_fail)): qui la soglia NON
      e' fissa uguale per tutte le modifiche. Una soglia fissa unica
      (es. 0.25) penalizza in modo fortemente asimmetrico le modifiche
      diverse da m6A: verificato sui dati, la fail_fraction MEDIANA e'
      0.08 per m6A ma 0.29-0.65 per Inosina/m5C/Pseudouridina -- una
      soglia assoluta tarata (anche implicitamente) su m6A scarterebbe
      oltre il 90% dei siti m5C. Questo e' coerente con la letteratura
      (vedi ModkitOpt, bioRxiv 2025): i modelli di chiamata per le
      modifiche piu' rare producono intrinsecamente piu' chiamate a
      bassa confidenza per singola lettura.

      Percio':
      * per m6A si usa una soglia FISSA (m6a_fail_fraction, default 0.25),
        che per questa modifica corrisponde comunque a un percentile alto
        (~95%) della sua distribuzione, quindi e' un filtro sensato in
        senso assoluto.
      * per le altre modifiche (m5C, Pseudouridina, Inosina) si usa una
        soglia ADATTIVA: il percentile `adaptive_percentile` (default 90)
        calcolato SEPARATAMENTE sulla distribuzione di fail_fraction di
        QUELLA modifica. Questo scarta il 10% peggiore di ciascuna
        modifica rispetto al proprio comportamento tipico, invece di
        applicare un righello assoluto tarato implicitamente su m6A.
    """
    total_cov = df["N_valid_cov"] + df["N_fail"]
    fail_fraction = df["N_fail"] / total_cov.replace(0, 1)
    before = len(df)

    cov_ok = df["N_valid_cov"] >= min_coverage

    is_m6a = df["mod_name"] == "m6A"
    m6a_ok = is_m6a & (fail_fraction <= m6a_fail_fraction)

    # soglia adattiva calcolata SOLO sulle righe delle altre modifiche
    other_fail = fail_fraction[~is_m6a]
    other_mod = df.loc[~is_m6a, "mod_name"]
    adaptive_thr_per_row = other_fail.groupby(other_mod, observed=True).transform(
        lambda s: s.quantile(adaptive_percentile / 100)
    )
    other_ok = pd.Series(False, index=df.index)
    other_ok.loc[~is_m6a] = other_fail <= adaptive_thr_per_row

    keep = cov_ok & (m6a_ok | other_ok)
    out = df[keep].copy()

    print(f"  Filtro qualita' (coverage>={min_coverage}; m6A: fail_fraction<={m6a_fail_fraction} fisso; "
          f"altre modifiche: fail_fraction<=percentile {adaptive_percentile:g} PROPRIO di ciascuna): "
          f"{before:,} -> {len(out):,} righe ({len(out)/before:.1%} tenute)")
    for mod in sorted(df["mod_name"].astype(str).unique()):
        n_before = (df["mod_name"] == mod).sum()
        n_after = (out["mod_name"] == mod).sum()
        if n_before > 0:
            print(f"    - {mod}: {n_before:,} -> {n_after:,} ({n_after/n_before:.1%} tenute)")

    # rete di sicurezza: rimuove le categorie non piu' osservate dopo il
    # filtro (basso rischio con le poche colonne rimaste categoriche, ma
    # costa pochissimo e previene ulteriori sorprese)
    for col in out.select_dtypes(include="category").columns:
        out[col] = out[col].cat.remove_unused_categories()

    return out


# --------------------------------------------------------------------------
# 4. CLASSIFICAZIONE PER % DI MODIFICA (quartili)
# --------------------------------------------------------------------------

MOD_LEVEL_BINS = [0, 25, 50, 75, 100.0001]  # 100.0001 per includere il 100% nell'ultimo bin
MOD_LEVEL_LABELS = ["0-25%", "25-50%", "50-75%", "75-100%"]


def classify_mod_level(df: pd.DataFrame) -> pd.DataFrame:
    df["mod_level"] = pd.cut(
        df["fraction_modified_pct"],
        bins=MOD_LEVEL_BINS,
        labels=MOD_LEVEL_LABELS,
        right=False,
        include_lowest=True,
    )
    return df


# --------------------------------------------------------------------------
# 5. ANNOTAZIONE DI REGIONE (esone/introne/5'UTR/3'UTR) - richiede GTF
# --------------------------------------------------------------------------

FEATURES_OF_INTEREST = {"gene", "transcript", "exon", "CDS", "five_prime_utr", "three_prime_utr"}


def _extract_attr(attr_field: str, key: str) -> str:
    # attr_field es: 'gene_id "ENSG..."; transcript_id "ENST..."; ...'
    k = f'{key} "'
    i = attr_field.find(k)
    if i == -1:
        return None
    i += len(k)
    j = attr_field.find('"', i)
    return attr_field[i:j]


class RegionAnnotator:
    """
    Costruisce, a partire da un GTF (Ensembl o GENCODE, anche .gz), un
    indice leggero per classificare ogni sito nella regione genica corretta.

    IMPORTANTE: la classificazione avviene su un SINGOLO trascritto
    rappresentativo per gene (non sull'unione di tutti i trascritti), per
    evitare che un esone alternativo di un'isoforma minoritaria faccia
    classificare come "esone" una posizione che nel trascritto di
    riferimento e' introne. Il trascritto scelto e', in ordine di
    preferenza:
      1) quello con tag "MANE_Select" (standard clinico, solo coding)
      2) quello con tag "Ensembl_canonical" (definito per ogni gene)
      3) in mancanza di entrambi, il trascritto piu' lungo del gene

    La classificazione si basa su gene_id + posizione, non sul nome del
    cromosoma, perche' modkit e il GTF possono usare notazioni diverse
    (es. "chr1" vs "1"). Il parsing e' fatto a mano (senza gffutils) per
    restare veloce anche su GTF dell'intero genoma con milioni di righe.
    """

    def __init__(self, gtf_path: Path):
        self.gene_to_tx, self.tx_index, self.gene_symbol, self.gene_coord_by_chrom = self._build_index(gtf_path)
        self._coord_trees = {}  # cache: cromosoma normalizzato -> IntervalTree (costruito on-demand)

    @staticmethod
    def _build_index(gtf_path: Path):
        opener = gzip.open if str(gtf_path).endswith(".gz") else open

        tx_span = {}          # transcript_id -> (start, end)
        tx_gene = {}           # transcript_id -> gene_id
        tx_strand = {}          # transcript_id -> '+'/'-'
        tx_feats = {}           # transcript_id -> [(start,end,type), ...]
        gene_best_tx = {}       # gene_id -> (priority, length, transcript_id)
        gene_symbol = {}        # gene_id -> gene_name (symbol), dalla riga "gene" del GTF
        gene_coord_by_chrom = {}  # cromosoma (nomenclatura Ensembl, es. "1","MT") -> [(start,end,gene_id), ...]

        def _priority(attr_field):
            if 'tag "MANE_Select"' in attr_field:
                return 2
            if 'tag "Ensembl_canonical"' in attr_field:
                return 1
            return 0

        n_lines = 0
        with opener(gtf_path, "rt") as f:
            for line in f:
                if line[0] == "#":
                    continue
                n_lines += 1
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                ftype = fields[2]
                if ftype not in FEATURES_OF_INTEREST:
                    continue
                attrs = fields[8]

                if ftype == "gene":
                    gene_id = _extract_attr(attrs, "gene_id")
                    gene_name = _extract_attr(attrs, "gene_name")
                    if gene_id is not None:
                        gene_symbol[gene_id] = gene_name if gene_name else gene_id
                        chrom = fields[0]
                        gstart, gend = int(fields[3]), int(fields[4])
                        gene_coord_by_chrom.setdefault(chrom, []).append((gstart, gend, gene_id))
                    continue

                tx_id = _extract_attr(attrs, "transcript_id")
                if tx_id is None:
                    continue
                start, end = int(fields[3]), int(fields[4])

                if ftype == "transcript":
                    gene_id = _extract_attr(attrs, "gene_id")
                    tx_span[tx_id] = (start, end)
                    tx_gene[tx_id] = gene_id
                    tx_strand[tx_id] = fields[6]
                    prio = _priority(attrs)
                    length = end - start
                    best = gene_best_tx.get(gene_id)
                    # preferisci priorita' piu' alta; a parita', trascritto piu' lungo
                    if best is None or (prio, length) > (best[0], best[1]):
                        gene_best_tx[gene_id] = (prio, length, tx_id)
                else:
                    tx_feats.setdefault(tx_id, []).append((start, end, ftype))

        n_tx_total = len(tx_span)
        gene_to_tx = {g: v[2] for g, v in gene_best_tx.items()}

        # IMPORTANTE per la memoria: teniamo in indice SOLO i trascritti
        # effettivamente scelti come rappresentativi di un gene (poche decine
        # di migliaia), non tutti i 500k+ trascritti del GTF. Altrimenti su
        # macchine con poca RAM il processo viene ucciso (OOM).
        used_tx_ids = set(gene_to_tx.values())
        tx_index = {
            tx_id: {"span": tx_span[tx_id], "strand": tx_strand[tx_id], "feats": tx_feats.get(tx_id, [])}
            for tx_id in used_tx_ids
        }
        del tx_span, tx_strand, tx_feats, tx_gene, gene_best_tx
        # ordiniamo per start: serve per la ricerca binaria nel lookup per coordinate
        for chrom in gene_coord_by_chrom:
            gene_coord_by_chrom[chrom].sort(key=lambda t: t[0])
        print(f"  GTF: {n_lines:,} righe rilevanti lette, {len(gene_to_tx):,} geni, "
              f"{n_tx_total:,} trascritti totali nel GTF ({len(tx_index):,} tenuti in indice), "
              f"{len(gene_symbol):,} gene symbol, {len(gene_coord_by_chrom)} cromosomi indicizzati per coordinate")
        return gene_to_tx, tx_index, gene_symbol, gene_coord_by_chrom

    @staticmethod
    def _normalize_chrom(chrom: str) -> str:
        """Converte la notazione del cromosoma usata da modkit (es. 'chr1',
        'chrM') in quella usata dal GTF Ensembl ('1', 'MT')."""
        c = str(chrom)
        if c.startswith("chr"):
            c = c[3:]
        if c in ("M", "MT"):
            c = "MT"
        return c

    def _get_tree(self, chrom_norm: str):
        if chrom_norm in self._coord_trees:
            return self._coord_trees[chrom_norm]
        from intervaltree import IntervalTree
        intervals = self.gene_coord_by_chrom.get(chrom_norm, [])
        tree = IntervalTree.from_tuples([(s, e + 1, gid) for s, e, gid in intervals])  # end esclusivo per IntervalTree
        self._coord_trees[chrom_norm] = tree
        return tree

    def lookup_by_coordinate(self, chrom: str, pos: int):
        """
        Cerca TUTTI i geni Ensembl (indipendentemente da cosa dica l'annotazione
        originale di modkit) che si sovrappongono a (chrom, pos). Utile per
        verificare/recuperare siti etichettati come "intergenica" o
        "gene_non_in_gtf": puo' darsi che un gene ci sia comunque, secondo
        Ensembl 115, anche se modkit (con una diversa versione/sorgente di
        annotazione) non lo aveva assegnato.
        Ritorna una lista di (gene_id, gene_symbol); vuota se nessun gene trovato.
        """
        chrom_norm = self._normalize_chrom(chrom)
        tree = self._get_tree(chrom_norm)
        hits = tree.at(pos)
        return [(iv.data, self.gene_symbol.get(iv.data, iv.data)) for iv in hits]

    def classify_site(self, pos: int, gene_id: str):
        """
        pos: coordinata 1-based del sito (colonna `end` del file MODKIT).
        Ritorna (region, transcript_id_usato). region e' una fra:
        5'UTR, 3'UTR, CDS, esone_ncRNA, introne (mappati con successo),
        intergenica (nessun gene secondo modkit stesso), gene_non_in_gtf
        (gene presente in modkit ma assente dal nostro GTF), o
        fuori_trascritto_canonico (gene trovato e presente nel GTF, ma la
        posizione cade FUORI dai confini genomici del trascritto canonico
        scelto come riferimento -- es. perche' modkit ha usato un'altra
        versione/isoforma con UTR piu' lunghe, o una diversa versione del
        gene model). Questi ultimi due casi sono siti che NON riusciamo a
        mappare in modo affidabile su un trascritto e vanno tenuti distinti
        dai veri siti intergenici.
        """
        if gene_id is None or gene_id == "." or pd.isna(gene_id):
            return "intergenica", None

        gene_id_base = str(gene_id).split(".")[0]  # rimuove versione ENSG...16 -> ENSG...
        tx_id = self.gene_to_tx.get(gene_id_base)
        if tx_id is None:
            return "gene_non_in_gtf", None

        tx = self.tx_index[tx_id]
        span = tx["span"]
        if not (span[0] <= pos <= span[1]):
            return "fuori_trascritto_canonico", tx_id

        hit_types = {t for (s, e, t) in tx["feats"] if s <= pos <= e}
        if "five_prime_utr" in hit_types:
            return "5'UTR", tx_id
        if "three_prime_utr" in hit_types:
            return "3'UTR", tx_id
        if "CDS" in hit_types:
            return "CDS", tx_id
        if "exon" in hit_types:
            return "esone_ncRNA", tx_id
        return "introne", tx_id


def annotate_regions(df: pd.DataFrame, gtf_path: Path):
    print("Indicizzazione del GTF (trascritto rappresentativo per gene)...")
    annot = RegionAnnotator(gtf_path)
    print("Classificazione della regione genica per ogni sito...")
    ends = df["end"].to_numpy()
    gene_ids = df["gene_id"].astype(str).to_numpy()  # da category a str solo per l'iterazione
    regions = [None] * len(df)
    tx_used = [None] * len(df)
    symbols = [None] * len(df)
    for i in range(len(df)):
        gid = gene_ids[i]
        r, t = annot.classify_site(int(ends[i]), gid)
        regions[i] = r
        tx_used[i] = t
        if gid not in (None, "."):
            symbols[i] = annot.gene_symbol.get(gid.split(".")[0])
    df["region"] = pd.Categorical(regions)
    # transcript_used e gene_symbol restano stringhe normali (NON category):
    # hanno decine di migliaia di valori distinti, stesso rischio di gene_id/
    # gene_name spiegato in load_one_file.
    df["transcript_used"] = tx_used
    df["gene_symbol"] = symbols
    return df, annot


# --------------------------------------------------------------------------
# 5bis-bis. ANNOTAZIONE MULTI-ISOFORMA (per coordinata, non per gene_id)
# --------------------------------------------------------------------------
# A differenza di annotate_regions() (che si fida del gene_id gia' presente
# nel file modkit, di provenienza/metodo sconosciuti, e forza ogni sito
# sull'UNICO trascritto canonico del gene), questa versione:
#   1) trova il gene per sola COORDINATA genomica, contro il nostro stesso
#      database costruito da Ensembl 115 (nessuna fiducia nel gene_id
#      originale del file modkit)
#   2) classifica il sito rispetto a OGNI isoforma mantenuta del gene
#      (fino a 15, vedi build_multiisoform3.py), non solo il trascritto
#      canonico -- un sito puo' quindi comparire su piu' righe, una per
#      isoforma con cui si sovrappone, con classificazioni anche diverse
#      (es. "introne" in un'isoforma, "CDS" in un'altra)

class LightMultiIsoformAnnotator:
    """
    Come la versione precedente, ma costruisce gli IntervalTree per gene e
    per feature UN CROMOSOMA ALLA VOLTA (non tutti insieme in anticipo).
    Con ~3.94M feature su tutto il genoma, tenere in memoria gli alberi di
    tutti i cromosomi contemporaneamente esaurisce la RAM disponibile; qui
    il picco di memoria resta limitato al cromosoma piu' grande, a patto
    che i siti vengano processati raggruppati per cromosoma (vedi
    annotate_regions_multi_isoform, che li ordina prima di classificarli).
    """

    def __init__(self, all_gene_spans_path: Path, transcript_spans_all_isoforms_path: Path,
                 region_features_all_isoforms_path: Path):
        import numpy as np

        # dati grezzi raggruppati per cromosoma, come ARRAY NUMPY (non liste
        # di tuple Python: per milioni di righe la differenza di memoria e'
        # enorme -- una tupla Python di 4 elementi pesa ~150-200 byte anche
        # con soli interi dentro, contro i 4-8 byte per elemento di un array
        # numpy tipizzato). NESSUN IntervalTree costruito qui.
        spans = pd.read_csv(all_gene_spans_path, sep="\t")
        self.gene_symbol = dict(zip(spans["gene_id"], spans["gene_symbol"]))
        gene_codes, self._gene_id_categories = pd.factorize(spans["gene_id"])
        self._gene_rows_by_chrom = {}
        tmp = pd.DataFrame({"chrom": spans["chrom"].values,
                             "start": spans["start"].values.astype("int32"),
                             "end": spans["end"].values.astype("int32"),
                             "gene_code": gene_codes.astype("int32")})
        for chrom, g in tmp.groupby("chrom", observed=True):
            self._gene_rows_by_chrom[chrom] = (g["start"].to_numpy(), g["end"].to_numpy(), g["gene_code"].to_numpy())
        print(f"  [multi-isoforma] {len(spans):,} geni raggruppati per cromosoma (da {all_gene_spans_path.name})")
        del spans, tmp, gene_codes

        iso_spans = pd.read_csv(transcript_spans_all_isoforms_path, sep="\t")
        self.isoforme_per_gene = {}
        self.iso_info = {}
        for row in iso_spans.itertuples(index=False):
            self.isoforme_per_gene.setdefault(row.gene_id, []).append(row.isoform_id)
            self.iso_info[row.isoform_id] = (row.start, row.end, row.strand, row.priorita, row.TSL)
        print(f"  [multi-isoforma] {len(iso_spans):,} isoforme caricate (da {transcript_spans_all_isoforms_path.name})")
        del iso_spans

        # isoform_id e region si ripetono milioni di volte in region_features
        # (isoform_id ~10-12 volte a testa, region solo 5 valori distinti):
        # tenerli come stringhe Python e' molto piu' pesante che codificarli
        # come piccoli interi + array di decodifica.
        #
        # Le feature vengono raggruppate PER ISOFORMA (non per cromosoma):
        # ogni isoforma ha in media ~12 feature, quindi ogni query confronta
        # array piccolissimi invece di scandire le ~400k feature di un intero
        # cromosoma ad ogni singolo sito -- la versione "un array per
        # cromosoma" creava e distruggeva continuamente array enormi,
        # frammentando la memoria fino ad esaurirla.
        feats = pd.read_csv(region_features_all_isoforms_path, sep="\t")
        region_codes, self._region_categories = pd.factorize(feats["region"])
        tmp = pd.DataFrame({"isoform_id": feats["isoform_id"].values,
                             "start": feats["start"].values.astype("int32"),
                             "end": feats["end"].values.astype("int32"),
                             "region_code": region_codes.astype("int8")})
        del feats, region_codes
        self._feat_arrays_by_isoform = {}
        for iso_id, g in tmp.groupby("isoform_id", observed=True, sort=False):
            self._feat_arrays_by_isoform[iso_id] = (g["start"].to_numpy(), g["end"].to_numpy(),
                                                      g["region_code"].to_numpy())
        print(f"  [multi-isoforma] {len(tmp):,} feature raggruppate per isoforma, come piccoli array numpy "
              f"(da {region_features_all_isoforms_path.name})")
        del tmp

        self._current_chrom = None
        self._current_gene_arrays = None

    def _prepara_cromosoma(self, chrom: str):
        """Seleziona (senza costruire alcuna struttura dati aggiuntiva) gli
        array numpy dei geni per il cromosoma richiesto (le feature, invece,
        sono gia' raggruppate per isoforma -- vedi __init__ -- e non
        dipendono dal cromosoma corrente)."""
        if chrom == self._current_chrom:
            return
        import gc
        gc.collect()

        self._current_chrom = chrom
        self._current_gene_arrays = self._gene_rows_by_chrom.get(chrom)

    def classify_site_multi(self, chrom: str, pos: int):
        """
        Ritorna una lista di dict, una per ogni isoforma con cui il sito si
        sovrappone: {gene_id, gene_symbol, isoform_id, region, priorita, TSL,
        n_isoforme_gene}. Se il sito non cade in nessun gene, o cade in un
        gene senza isoforme valide, o cade nel gene ma fuori da tutte le
        isoforme mantenute, ritorna una lista con UNA riga che documenta il
        caso (region = "intergenica" / "gene_senza_isoforme_valide" /
        "fuori_isoforme_mantenute").

        NOTA: per efficienza, i siti vanno passati raggruppati per
        cromosoma (vedi annotate_regions_multi_isoform).
        """
        import numpy as np
        self._prepara_cromosoma(chrom)

        if self._current_gene_arrays is None:
            return [{"gene_id": None, "gene_symbol": None, "isoform_id": None,
                      "region": "intergenica", "priorita": None, "TSL": None, "n_isoforme_gene": 0}]
        g_starts, g_ends, g_codes = self._current_gene_arrays
        gene_mask = (g_starts <= pos) & (pos <= g_ends)
        if not gene_mask.any():
            return [{"gene_id": None, "gene_symbol": None, "isoform_id": None,
                      "region": "intergenica", "priorita": None, "TSL": None, "n_isoforme_gene": 0}]

        # se piu' geni si sovrappongono alla posizione, il piu' corto (piu' specifico)
        idxs = np.nonzero(gene_mask)[0]
        lengths = g_ends[idxs] - g_starts[idxs]
        best_idx = idxs[np.argmin(lengths)]
        gene_id = self._gene_id_categories[g_codes[best_idx]]
        gene_symbol = self.gene_symbol.get(gene_id, gene_id)

        isoforme = self.isoforme_per_gene.get(gene_id, [])
        if not isoforme:
            return [{"gene_id": gene_id, "gene_symbol": gene_symbol, "isoform_id": None,
                      "region": "gene_senza_isoforme_valide", "priorita": None, "TSL": None, "n_isoforme_gene": 0}]

        # per ciascuna isoforma del gene (poche unita', al massimo 15),
        # confronto vettoriale SOLO sulle sue feature (poche decine di
        # elementi) -- non su tutte le feature del cromosoma.
        righe = []
        for iso_id in isoforme:
            iso_start, iso_end, strand, prio, tsl = self.iso_info[iso_id]
            if not (iso_start <= pos <= iso_end):
                continue
            region = "introne"
            arrays = self._feat_arrays_by_isoform.get(iso_id)
            if arrays is not None:
                f_starts, f_ends, f_region_codes = arrays
                feat_mask = (f_starts <= pos) & (pos <= f_ends)
                if feat_mask.any():
                    region = self._region_categories[f_region_codes[feat_mask][0]]
            righe.append({"gene_id": gene_id, "gene_symbol": gene_symbol, "isoform_id": iso_id,
                           "region": region, "priorita": prio, "TSL": tsl, "n_isoforme_gene": len(isoforme)})

        if not righe:
            return [{"gene_id": gene_id, "gene_symbol": gene_symbol, "isoform_id": None,
                      "region": "fuori_isoforme_mantenute", "priorita": None, "TSL": None, "n_isoforme_gene": len(isoforme)}]
        return righe


def annotate_regions_multi_isoform(df: pd.DataFrame, all_gene_spans_path: Path,
                                    transcript_spans_all_isoforms_path: Path,
                                    region_features_all_isoforms_path: Path,
                                    output_path: Path, chunk_size: int = 30_000) -> Path:
    """
    Versione "lunga" dell'annotazione: per ogni sito, una riga per ogni
    isoforma con cui si sovrappone. L'indice originale del sito viene
    preservato in 'site_idx' per poter riaggregare in seguito.

    Scrive DIRETTAMENTE su disco in blocchi (invece di accumulare tutto in
    memoria e ritornare un unico DataFrame): con ~2.6M siti e un fattore di
    espansione di 4-5x, il DataFrame finale avrebbe 12-13M righe e la lista
    di dict intermedia usata per costruirlo esaurisce la RAM disponibile.
    Ritorna il path del file scritto (CSV).
    """
    print("Caricamento tabelle multi-isoforma...")
    annot = LightMultiIsoformAnnotator(all_gene_spans_path, transcript_spans_all_isoforms_path,
                                        region_features_all_isoforms_path)
    print(f"Classificazione multi-isoforma di ogni sito, a blocchi di {chunk_size:,} "
          f"(scrittura diretta su {output_path}, puo' richiedere qualche minuto)...")

    # IMPORTANTE: ordina per cromosoma prima di classificare. Gli
    # IntervalTree vengono ricostruiti ad ogni cambio di cromosoma (vedi
    # _prepara_cromosoma): processare i siti alla rinfusa costringerebbe a
    # ricostruirli quasi ad ogni sito, vanificando l'ottimizzazione.
    df_reset = (df.reset_index(drop=True)
                  .drop(columns=["gene_id", "gene_name", "gene_type"], errors="ignore")
                  .sort_values("chrom", kind="stable")
                  .reset_index(drop=True))
    chroms = df_reset["chrom"].astype(str).to_numpy()
    ends = df_reset["end"].to_numpy()
    n = len(df_reset)

    import gc
    import resource
    first_chunk = True
    for chunk_start in range(0, n, chunk_size):
        chunk_end = min(chunk_start + chunk_size, n)
        expanded_rows = []
        site_idx_col = []
        for i in range(chunk_start, chunk_end):
            righe = annot.classify_site_multi(chroms[i], int(ends[i]))
            for r in righe:
                expanded_rows.append(r)
                site_idx_col.append(i)

        expanded = pd.DataFrame(expanded_rows)
        expanded["site_idx"] = site_idx_col
        del expanded_rows, site_idx_col
        chunk_out = expanded.merge(df_reset.iloc[chunk_start:chunk_end], left_on="site_idx",
                                    right_index=True, how="left")
        chunk_out.to_csv(output_path, sep="\t", index=False, mode="w" if first_chunk else "a",
                          header=first_chunk)
        n_righe_chunk = len(chunk_out)
        del expanded, chunk_out
        gc.collect()
        first_chunk = False
        ram_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        print(f"  {chunk_end:,}/{n:,} siti processati ({n_righe_chunk:,} righe scritte in questo blocco) "
              f"[cromosoma corrente: {chroms[chunk_end-1] if chunk_end<=n else '?'}, RAM: {ram_mb:.0f} MB]")

    print(f"Fatto. Output: {output_path}")
    return output_path


class LightCoordinateIndex:
    """
    Indice leggero (IntervalTree per cromosoma) di TUTTI i geni Ensembl
    (non solo i trascritti canonici), costruito da all_gene_spans.tsv
    (~79k righe, pochi MB) invece che dal GTF completo. Usato solo per il
    recupero per coordinate dei siti "intergenica"/"gene_non_in_gtf" -- qui
    non serve la classificazione fine in regione, solo sapere se ESISTE un
    gene che si sovrappone alla posizione.
    """

    def __init__(self, all_gene_spans_path: Path):
        from intervaltree import IntervalTree
        spans = pd.read_csv(all_gene_spans_path, sep="\t")
        self.trees = {}
        self.gene_symbol = {}
        for row in spans.itertuples(index=False):
            self.trees.setdefault(row.chrom, IntervalTree()).addi(row.start, row.end + 1, row.gene_id)
            self.gene_symbol[row.gene_id] = row.gene_symbol
        print(f"  [indice coordinate] {len(spans):,} geni indicizzati su {len(self.trees)} cromosomi "
              f"(da {all_gene_spans_path.name})")

    def lookup_by_coordinate(self, chrom: str, pos: int):
        tree = self.trees.get(str(chrom))
        if not tree:
            return []
        hits = tree.at(pos)
        return [(iv.data, self.gene_symbol.get(iv.data, iv.data)) for iv in hits]


def recover_unmapped_by_coordinate(df: pd.DataFrame, coord_index: "LightCoordinateIndex", outdir: Path):
    """
    Per i siti classificati come "intergenica" (gene_id='.' secondo modkit)
    o "gene_non_in_gtf" (gene_id presente in modkit ma assente dal nostro
    GTF), riprova un lookup INDIPENDENTE per sole coordinate genomiche
    (chrom+posizione) contro TUTTI i geni di Ensembl 115 -- non solo quello
    che modkit aveva assegnato. Se un gene Ensembl si sovrappone comunque
    alla posizione, il sito viene "recuperato" (etichettato come tale, con
    il/i gene_id trovati), anche se la classificazione fine in regione
    (5'UTR/CDS/...) non viene rifatta (richiederebbe scegliere un trascritto
    per un gene che modkit non aveva assegnato: lasciato a un'analisi
    successiva mirata, se utile).
    """
    target_mask = df["region"].isin(["intergenica", "gene_non_in_gtf"])
    sub = df[target_mask]
    print(f"  [recupero coordinate] {len(sub):,} siti da verificare (intergenica + gene_non_in_gtf)...")
    if sub.empty:
        return None

    chroms = sub["chrom"].astype(str).to_numpy()
    ends = sub["end"].to_numpy()
    found_gene_ids = [None] * len(sub)
    found_gene_symbols = [None] * len(sub)
    n_hits = 0
    for i in range(len(sub)):
        hits = coord_index.lookup_by_coordinate(chroms[i], int(ends[i]))
        if hits:
            n_hits += 1
            found_gene_ids[i] = ",".join(h[0] for h in hits)
            found_gene_symbols[i] = ",".join(h[1] for h in hits)

    result = sub.copy()
    result["gene_id_recuperato"] = found_gene_ids
    result["gene_symbol_recuperato"] = found_gene_symbols
    result["recuperato"] = result["gene_id_recuperato"].notna()

    print(f"  [recupero coordinate] {n_hits:,}/{len(sub):,} ({n_hits/len(sub):.1%}) siti RECUPERATI "
          f"(un gene Ensembl 115 si sovrappone comunque alla posizione)")

    fname = outdir / "10_siti_recuperati_per_coordinate.csv"
    result.to_csv(fname, index=False)
    print(f"  [recupero coordinate] salvato {fname}")

    # riepilogo per modifica e condizione: quanti erano "persi" e quanti recuperabili
    summary = (
        result.groupby(["mod_name", "condition", "region"], observed=True)["recuperato"]
        .agg(n_siti="size", n_recuperati="sum")
        .reset_index()
    )
    summary["pct_recuperati"] = summary["n_recuperati"] / summary["n_siti"]
    fname_sum = outdir / "10_riepilogo_recupero_per_coordinate.csv"
    summary.to_csv(fname_sum, index=False)
    print(f"  [recupero coordinate] riepilogo per mod/condizione -> {fname_sum}")

    return result


# --------------------------------------------------------------------------
# 5octies. METAGENE SU ELEMENTI ALU (per l'editing Inosina/A-to-I)
# --------------------------------------------------------------------------

class AluAnnotator:
    """
    Indice per-cromosoma (IntervalTree) degli elementi Alu (o di qualunque
    elemento ripetuto fornito in formato bed a 6-7 colonne:
    chrom,start,end,name,score,strand[,famiglia]), per calcolare la
    posizione RELATIVA (0-1, strand-aware) di un sito dentro l'elemento
    Alu che lo contiene. Utile perche' l'editing A-to-I (Inosina) avviene
    predominantemente dentro sequenze Alu, spesso con un bias posizionale
    legato alla struttura secondaria a doppio filamento che si forma tra
    Alu invertiti.
    """

    def __init__(self, bed_path: Path):
        self.trees = self._build_index(bed_path)

    @staticmethod
    def _build_index(bed_path: Path):
        from intervaltree import IntervalTree
        trees = {}
        opener = gzip.open if str(bed_path).endswith(".gz") else open
        n = 0
        with opener(bed_path, "rt") as f:
            for line in f:
                fields = line.rstrip("\n\r").split("\t")
                if len(fields) < 6:
                    continue
                chrom, start, end, name, score, strand = fields[:6]
                start, end = int(start), int(end)
                trees.setdefault(chrom, IntervalTree()).addi(start, end + 1, (start, end, strand, name))
                n += 1
        print(f"  [Alu] {n:,} elementi indicizzati su {len(trees)} cromosomi (da {bed_path.name})")
        return trees

    def relative_position(self, chrom: str, pos: int):
        """Ritorna (posizione_relativa 0-1, nome_elemento) se pos cade dentro
        un elemento Alu, altrimenti None. Se piu' elementi si sovrappongono,
        usa il piu' corto (in genere il piu' specifico)."""
        tree = self.trees.get(chrom)
        if not tree:
            return None
        hits = tree.at(pos)
        if not hits:
            return None
        iv = min(hits, key=lambda h: h.end - h.begin)
        start, end, strand, name = iv.data
        span = end - start
        if span == 0:
            return None
        rel = (pos - start) / span if strand == "+" else (end - pos) / span
        return rel, name


def compute_alu_positions(df: pd.DataFrame, alu: "AluAnnotator") -> pd.DataFrame:
    """Aggiunge `in_alu` (bool) e `alu_rel_x` (0-1, NaN se fuori Alu) al df
    (pensato per un df gia' filtrato su mod_name=='Inosina')."""
    df = df.copy()
    in_alu, rel_x = [False] * len(df), [float("nan")] * len(df)
    chroms = df["chrom"].astype(str).to_numpy()
    ends = df["end"].to_numpy()
    for i in range(len(df)):
        res = alu.relative_position(chroms[i], int(ends[i]))
        if res is not None:
            in_alu[i] = True
            rel_x[i] = res[0]
    df["in_alu"] = in_alu
    df["alu_rel_x"] = rel_x
    return df


def summary_alu_overlap(df_ino: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    """% di siti Inosina che cadono dentro un elemento Alu, per campione.
    In letteratura ci si aspetta che la stragrande maggioranza dell'editing
    A-to-I nel trascrittoma umano avvenga in Alu: un valore molto piu' basso
    qui potrebbe segnalare rumore/falsi positivi nella chiamata Inosina."""
    tab = df_ino.groupby(["condition", "replicate"]).agg(
        n_siti=("in_alu", "size"), n_in_alu=("in_alu", "sum")
    ).reset_index()
    tab["pct_in_alu"] = tab["n_in_alu"] / tab["n_siti"]
    fname = outdir / "11_inosina_overlap_alu.csv"
    tab.to_csv(fname, index=False)
    print(f"  [Alu] overlap Inosina/Alu per campione -> {fname}")
    for _, r in tab.sort_values(["condition", "replicate"]).iterrows():
        print(f"    - {r['condition']}_{r['replicate']}: {r['pct_in_alu']:.1%} dei siti Inosina in Alu "
              f"({int(r['n_in_alu']):,}/{int(r['n_siti']):,})")
    return tab


def plot_alu_metagene(df_ino: pd.DataFrame, outdir: Path, conditions=None, n_bins=50, mode="density"):
    """Metagene ristretto agli elementi Alu: posizione relativa 0-1 dentro
    l'Alu (invece che 5'UTR-CDS-3'UTR), una linea per campione."""
    sub = df_ino[df_ino["in_alu"]].copy()
    if conditions:
        sub = sub[sub["condition"].isin(conditions)]
    if sub.empty:
        print("  [Alu metagene] nessun sito Inosina in Alu da plottare")
        return
    sub["sample_label"] = sub["condition"].astype(str) + "_" + sub["replicate"].astype(str)
    bins = np.linspace(0, 1, n_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    fig, ax = plt.subplots(figsize=(9, 5.5))
    conditions_seen = sorted(sub["condition"].unique())
    cmap = plt.get_cmap("tab10")
    color_by_cond = {c: cmap(i % 10) for i, c in enumerate(conditions_seen)}
    linestyles = ["-", "--", "-.", ":"]
    for cond in conditions_seen:
        reps = sorted(sub[sub["condition"] == cond]["replicate"].unique())
        for i, rep in enumerate(reps):
            vals = sub[(sub["condition"] == cond) & (sub["replicate"] == rep)]["alu_rel_x"].values
            if len(vals) == 0:
                continue
            if mode == "density":
                y, _ = np.histogram(vals, bins=bins, density=True)
            else:
                y, _ = np.histogram(vals, bins=bins, density=False)
            ax.plot(bin_centers, y, color=color_by_cond[cond], linestyle=linestyles[i % len(linestyles)],
                     label=f"{cond}_{rep} (n={len(vals):,})", linewidth=1.8)
    ax.set_xlabel("Posizione relativa dentro l'elemento Alu (5' -> 3', strand-aware)")
    ax.set_ylabel("Densita' di siti" if mode == "density" else "N. siti")
    ax.set_title(f"Metagene su elementi Alu - Inosina ({mode})")
    ax.legend(fontsize=8)
    plt.tight_layout()
    fname = outdir / f"metagene_alu_Inosina_{mode}.png"
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  [Alu metagene] salvato {fname}")


# Categorie di `region` che rappresentano una mappatura riuscita su una
# feature del trascritto canonico, vs. categorie che rappresentano un
# fallimento/impossibilita' di mappatura (da tenere sempre distinte).
REGIONS_MAPPATE = {"5'UTR", "CDS", "esone_ncRNA", "introne", "3'UTR"}
REGIONS_NON_MAPPATE = {"intergenica", "gene_non_in_gtf", "fuori_trascritto_canonico"}

_MOTIVO_NON_MAPPATO = {
    "intergenica": "nessun gene sovrapposto secondo l'annotazione originale di modkit (gene_id='.')",
    "gene_non_in_gtf": "gene_id presente in modkit ma non trovato nel GTF Ensembl usato per l'annotazione "
                        "(probabile mismatch di versione/ID, o gene ritirato/fuso)",
    "fuori_trascritto_canonico": "gene trovato e presente nel GTF, ma la posizione cade FUORI dai confini "
                                  "genomici del trascritto canonico scelto (probabile isoforma alternativa "
                                  "con UTR piu' lunghe/diverse rispetto al trascritto MANE/canonico)",
}


def summary_mapping_status(df: pd.DataFrame) -> pd.DataFrame:
    """
    Quanti siti sono stati mappati con successo su una regione del
    trascritto canonico (5'UTR/CDS/esone_ncRNA/introne/3'UTR) vs. quanti no,
    per modifica e per la causa esatta del mancato mapping. Utile per capire
    quanto ci si puo' fidare delle analisi metagene/regione (che usano solo
    i siti mappati) e se conviene rivedere l'annotazione per una modifica
    specifica.
    """
    if "region" not in df.columns:
        return pd.DataFrame()
    tab = df.groupby(["mod_name", "region"], observed=True).size().reset_index(name="n_siti")
    tab["stato"] = tab["region"].apply(lambda r: "mappato" if r in REGIONS_MAPPATE else "NON mappato")
    tab["motivo"] = tab["region"].map(_MOTIVO_NON_MAPPATO).fillna("mappato con successo su una feature del trascritto")
    tab["pct_sul_totale_modifica"] = tab.groupby("mod_name")["n_siti"].transform(lambda s: s / s.sum())
    tab = tab.sort_values(["mod_name", "n_siti"], ascending=[True, False])

    # riepilogo aggregato mappato/non mappato per modifica, stampato a schermo
    agg = tab.groupby(["mod_name", "stato"])["n_siti"].sum().reset_index()
    agg["pct"] = agg.groupby("mod_name")["n_siti"].transform(lambda s: s / s.sum())
    print("  [mappatura] riepilogo mappato vs non mappato per modifica:")
    for mod in agg["mod_name"].unique():
        sub = agg[agg["mod_name"] == mod]
        parts = [f"{r['stato']}={r['pct']:.1%}" for _, r in sub.iterrows()]
        print(f"    - {mod}: " + ", ".join(parts))

    return tab


def plot_mapping_status(tab: pd.DataFrame, outdir: Path):
    if tab.empty:
        return
    pivot = tab.pivot(index="mod_name", columns="region", values="pct_sul_totale_modifica").fillna(0)
    region_order = ["5'UTR", "CDS", "esone_ncRNA", "introne", "3'UTR",
                     "fuori_trascritto_canonico", "gene_non_in_gtf", "intergenica"]
    cols = [c for c in region_order if c in pivot.columns] + [c for c in pivot.columns if c not in region_order]
    pivot = pivot[cols]
    # colori: verdi per le categorie "mappate", rossi/grigi per le "non mappate"
    color_map = {
        "5'UTR": "#2166ac", "CDS": "#4393c3", "esone_ncRNA": "#92c5de", "introne": "#d1e5f0", "3'UTR": "#053061",
        "fuori_trascritto_canonico": "#f4a582", "gene_non_in_gtf": "#d6604d", "intergenica": "#b2182b",
    }
    colors = [color_map.get(c, "#999999") for c in pivot.columns]
    ax = pivot.plot(kind="barh", stacked=True, figsize=(10, 0.6 * len(pivot) + 2), color=colors)
    ax.set_xlabel("Frazione dei siti")
    ax.set_title("Stato di mappatura sul trascritto canonico, per modifica\n(blu = mappato con successo, rosso/arancio = NON mappato)")
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    fname = outdir / "stato_mappatura_per_modifica.png"
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  [mappatura] salvato {fname}")


# --------------------------------------------------------------------------
# 5bis. METAGENE PLOT (posizione relativa sul trascritto maturo: 5'UTR/CDS/3'UTR)
# --------------------------------------------------------------------------

# mappa il nome di regione (colonna `region`) al featuretype GTF corrispondente
REGION_TO_FEATURETYPE = {"5'UTR": "5'UTR", "CDS": "CDS", "3'UTR": "3'UTR"}
# (nella versione leggera le region_features.tsv usano gia' le etichette
# "coarse" 5'UTR/CDS/3'UTR invece dei nomi grezzi del GTF five_prime_utr/
# CDS/three_prime_utr: qui la mappa e' quindi un'identita').
# posizione di partenza di ciascun segmento sull'asse del metagene (0-1, 1-2, 2-3)
REGION_TO_SEGMENT_BASE = {"5'UTR": 0, "CDS": 1, "3'UTR": 2}


def _relative_position_in_segment(tx_info: dict, pos: int, featuretype: str):
    """
    Calcola la posizione relativa (0-1) di `pos` all'interno del segmento
    (5'UTR, CDS o 3'UTR) del trascritto, seguendo le coordinate "spliced"
    (cioe' sul trascritto maturo, saltando gli introni) e rispettando lo
    strand. Ritorna None se il sito non ricade in nessun intervallo del
    segmento richiesto.
    """
    intervals = sorted(
        [(s, e) for (s, e, t) in tx_info["feats"] if t == featuretype],
        key=lambda iv: iv[0],
    )
    if not intervals:
        return None

    total_len = sum(e - s + 1 for s, e in intervals)
    if total_len == 0:
        return None

    ordered = intervals if tx_info["strand"] == "+" else intervals[::-1]

    offset = 0
    for s, e in ordered:
        if s <= pos <= e:
            offset += (pos - s) if tx_info["strand"] == "+" else (e - pos)
            return offset / total_len
        offset += (e - s + 1)
    return None


def compute_metagene_positions(df: pd.DataFrame, annot: "LightRegionAnnotator") -> pd.DataFrame:
    """
    Aggiunge la colonna `metagene_x` (0-3: 0-1 = 5'UTR, 1-2 = CDS, 2-3 = 3'UTR)
    solo per i siti che ricadono in 5'UTR/CDS/3'UTR di un trascritto codificante.
    I siti introne/esone_ncRNA/intergenica/gene_non_in_gtf non hanno un
    metagene definito (restano NaN) e vanno esclusi dal plot.
    """
    xs = []
    for region, end, tx_id in zip(df["region"].astype(str).values, df["end"].values, df["transcript_used"].astype(str).values):
        featuretype = REGION_TO_FEATURETYPE.get(region)
        if featuretype is None or tx_id is None or tx_id == "None":
            xs.append(float("nan"))
            continue
        tx_info = annot.tx_index.get(tx_id)
        if tx_info is None:
            xs.append(float("nan"))
            continue
        rel = _relative_position_in_segment(tx_info, int(end), featuretype)
        if rel is None:
            xs.append(float("nan"))
        else:
            xs.append(REGION_TO_SEGMENT_BASE[region] + rel)
    df["metagene_x"] = xs
    return df


def plot_metagene(df: pd.DataFrame, mod_name: str, outdir: Path, conditions=None, n_bins=60, mode="density", subset_label=None):
    """
    Metagene plot (5'UTR-CDS-3'UTR normalizzati) per una modifica, con una
    linea per ciascun campione (condizione_replicato).

    mode:
      - "density": densita' normalizzata (la forma della distribuzione,
        confrontabile anche tra campioni con n. siti molto diverso)
      - "counts":  conteggio assoluto di siti per bin (mostra anche le
        differenze di ampiezza/numerosita' tra campioni)
      - "pct_mod": % di modifica MEDIA (fraction_modified_pct) dei siti
        in ciascun bin, invece del numero di siti

    subset_label: se fornito, viene aggiunto al titolo e al nome del file
    (es. "solo geni ZNF/ZFP") per generare metagene ristretti a un
    sottoinsieme di geni senza sovrascrivere il metagene generale.
    """
    sub = df[(df["mod_name"] == mod_name) & df["metagene_x"].notna()].copy()
    if conditions:
        sub = sub[sub["condition"].isin(conditions)]
    if sub.empty:
        print(f"  [metagene] nessun sito con posizione definita per {mod_name}, salto il plot ({mode})")
        return

    sub["sample_label"] = sub["condition"] + "_" + sub["replicate"]
    bins = np.linspace(0, 3, n_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2

    fig, ax = plt.subplots(figsize=(10, 5.5))
    conditions_seen = sorted(sub["condition"].unique())
    cmap = plt.get_cmap("tab10")
    color_by_cond = {c: cmap(i % 10) for i, c in enumerate(conditions_seen)}
    linestyles = ["-", "--", "-.", ":"]

    for cond in conditions_seen:
        reps = sorted(sub[sub["condition"] == cond]["replicate"].unique())
        for i, rep in enumerate(reps):
            sample = sub[(sub["condition"] == cond) & (sub["replicate"] == rep)]
            vals = sample["metagene_x"].values
            if len(vals) == 0:
                continue

            if mode == "density":
                y, _ = np.histogram(vals, bins=bins, density=True)
            elif mode == "counts":
                y, _ = np.histogram(vals, bins=bins, density=False)
            elif mode == "pct_mod":
                bin_idx = np.digitize(vals, bins) - 1
                bin_idx = np.clip(bin_idx, 0, n_bins - 1)
                pct_vals = sample["fraction_modified_pct"].values
                y = np.full(n_bins, np.nan)
                for b in range(n_bins):
                    mask = bin_idx == b
                    if mask.any():
                        y[b] = pct_vals[mask].mean()
            else:
                raise ValueError(f"mode sconosciuta: {mode}")

            ax.plot(
                bin_centers, y,
                color=color_by_cond[cond],
                linestyle=linestyles[i % len(linestyles)],
                label=f"{cond}_{rep} (n={len(vals):,})",
                linewidth=1.8,
            )

    ax.axvline(1, color="grey", linestyle=":", linewidth=1)
    ax.axvline(2, color="grey", linestyle=":", linewidth=1)
    ax.set_xlim(0, 3)
    ax.set_xticks([0.5, 1.5, 2.5])
    ax.set_xticklabels(["5'UTR", "CDS", "3'UTR"])

    ylabel = {
        "density": "Densita' di siti (normalizzata per campione)",
        "counts": "N. siti (conteggio assoluto per bin)",
        "pct_mod": "% modifica media dei siti nel bin",
    }[mode]
    suffix = {"density": "densita", "counts": "conteggi", "pct_mod": "pct_modifica"}[mode]
    title_suffix = f" — {subset_label}" if subset_label else ""
    fname_suffix = f"_{subset_label.replace('/', '-').replace(' ', '_')}" if subset_label else ""

    ax.set_ylabel(ylabel)
    ax.set_title(f"Metagene plot - {mod_name} ({ylabel.split('(')[0].strip()}){title_suffix}")
    ax.legend(fontsize=8, loc="upper right")
    plt.tight_layout()
    fname = outdir / f"metagene_{mod_name}_{suffix}{fname_suffix}.png"
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  [metagene] salvato {fname}")


# --------------------------------------------------------------------------
# 5ter. TEST STATISTICO SU UNA FINESTRA DEL METAGENE (es. attorno allo stop codon)
# --------------------------------------------------------------------------

def test_metagene_window(
    df: pd.DataFrame, mod_name: str, outdir: Path,
    x_min: float, x_max: float, cond_a: str, cond_b: str,
):
    """
    Confronta, con un test di Mann-Whitney U (non parametrico, non assume
    normalita'), la % di modifica (fraction_modified_pct) dei siti che
    ricadono nella finestra [x_min, x_max] del metagene, tra due condizioni.
    Aggrega i replicati di ciascuna condizione. Salva un CSV riassuntivo.
    """
    from scipy.stats import mannwhitneyu

    sub = df[(df["mod_name"] == mod_name) & df["metagene_x"].between(x_min, x_max)]

    a = sub[sub["condition"] == cond_a]["fraction_modified_pct"].dropna().values
    b = sub[sub["condition"] == cond_b]["fraction_modified_pct"].dropna().values

    if len(a) == 0 or len(b) == 0:
        print(f"  [test] finestra [{x_min},{x_max}] vuota per {cond_a} o {cond_b}, salto il test")
        return None

    stat, pval = mannwhitneyu(a, b, alternative="two-sided")
    # effect size: rank-biserial correlation (r = 1 - 2U/(n1*n2))
    r_rb = 1 - (2 * stat) / (len(a) * len(b))

    result = pd.DataFrame([{
        "mod_name": mod_name,
        "finestra_metagene": f"[{x_min},{x_max}]",
        "condizione_A": cond_a, "n_A": len(a), "pct_mod_mediana_A": float(np.median(a)),
        "condizione_B": cond_b, "n_B": len(b), "pct_mod_mediana_B": float(np.median(b)),
        "mannwhitney_U": float(stat), "p_value": float(pval),
        "rank_biserial_r": float(r_rb),
    }])
    fname = outdir / f"test_{mod_name}_{cond_a}_vs_{cond_b}_finestra_{x_min}-{x_max}.csv"
    result.to_csv(fname, index=False)
    print(f"  [test] {mod_name} {cond_a} (mediana={np.median(a):.2f}%, n={len(a)}) vs "
          f"{cond_b} (mediana={np.median(b):.2f}%, n={len(b)}): p={pval:.3g}, r={r_rb:.3f} -> {fname}")
    return result


# --------------------------------------------------------------------------
# 5quater. CONFRONTO SITO-PER-SITO (PAIRED) TRA DUE CONDIZIONI
# --------------------------------------------------------------------------

def paired_site_comparison(df: pd.DataFrame, mod_name: str, cond_a: str, cond_b: str, outdir: Path):
    """
    Confronto sito-per-sito (stesso chrom+start+end+gene_id) tra due
    condizioni per una data modifica. I replicati di ciascuna condizione
    vengono prima mediati per ottenere un unico valore di % di modifica
    per sito e condizione; vengono tenuti solo i siti presenti (con
    copertura, quindi chiamati da modkit) in ENTRAMBE le condizioni.

    Calcola:
      - delta = % modifica (cond_b) - % modifica (cond_a) per ogni sito
      - test di Wilcoxon signed-rank (paired, non parametrico) su tutti i
        siti appaiati: e' un test piu' potente del Mann-Whitney "a finestra"
        perche' sfrutta l'accoppiamento sito-per-sito invece di confrontare
        due popolazioni indipendenti.

    Salva un CSV con tutti i siti appaiati (ordinato per delta, cosi' i
    geni con la maggiore perdita/guadagno di modifica sono in cima/fondo)
    e due grafici: uno scatter (esagoni, per gestire milioni di punti) e
    un istogramma dei delta.
    """
    from scipy.stats import wilcoxon

    site_cols = ["chrom", "start", "end", "gene_id"]
    sub = df[(df["mod_name"] == mod_name) & df["condition"].isin([cond_a, cond_b])]
    if sub.empty:
        print(f"  [paired] nessun dato per {mod_name} ({cond_a}/{cond_b})")
        return None

    per_site_cond = (
        sub.groupby(site_cols + ["condition"], observed=True)
        .agg(pct_mod=("fraction_modified_pct", "mean"), cov=("N_valid_cov", "mean"))
        .reset_index()
    )
    wide = per_site_cond.pivot(index=site_cols, columns="condition", values=["pct_mod", "cov"])
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()

    col_a, col_b = f"pct_mod_{cond_a}", f"pct_mod_{cond_b}"
    if col_a not in wide.columns or col_b not in wide.columns:
        print(f"  [paired] condizioni {cond_a}/{cond_b} non trovate per {mod_name}")
        return None

    paired = wide.dropna(subset=[col_a, col_b]).copy()
    if paired.empty:
        print(f"  [paired] nessun sito comune tra {cond_a} e {cond_b} per {mod_name}")
        return None

    paired["delta"] = paired[col_b] - paired[col_a]

    # aggiunge gene_symbol/gene_name/gene_type/region per leggibilita' (dal trascritto canonico)
    extra_cols = [c for c in ["gene_symbol", "gene_name", "gene_type", "region", "rna_class"] if c in df.columns]
    if extra_cols:
        lookup = df.drop_duplicates(subset=["gene_id"])[["gene_id"] + extra_cols]
        paired = paired.merge(lookup, on="gene_id", how="left")

    stat, pval = wilcoxon(
        paired[col_a].to_numpy(dtype=float), paired[col_b].to_numpy(dtype=float), zero_method="wilcox"
    )
    n = len(paired)
    median_delta = float(paired["delta"].median())
    frac_down = float((paired["delta"] < 0).mean())
    frac_up = float((paired["delta"] > 0).mean())

    print(f"  [paired] {mod_name}: {n:,} siti appaiati {cond_a} vs {cond_b} | "
          f"delta mediano={median_delta:+.2f}pp | "
          f"{frac_down:.1%} in calo / {frac_up:.1%} in aumento in {cond_b} | "
          f"Wilcoxon p={pval:.3g}")

    fname = outdir / f"siti_appaiati_{mod_name}_{cond_a}_vs_{cond_b}.csv"
    paired.sort_values("delta").to_csv(fname, index=False)
    print(f"  [paired] salvato {fname}")

    stats_row = pd.DataFrame([{
        "mod_name": mod_name, "condizione_A": cond_a, "condizione_B": cond_b,
        "n_siti_appaiati": n, "delta_mediano_pp": median_delta,
        "frazione_siti_in_calo": frac_down, "frazione_siti_in_aumento": frac_up,
        "wilcoxon_stat": float(stat), "p_value": float(pval),
    }])
    stats_row.to_csv(outdir / f"siti_appaiati_{mod_name}_{cond_a}_vs_{cond_b}_stats.csv", index=False)

    _plot_paired_scatter(paired, col_a, col_b, mod_name, cond_a, cond_b, outdir)
    _plot_paired_delta_hist(paired, mod_name, cond_a, cond_b, median_delta, pval, outdir)
    return paired


def _plot_paired_scatter(paired, col_a, col_b, mod_name, cond_a, cond_b, outdir):
    fig, ax = plt.subplots(figsize=(6.5, 6))
    hb = ax.hexbin(paired[col_a], paired[col_b], gridsize=60, cmap="viridis", mincnt=1, bins="log")
    ax.plot([0, 100], [0, 100], color="red", linestyle="--", linewidth=1, label="y = x (nessuna differenza)")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_xlabel(f"% modifica in {cond_a}")
    ax.set_ylabel(f"% modifica in {cond_b}")
    ax.set_title(f"{mod_name}: confronto sito-per-sito ({len(paired):,} siti)")
    ax.legend(fontsize=8, loc="upper left")
    fig.colorbar(hb, ax=ax, label="N. siti (log)")
    plt.tight_layout()
    fname = outdir / f"scatter_appaiato_{mod_name}_{cond_a}_vs_{cond_b}.png"
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  [paired] salvato {fname}")


def _plot_paired_delta_hist(paired, mod_name, cond_a, cond_b, median_delta, pval, outdir):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(paired["delta"], bins=100, color="#4C72B0", edgecolor="none")
    ax.axvline(0, color="black", linewidth=1)
    ax.axvline(median_delta, color="red", linestyle="--", linewidth=1.5,
               label=f"mediana = {median_delta:+.2f}pp")
    ax.set_xlabel(f"delta % modifica ({cond_b} - {cond_a})")
    ax.set_ylabel("N. siti")
    ax.set_title(f"{mod_name}: distribuzione dei delta sito-per-sito ({cond_b} vs {cond_a})\n"
                 f"Wilcoxon p={pval:.3g}, n={len(paired):,}")
    ax.legend(fontsize=9)
    plt.tight_layout()
    fname = outdir / f"delta_hist_{mod_name}_{cond_a}_vs_{cond_b}.png"
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  [paired] salvato {fname}")


# --------------------------------------------------------------------------
# 5quinquies. CLASSIFICAZIONE IN CLASSI DI RNA (mRNA/tRNA/rRNA/non-coding...)
# --------------------------------------------------------------------------

# Mappa il gene_type (biotipo Ensembl/GENCODE) a una classe di RNA piu'
# leggibile/aggregata. Il biotipo "." indica siti senza gene sovrapposto.
GENE_TYPE_TO_RNA_CLASS = {
    ".": "nessun gene (intergenico)",
    "protein_coding": "mRNA",
    "lncRNA": "lncRNA",
    "miRNA": "miRNA",
    "snoRNA": "snoRNA",
    "snRNA": "snRNA",
    "misc_RNA": "misc_RNA",
    "rRNA": "rRNA",
    "rRNA_pseudogene": "rRNA (pseudogene)",
    "Mt_rRNA": "rRNA mitocondriale",
    "Mt_tRNA": "tRNA mitocondriale",
    "TEC": "TEC (da confermare)",
    "artifact": "artefatto tecnico",
}


def classify_rna_class(gene_type) -> str:
    if gene_type is None or pd.isna(gene_type):
        return "nessun gene (intergenico)"
    gene_type = str(gene_type)
    if gene_type in GENE_TYPE_TO_RNA_CLASS:
        return GENE_TYPE_TO_RNA_CLASS[gene_type]
    if "pseudogene" in gene_type:
        return "pseudogene"
    if gene_type.startswith("IG_"):
        return "gene immunoglobulinico"
    if gene_type.startswith("TR_"):
        return "gene recettore-T"
    return "altro non-coding"


def add_rna_class(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["rna_class"] = pd.Categorical(df["gene_type"].astype(str).map(classify_rna_class))
    return df


def summary_by_rna_class(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["mod_name", "rna_class"], observed=True)
        .size()
        .reset_index(name="n_siti")
        .sort_values(["mod_name", "n_siti"], ascending=[True, False])
    )


def plot_rna_class(summary_rna: pd.DataFrame, outdir: Path):
    if summary_rna.empty:
        return
    pivot = summary_rna.pivot(index="rna_class", columns="mod_name", values="n_siti").fillna(0)
    pivot = pivot.loc[pivot.sum(axis=1).sort_values(ascending=False).index]
    ax = pivot.plot(kind="bar", figsize=(11, 5.5), logy=True)
    ax.set_ylabel("N. siti (scala log)")
    ax.set_xlabel("Classe di RNA (da gene_type)")
    ax.set_title("Distribuzione dei siti modificati per classe di RNA")
    plt.xticks(rotation=40, ha="right")
    plt.tight_layout()
    fname = outdir / "distribuzione_classi_rna.png"
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  [rna_class] salvato {fname}")


# --------------------------------------------------------------------------
# 5sexies. ARRICCHIMENTO DEL SOTTOINSIEME "RESPONSIVO" (delta forte tra due condizioni)
# --------------------------------------------------------------------------

def _benjamini_hochberg(pvals):
    """Correzione BH per test multipli. Ritorna i p-value corretti (q-value)."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    # rende la sequenza monotona (dal basso verso l'alto, prendendo il minimo cumulato da destra)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    qvals = np.empty(n)
    qvals[order] = np.clip(ranked, 0, 1)
    return qvals


def enrich_responsive_subset(
    paired: pd.DataFrame, mod_name: str, outdir: Path, delta_threshold: float = -30.0,
    category_cols=("region", "gene_type", "rna_class"),
):
    """
    A partire dalla tabella di siti appaiati (colonna `delta` = cond_B - cond_A),
    isola il sottoinsieme "responsivo" (delta < delta_threshold, cioe' una forte
    perdita di modifica) e testa, per ciascuna categoria (regione genica,
    gene_type, ...), se e' sovra/sotto-rappresentata rispetto al resto dei
    siti testati, con un test esatto di Fisher (2x2) + correzione BH.

    Salva:
      - CSV con i soli siti responsivi (per uso downstream, es. liste geniche)
      - CSV con i risultati dell'arricchimento per ciascuna categoria testata
      - CSV con la classifica dei geni per frazione di siti responsivi
    """
    from scipy.stats import fisher_exact

    responsive = paired[paired["delta"] < delta_threshold].copy()
    background = paired[paired["delta"] >= delta_threshold]
    n_resp, n_bg = len(responsive), len(background)
    print(f"  [arricchimento] {mod_name}: {n_resp:,} siti responsivi (delta<{delta_threshold}) "
          f"su {len(paired):,} totali ({n_resp/len(paired):.1%})")

    responsive.to_csv(outdir / f"sottoinsieme_responsivo_{mod_name}_delta_lt_{delta_threshold}.csv", index=False)

    enrich_rows = []
    for col in category_cols:
        if col not in paired.columns:
            continue
        for cat in paired[col].dropna().unique():
            a = (responsive[col] == cat).sum()               # responsivi in categoria
            b = n_resp - a                                     # responsivi fuori categoria
            c = (background[col] == cat).sum()                 # background in categoria
            d = n_bg - c                                        # background fuori categoria
            if a + c < 5:  # categoria troppo rara per un test sensato
                continue
            odds_ratio, pval = fisher_exact([[a, b], [c, d]], alternative="two-sided")
            enrich_rows.append({
                "categoria_tipo": col, "categoria_valore": cat,
                "n_responsivi_in_categoria": a, "n_responsivi_totali": n_resp,
                "n_background_in_categoria": c, "n_background_totali": n_bg,
                "pct_responsivi": a / n_resp if n_resp else np.nan,
                "pct_background": c / n_bg if n_bg else np.nan,
                "odds_ratio": odds_ratio, "p_value": pval,
            })

    enrich_df = pd.DataFrame(enrich_rows)
    if not enrich_df.empty:
        enrich_df["q_value_BH"] = _benjamini_hochberg(enrich_df["p_value"].values)
        enrich_df = enrich_df.sort_values("p_value")
        fname = outdir / f"arricchimento_{mod_name}_delta_lt_{delta_threshold}.csv"
        enrich_df.to_csv(fname, index=False)
        print(f"  [arricchimento] salvato {fname}")
        sig = enrich_df[enrich_df["q_value_BH"] < 0.05]
        if not sig.empty:
            print(f"  [arricchimento] categorie significative (q<0.05):")
            for _, r in sig.iterrows():
                direction = "sovra" if r["pct_responsivi"] > r["pct_background"] else "sotto"
                print(f"    - {r['categoria_tipo']}={r['categoria_valore']}: {direction}-rappresentata "
                      f"({r['pct_responsivi']:.1%} vs {r['pct_background']:.1%}), q={r['q_value_BH']:.2e}")

    # classifica geni per frazione di siti responsivi (solo geni con >=3 siti testati)
    gene_col = "gene_symbol" if "gene_symbol" in paired.columns else "gene_id"
    gene_stats = (
        paired.groupby(["gene_id", gene_col])
        .agg(n_siti_totali=("delta", "size"), n_siti_responsivi=("delta", lambda s: (s < delta_threshold).sum()),
             delta_medio=("delta", "mean"))
        .reset_index()
    )
    gene_stats = gene_stats[gene_stats["n_siti_totali"] >= 3].copy()
    gene_stats["frazione_responsiva"] = gene_stats["n_siti_responsivi"] / gene_stats["n_siti_totali"]
    gene_stats = gene_stats.sort_values(["frazione_responsiva", "n_siti_responsivi"], ascending=False)
    fname_genes = outdir / f"geni_classificati_{mod_name}_delta_lt_{delta_threshold}.csv"
    gene_stats.to_csv(fname_genes, index=False)
    print(f"  [arricchimento] classifica geni salvata in {fname_genes} "
          f"({len(gene_stats):,} geni con >=3 siti testati)")

    return responsive, enrich_df, gene_stats


# --------------------------------------------------------------------------
# 5sexies-ter. ARRICCHIMENTO DI FAMIGLIE GENICHE (regex su gene_symbol) + REGIONE PER SOTTOINSIEME
# --------------------------------------------------------------------------

def family_enrichment(gene_stats: pd.DataFrame, family_name: str, pattern: str, outdir: Path,
                       mod_name: str, delta_threshold: float, responsive_col="n_siti_responsivi"):
    """
    Test di Fisher esatto: i geni la cui sigla (gene_symbol) fa match con
    `pattern` (es. famiglia zinc-finger ZNF/ZFP) sono sovra-rappresentati
    tra i geni "responsivi" (>=1 sito con delta < soglia) rispetto al
    resto dei geni testati?
    """
    from scipy.stats import fisher_exact

    gene_col = "gene_symbol" if "gene_symbol" in gene_stats.columns else "gene_id"
    is_family = gene_stats[gene_col].astype(str).str.match(pattern, case=False)
    is_responsive = gene_stats[responsive_col] >= 1

    a = int((is_family & is_responsive).sum())
    b = int((is_family & ~is_responsive).sum())
    c = int((~is_family & is_responsive).sum())
    d = int((~is_family & ~is_responsive).sum())
    odds, pval = fisher_exact([[a, b], [c, d]], alternative="two-sided")

    result = pd.DataFrame([{
        "famiglia": family_name, "pattern_regex": pattern, "modifica": mod_name,
        "n_famiglia_responsivi": a, "n_famiglia_totali": a + b,
        "pct_famiglia_responsiva": a / (a + b) if (a + b) else np.nan,
        "n_altri_responsivi": c, "n_altri_totali": c + d,
        "pct_altri_responsiva": c / (c + d) if (c + d) else np.nan,
        "odds_ratio": odds, "p_value": pval,
    }])
    fname = outdir / f"arricchimento_famiglia_{family_name.split('/')[0]}_{mod_name}_delta_lt_{delta_threshold}.csv"
    result.to_csv(fname, index=False)
    print(f"  [famiglia genica] {family_name}: {a}/{a+b} ({a/(a+b) if a+b else 0:.1%}) responsivi vs "
          f"{c}/{c+d} ({c/(c+d) if c+d else 0:.1%}) nel resto -> odds ratio={odds:.2f}, p={pval:.3e} -> {fname}")
    return result


def region_distribution_subset(df: pd.DataFrame, mod_name: str, gene_symbol_pattern: str, outdir: Path,
                                subset_label: str, conditions=None):
    """
    Distribuzione dei siti per regione genica (5'UTR/CDS/3'UTR/introne/...),
    ristretta ai soli geni il cui gene_symbol fa match con
    `gene_symbol_pattern` (es. famiglia ZNF/ZFP), per condizione.
    """
    if "region" not in df.columns or "gene_symbol" not in df.columns:
        print("  [regione sottoinsieme] servono le colonne region/gene_symbol (--gtf), salto")
        return None

    sub = df[(df["mod_name"] == mod_name) & df["gene_symbol"].astype(str).str.match(gene_symbol_pattern, case=False)]
    if conditions:
        sub = sub[sub["condition"].isin(conditions)]
    if sub.empty:
        print(f"  [regione sottoinsieme] nessun sito {mod_name} per {subset_label}")
        return None

    tab = (
        sub.groupby(["condition", "region"], observed=True).size()
        .reset_index(name="n_siti")
    )
    tab_pct = tab.copy()
    tab_pct["pct_sul_totale_condizione"] = tab_pct.groupby("condition")["n_siti"].transform(lambda s: s / s.sum())

    fname = outdir / f"regione_{mod_name}_{subset_label.replace('/', '-').replace(' ', '_')}.csv"
    tab_pct.to_csv(fname, index=False)
    print(f"  [regione sottoinsieme] {subset_label}, {mod_name}: {sub['gene_symbol'].nunique()} geni, "
          f"{len(sub)} siti totali -> {fname}")

    pivot = tab_pct.pivot(index="condition", columns="region", values="pct_sul_totale_condizione").fillna(0)
    region_order = ["5'UTR", "CDS", "esone_ncRNA", "introne", "3'UTR", "intergenica", "gene_non_in_gtf"]
    cols = [c for c in region_order if c in pivot.columns] + [c for c in pivot.columns if c not in region_order]
    pivot = pivot[cols]
    ax = pivot.plot(kind="bar", stacked=True, figsize=(8, 5))
    ax.set_ylabel("Frazione dei siti (per condizione)")
    ax.set_title(f"Distribuzione per regione — {mod_name}, {subset_label}")
    plt.xticks(rotation=0)
    plt.legend(title="Regione", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    fname_png = outdir / f"regione_{mod_name}_{subset_label.replace('/', '-').replace(' ', '_')}.png"
    plt.savefig(fname_png, dpi=150)
    plt.close()
    print(f"  [regione sottoinsieme] salvato {fname_png}")
    return tab_pct


# --------------------------------------------------------------------------
# 5sexies-bis. ARRICCHIMENTO GO (Biological Process / Molecular Function / Cellular Component)
# --------------------------------------------------------------------------

def load_gmt(gmt_path: Path) -> dict:
    """Legge un file .gmt (formato MSigDB/Enrichr: termine\\turl_o_descrizione\\tgene1\\tgene2...)."""
    sets = {}
    with open(gmt_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term = parts[0]
            genes = set(g.upper() for g in parts[2:] if g)
            if genes:
                sets[term] = genes
    return sets


def go_enrichment(
    foreground_symbols, background_symbols, gmt_path: Path,
    min_term_size: int = 5, max_term_size: int = 1000,
):
    """
    Test ipergeometrico (equivalente a Fisher esatto 2x2 unilaterale) per
    ciascun termine GO del gmt: e' sovra-rappresentato nel foreground
    (es. i geni "responsivi") rispetto al background (es. tutti i geni
    testati)? Corregge per test multipli con Benjamini-Hochberg.

    foreground_symbols/background_symbols: iterabili di gene symbol
    (il foreground DEVE essere un sottoinsieme del background).
    """
    from scipy.stats import hypergeom

    gene_sets = load_gmt(gmt_path)
    fg = set(s.upper() for s in foreground_symbols)
    bg = set(s.upper() for s in background_symbols) | fg  # il foreground e' sempre incluso nel background
    n_fg, n_bg = len(fg), len(bg)

    rows = []
    for term, genes in gene_sets.items():
        term_in_bg = genes & bg
        if not (min_term_size <= len(term_in_bg) <= max_term_size):
            continue
        term_in_fg = genes & fg
        k = len(term_in_fg)
        if k == 0:
            continue
        M, n, N = n_bg, len(term_in_bg), n_fg
        # P(X >= k) con X ~ Ipergeometrica(M, n, N)
        pval = hypergeom.sf(k - 1, M, n, N)
        rows.append({
            "termine_GO": term, "n_geni_foreground_nel_termine": k,
            "n_geni_totali_nel_termine": n, "n_foreground": N, "n_background": M,
            "fold_enrichment": (k / N) / (n / M) if N and M else np.nan,
            "geni": ",".join(sorted(term_in_fg)),
            "p_value": pval,
        })

    result = pd.DataFrame(rows)
    if not result.empty:
        result["q_value_BH"] = _benjamini_hochberg(result["p_value"].values)
        result = result.sort_values("p_value")
    return result


# --------------------------------------------------------------------------
# 5septies. QUERY SU GENI SPECIFICI (es. enzimi di interesse) + PATHWAY MIRROR
# --------------------------------------------------------------------------

GO_GMT_URLS = {
    "c5.bp.v7.0.symbols.gmt": "https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/msigdb/c5.bp.v7.0.symbols.gmt",
    "c5.mf.v7.0.symbols.gmt": "https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/msigdb/c5.mf.v7.0.symbols.gmt",
    "c5.cc.v7.0.symbols.gmt": "https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/msigdb/c5.cc.v7.0.symbols.gmt",
}

# Set genici aggiuntivi (stesso mirror): KEGG, Reactome, Hallmark, e bersagli
# noti dei fattori di trascrizione (TFT, TRANSFAC-based). Usati per capire
# in quali pathway/processi ricadono i geni responsivi, e quali geni sono
# regolati dai TF piu' responsivi.
PATHWAY_GMT_URLS = {
    "c2.cp.kegg.v7.0.symbols.gmt": "https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/msigdb/c2.cp.kegg.v7.0.symbols.gmt",
    "c2.cp.reactome.v7.0.symbols.gmt": "https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/msigdb/c2.cp.reactome.v7.0.symbols.gmt",
    "h.all.v7.0.symbols.gmt": "https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/msigdb/h.all.v7.0.symbols.gmt",
    "c3.tft.v7.0.symbols.gmt": "https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/msigdb/c3.tft.v7.0.symbols.gmt",
}


def ensure_gmts(go_dir: Path, urls: dict) -> bool:
    """Scarica (se assenti) i .gmt indicati in `urls` dal mirror GitHub del
    Broad Institute. Richiede una connessione verso raw.githubusercontent.com;
    se non disponibile, ritorna False e l'arricchimento relativo viene saltato."""
    import urllib.request

    go_dir.mkdir(parents=True, exist_ok=True)
    ok = True
    for fname, url in urls.items():
        fpath = go_dir / fname
        if fpath.exists() and fpath.stat().st_size > 0:
            continue
        try:
            print(f"  [gmt] scarico {fname}...")
            urllib.request.urlretrieve(url, fpath)
        except Exception as e:
            print(f"  [gmt] impossibile scaricare {fname} ({e}); relativo arricchimento saltato")
            ok = False
    return ok


def ensure_go_gmts(go_dir: Path) -> bool:
    """Retrocompatibilita': scarica solo i 3 file GO BP/MF/CC."""
    return ensure_gmts(go_dir, GO_GMT_URLS)


def run_gene_set_enrichment(foreground_symbols, background_symbols, go_dir: Path, outdir: Path,
                             label: str, files: dict, sig_threshold: float = 0.05):
    """Versione generica di run_go_enrichment_bp_mf_cc: esegue go_enrichment
    per ciascuna voce di `files` ({etichetta: nome_file.gmt}) e salva un CSV
    per ciascuna. Usata sia per GO (BP/MF/CC) sia per pathway (KEGG/Reactome/
    Hallmark/TFT)."""
    results = {}
    for aspect, fname in files.items():
        gmt_path = go_dir / fname
        if not gmt_path.exists():
            print(f"  [{aspect}] file {gmt_path} non trovato, salto")
            continue
        res = go_enrichment(foreground_symbols, background_symbols, gmt_path)
        results[aspect] = res
        fout = outdir / f"{aspect}_{label}.csv"
        res.to_csv(fout, index=False)
        n_sig = (res["q_value_BH"] < sig_threshold).sum() if not res.empty else 0
        print(f"  [{aspect}] {len(res)} termini testati, {n_sig} con q<{sig_threshold} -> {fout}")
        if n_sig > 0:
            for _, r in res[res["q_value_BH"] < sig_threshold].head(15).iterrows():
                print(f"    - {r['termine_GO']}: {r['n_geni_foreground_nel_termine']}/{r['n_foreground']} geni, "
                      f"fold={r['fold_enrichment']:.1f}, q={r['q_value_BH']:.2e}")
    return results


def run_go_enrichment_bp_mf_cc(foreground_symbols, background_symbols, go_dir: Path, outdir: Path, label: str, sig_threshold: float = 0.05):
    """Esegue go_enrichment per BP, MF, CC (se i relativi .gmt sono presenti in go_dir) e salva 3 CSV."""
    files = {"GO_BP": "c5.bp.v7.0.symbols.gmt", "GO_MF": "c5.mf.v7.0.symbols.gmt", "GO_CC": "c5.cc.v7.0.symbols.gmt"}
    return run_gene_set_enrichment(foreground_symbols, background_symbols, go_dir, outdir, label, files, sig_threshold)


def run_pathway_enrichment(foreground_symbols, background_symbols, go_dir: Path, outdir: Path, label: str, sig_threshold: float = 0.05):
    """Esegue l'arricchimento su KEGG, Reactome, Hallmark e TFT (bersagli
    noti di fattori di trascrizione) per capire in quali pathway/processi
    ricadono i geni del foreground (es. i geni m5C-responsivi)."""
    files = {
        "KEGG": "c2.cp.kegg.v7.0.symbols.gmt",
        "REACTOME": "c2.cp.reactome.v7.0.symbols.gmt",
        "HALLMARK": "h.all.v7.0.symbols.gmt",
        "TFT_bersagli_TF": "c3.tft.v7.0.symbols.gmt",
    }
    return run_gene_set_enrichment(foreground_symbols, background_symbols, go_dir, outdir, label, files, sig_threshold)


def tf_target_lookup(tf_symbols, tft_gmt_path: Path, outdir: Path, label: str):
    """
    Per una lista di gene symbol di fattori di trascrizione, cerca nel
    database TFT (TRANSFAC, motivi noti) i termini che corrispondono a quel
    TF (es. 'RREB1_01') e riporta i geni bersaglio noti. Non tutti i TF
    hanno un motivo caratterizzato in questo database (in particolare la
    famiglia zinc-finger KRAB e' storicamente poco coperta): quelli assenti
    vengono segnalati esplicitamente, non silenziati.
    """
    if not tft_gmt_path.exists():
        print(f"  [target TF] file {tft_gmt_path} non trovato, salto")
        return None

    sets = load_gmt(tft_gmt_path)
    rows = []
    found_tfs, missing_tfs = [], []
    for tf in tf_symbols:
        matches = [t for t in sets if t.upper() == tf.upper() or t.upper().startswith(tf.upper() + "_")]
        if not matches:
            missing_tfs.append(tf)
            continue
        found_tfs.append(tf)
        for term in matches:
            targets = sorted(sets[term])
            rows.append({
                "TF": tf, "termine_TFT": term, "n_geni_bersaglio": len(targets),
                "geni_bersaglio": ",".join(targets),
            })

    result = pd.DataFrame(rows)
    fname = outdir / f"target_TF_{label}.csv"
    result.to_csv(fname, index=False)
    print(f"  [target TF] {len(found_tfs)}/{len(tf_symbols)} TF con bersagli noti in TFT (TRANSFAC) -> {fname}")
    if missing_tfs:
        print(f"  [target TF] NESSUN motivo noto in TFT per: {missing_tfs}")
    return result


def report_genes_of_interest(df: pd.DataFrame, gene_symbols, mod_name: str, outdir: Path):
    """
    Estrae tutti i siti di una data modifica per una lista di geni (per
    gene_symbol), con % di modifica per campione (condizione+replicato),
    utile per controllare "a mano" lo stato di modifica di geni specifici
    (es. gli enzimi scrittori/rimuoventi coinvolti nel pathway studiato).
    """
    if "gene_symbol" not in df.columns:
        print("  [geni di interesse] colonna gene_symbol assente (serve --gtf), salto")
        return None

    sub = df[(df["mod_name"] == mod_name) & df["gene_symbol"].isin(gene_symbols)].copy()
    if sub.empty:
        print(f"  [geni di interesse] nessun sito {mod_name} trovato per {gene_symbols}")
        return None

    cols = ["gene_symbol", "gene_id", "condition", "replicate", "chrom", "start", "end",
            "region", "fraction_modified_pct", "N_valid_cov"]
    cols = [c for c in cols if c in sub.columns]
    out = sub[cols].sort_values(["gene_symbol", "condition", "replicate", "start"])

    fname = outdir / (f"geni_interesse_{mod_name}_" + "_".join(gene_symbols) + ".csv")
    out.to_csv(fname, index=False)
    print(f"  [geni di interesse] {len(out)} righe ({sub['gene_symbol'].nunique()} geni trovati su {len(gene_symbols)} richiesti) -> {fname}")

    found = sorted(sub["gene_symbol"].unique())
    missing = sorted(set(gene_symbols) - set(found))
    if missing:
        print(f"  [geni di interesse] NON trovati (nessun sito {mod_name} chiamato): {missing}")

    # riepilogo per gene x condizione (media replicati)
    summary = (
        sub.groupby(["gene_symbol", "condition"])
        .agg(n_siti=("fraction_modified_pct", "size"), pct_mod_medio=("fraction_modified_pct", "mean"))
        .reset_index()
    )
    fname_sum = outdir / (f"geni_interesse_{mod_name}_riepilogo.csv")
    summary.to_csv(fname_sum, index=False)
    print(f"  [geni di interesse] riepilogo per gene/condizione -> {fname_sum}")

    # coverage grezza: TUTTE le posizioni chiamate (qualunque modifica) per
    # questi geni, per distinguere un vero calo biologico della modifica da
    # un semplice calo di coverage/espressione del trascritto
    cov_all = df[df["gene_symbol"].isin(gene_symbols)].copy()
    cov_summary = (
        cov_all.groupby(["gene_symbol", "condition", "replicate"])
        .size().reset_index(name="n_posizioni_tutte_le_mod")
    )
    fname_cov = outdir / (f"geni_interesse_{mod_name}_coverage.csv")
    cov_summary.to_csv(fname_cov, index=False)
    print(f"  [geni di interesse] coverage grezza (tutte le mod) -> {fname_cov}")

    _plot_genes_of_interest(sub, mod_name, outdir, cov_all=cov_all)
    return out


def _plot_genes_of_interest(sub: pd.DataFrame, mod_name: str, outdir: Path, cov_all: pd.DataFrame = None):
    sub = sub.copy()
    sub["sample_label"] = sub["condition"].astype(str) + "_" + sub["replicate"].astype(str)
    genes = sorted(sub["gene_symbol"].unique())
    # se abbiamo la coverage grezza, includiamo anche i geni presenti li' ma
    # senza siti per mod_name (altrimenti spariscono dal plot senza spiegazione)
    if cov_all is not None:
        genes = sorted(set(genes) | set(cov_all["gene_symbol"].unique()))
    n = len(genes)
    has_cov = cov_all is not None and not cov_all.empty
    nrows = 2 if has_cov else 1
    fig, axes = plt.subplots(nrows, n, figsize=(3.2 * n, 4.5 * nrows), sharey="row", squeeze=False)

    for j, gene in enumerate(genes):
        ax = axes[0][j]
        g = sub[sub["gene_symbol"] == gene]
        if not g.empty:
            ax.scatter(g["sample_label"], g["fraction_modified_pct"], alpha=0.6, s=25)
            means = g.groupby("sample_label")["fraction_modified_pct"].mean()
            ax.scatter(means.index, means.values, color="red", marker="_", s=400, linewidths=2, zorder=5)
        else:
            ax.text(0.5, 0.5, f"nessun sito {mod_name}", ha="center", va="center", transform=ax.transAxes, fontsize=8, color="grey")
        ax.set_title(gene, fontsize=10)
        ax.tick_params(axis="x", rotation=60)
        ax.set_ylim(0, 100)
    axes[0][0].set_ylabel(f"% modifica ({mod_name})")

    if has_cov:
        cov_all = cov_all.copy()
        cov_all["sample_label"] = cov_all["condition"].astype(str) + "_" + cov_all["replicate"].astype(str)
        for j, gene in enumerate(genes):
            ax = axes[1][j]
            c = cov_all[cov_all["gene_symbol"] == gene].groupby("sample_label").size()
            # ordina i campioni come nel pannello sopra, se possibile
            all_samples = sorted(set(sub["sample_label"]) | set(cov_all["sample_label"]))
            c = c.reindex(all_samples, fill_value=0)
            ax.bar(c.index, c.values, color="#888888")
            ax.tick_params(axis="x", rotation=60)
            ax.set_xlabel("")
        axes[1][0].set_ylabel("N. posizioni chiamate\n(TUTTE le modifiche)")

    fig.suptitle(f"Geni di interesse: {mod_name} (sopra, punti=siti singoli/media rossa) "
                 + ("e coverage grezza totale (sotto)" if has_cov else ""))
    plt.tight_layout()
    fname = outdir / f"geni_interesse_{mod_name}_plot.png"
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  [geni di interesse] salvato {fname}")


# --------------------------------------------------------------------------
# 6. TABELLE RIASSUNTIVE
# --------------------------------------------------------------------------

def summary_by_sample_mod(df: pd.DataFrame) -> pd.DataFrame:
    """Numero totale di siti modificati per campione e per tipo di modifica."""
    return (
        df.groupby(["condition", "replicate", "mod_name"])
        .size()
        .reset_index(name="n_siti")
        .sort_values(["condition", "replicate", "mod_name"])
    )


def summary_coverage(df: pd.DataFrame) -> pd.DataFrame:
    """
    Coverage (N_valid_cov, cioe' la profondita' di lettura passante-filtro
    alla posizione) per campione: totale (somma, indicativo del "volume" di
    dati/posizioni chiamate) E mediana per sito (indicativo della profondita'
    di sequenziamento media, indipendente dal numero di posizioni). Se un
    campione ha piu' siti totali ma la STESSA profondita' mediana per sito
    delle altre, l'eccesso di siti non e' spiegato da maggiore profondita'
    di sequenziamento.
    """
    tab = (
        df.groupby(["condition", "replicate"])
        .agg(
            n_posizioni_totali=("N_valid_cov", "size"),
            cov_totale=("N_valid_cov", "sum"),
            cov_mediana_per_sito=("N_valid_cov", "median"),
            cov_media_per_sito=("N_valid_cov", "mean"),
        )
        .reset_index()
        .sort_values(["condition", "replicate"])
    )
    return tab


def plot_coverage(tab: pd.DataFrame, outdir: Path):
    tab = tab.copy()
    tab["sample_label"] = tab["condition"].astype(str) + "_" + tab["replicate"].astype(str)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].bar(tab["sample_label"], tab["n_posizioni_totali"], color="#4C72B0")
    axes[0].set_ylabel("N. posizioni totali chiamate (tutte le mod)")
    axes[0].set_title("Volume di dati per campione")
    axes[0].tick_params(axis="x", rotation=60)

    axes[1].bar(tab["sample_label"], tab["cov_mediana_per_sito"], color="#55A868")
    axes[1].set_ylabel("Coverage mediana per sito (N_valid_cov)")
    axes[1].set_title("Profondita' di sequenziamento per campione\n(indipendente dal n. di posizioni)")
    axes[1].tick_params(axis="x", rotation=60)

    plt.tight_layout()
    fname = outdir / "coverage_per_campione.png"
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  [coverage] salvato {fname}")


def summary_by_condition_mod(df: pd.DataFrame) -> pd.DataFrame:
    """Media tra replicati del numero di siti per condizione e modifica."""
    per_sample = summary_by_sample_mod(df)
    return (
        per_sample.groupby(["condition", "mod_name"])["n_siti"]
        .agg(["mean", "std", "count"])
        .reset_index()
        .rename(columns={"mean": "n_siti_medio", "std": "n_siti_std", "count": "n_replicati"})
    )


def summary_by_level(df: pd.DataFrame) -> pd.DataFrame:
    """Distribuzione dei siti per livello di modifica (basso/medio/alto),
    per condizione E replicato separatamente."""
    return (
        df.groupby(["condition", "replicate", "mod_name", "mod_level"])
        .size()
        .reset_index(name="n_siti")
    )


def summary_by_gene_type(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby(["mod_name", "gene_type"])
        .size()
        .reset_index(name="n_siti")
        .sort_values("n_siti", ascending=False)
    )


# --------------------------------------------------------------------------
# 7. TABELLE RIASSUNTIVE (continua)
# --------------------------------------------------------------------------

def summary_by_region(df: pd.DataFrame) -> pd.DataFrame:
    """Distribuzione dei siti per regione genica, per condizione E
    replicato separatamente (oltre a un riepilogo aggregato)."""
    if "region" not in df.columns:
        return pd.DataFrame()
    return (
        df.groupby(["condition", "replicate", "mod_name", "region"])
        .size()
        .reset_index(name="n_siti")
        .sort_values(["condition", "replicate", "mod_name", "n_siti"], ascending=[True, True, True, False])
    )


def per_gene_status(df: pd.DataFrame) -> pd.DataFrame:
    """
    Tabella "stato di modifica per gene": per ciascun gene, elenca le
    modifiche osservate, quante volte, la percentuale media di modifica
    e (se disponibile) le regioni interessate.
    """
    agg = {
        "chrom": "first",
        "gene_type": "first",
        "fraction_modified_pct": ["mean", "max", "count"],
    }
    grouped = df.groupby(["gene_id", "gene_name", "mod_name"]).agg(agg)
    grouped.columns = ["chrom", "gene_type", "pct_mod_medio", "pct_mod_max", "n_siti"]
    grouped = grouped.reset_index()

    if "region" in df.columns:
        regions = (
            df.groupby(["gene_id", "mod_name"])["region"]
            .apply(lambda s: ",".join(sorted(set(s.dropna()))))
            .reset_index()
            .rename(columns={"region": "regioni"})
        )
        grouped = grouped.merge(regions, on=["gene_id", "mod_name"], how="left")

    if "transcript_used" in df.columns:
        tx = df.groupby(["gene_id"])["transcript_used"].first().reset_index()
        grouped = grouped.merge(tx, on="gene_id", how="left")

    if "gene_symbol" in df.columns:
        gs = df.groupby(["gene_id"])["gene_symbol"].first().reset_index()
        grouped = grouped.merge(gs, on="gene_id", how="left")

    return grouped.sort_values(["gene_name", "mod_name"])


# --------------------------------------------------------------------------
# 7. GRAFICI DIAGNOSTICI
# --------------------------------------------------------------------------

def plot_counts_by_condition(summary_cond: pd.DataFrame, outdir: Path):
    pivot = summary_cond.pivot(index="condition", columns="mod_name", values="n_siti_medio")
    ax = pivot.plot(kind="bar", figsize=(9, 5))
    ax.set_ylabel("N. siti modificati (media replicati)")
    ax.set_title("Siti modificati per condizione e tipo di modifica")
    plt.tight_layout()
    plt.savefig(outdir / "siti_per_condizione_e_modifica.png", dpi=150)
    plt.close()


def plot_level_distribution(summary_level: pd.DataFrame, outdir: Path):
    """Un plot per modifica; asse x = campione (condizione_replicato),
    cosi' ogni replicato e' visibile separatamente."""
    df = summary_level.copy()
    df["sample_label"] = df["condition"] + "_" + df["replicate"]
    for mod in df["mod_name"].unique():
        sub = df[df["mod_name"] == mod]
        pivot = sub.pivot(index="sample_label", columns="mod_level", values="n_siti").fillna(0)
        pivot = pivot[[c for c in MOD_LEVEL_LABELS if c in pivot.columns]]
        pivot = pivot.sort_index()
        ax = pivot.plot(kind="bar", stacked=True, figsize=(10, 5))
        ax.set_ylabel("N. siti")
        ax.set_xlabel("Campione (condizione_replicato)")
        ax.set_title(f"Distribuzione livello di modifica - {mod} (per replicato)")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig(outdir / f"livello_modifica_{mod}.png", dpi=150)
        plt.close()


def plot_region_distribution(summary_region: pd.DataFrame, outdir: Path):
    """Un plot per modifica; asse x = campione (condizione_replicato),
    stack per regione genica, cosi' ogni replicato e' visibile separatamente."""
    if summary_region.empty:
        return
    df = summary_region.copy()
    df["sample_label"] = df["condition"] + "_" + df["replicate"]
    region_order = ["5'UTR", "CDS", "esone_ncRNA", "introne", "3'UTR", "intergenica", "gene_non_in_gtf"]
    for mod in df["mod_name"].unique():
        sub = df[df["mod_name"] == mod]
        pivot = sub.pivot(index="sample_label", columns="region", values="n_siti").fillna(0)
        cols = [c for c in region_order if c in pivot.columns] + [c for c in pivot.columns if c not in region_order]
        pivot = pivot[cols].sort_index()
        ax = pivot.plot(kind="bar", stacked=True, figsize=(10, 5))
        ax.set_ylabel("N. siti")
        ax.set_xlabel("Campione (condizione_replicato)")
        ax.set_title(f"Distribuzione per regione genica - {mod} (per replicato)")
        plt.xticks(rotation=45, ha="right")
        plt.legend(title="Regione", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()
        plt.savefig(outdir / f"distribuzione_regioni_{mod}.png", dpi=150)
        plt.close()


# --------------------------------------------------------------------------
# 8. MAIN
# --------------------------------------------------------------------------

def main():
    # --------------------------------------------------------------------
    # PERCORSI DI DEFAULT (modifica qui se cambi macchina/struttura cartelle;
    # restano comunque sovrascrivibili da riga di comando con --input-dir,
    # --outdir, ecc. se serve puntare altrove in un run specifico)
    # --------------------------------------------------------------------
    DEFAULT_INPUT_DIR = Path("data/modkit")
    DEFAULT_OUTPUT_DIR = Path("results/modkit")
    # le tabelle leggere di annotazione si assumono nella stessa cartella dei tsv di input
    DEFAULT_TRANSCRIPT_SPANS = DEFAULT_INPUT_DIR / "mane_transcript_spans.tsv"
    DEFAULT_REGION_FEATURES = DEFAULT_INPUT_DIR / "mane_region_features.tsv"
    DEFAULT_ALL_GENE_SPANS = DEFAULT_INPUT_DIR / "all_gene_spans.tsv"

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR, help="cartella con i file *_annotated.tsv")
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTPUT_DIR, help="cartella di output per tabelle e grafici")
    ap.add_argument("--transcript-spans", type=Path, default=None,
                     help="tabella leggera pre-calcolata: un trascritto canonico per gene (chrom/start/end/strand). "
                          "Default: mane_transcript_spans.tsv dentro --input-dir")
    ap.add_argument("--region-features", type=Path, default=None,
                     help="tabella leggera pre-calcolata: intervalli 5'UTR/CDS/3'UTR/esone_ncRNA/introne per trascritto canonico. "
                          "Default: mane_region_features.tsv dentro --input-dir")
    ap.add_argument("--all-gene-spans", type=Path, default=None,
                     help="tabella leggera pre-calcolata: TUTTI i geni (per il recupero per coordinate). "
                          "Default: all_gene_spans.tsv dentro --input-dir")
    ap.add_argument("--skip-annotation", action="store_true", default=False,
                     help="salta del tutto l'annotazione di regione/gene_symbol (nessuna delle 3 tabelle sopra viene letta)")
    ap.add_argument("--min-coverage", type=int, default=10,
                     help="copertura minima (N_valid_cov) per tenere un sito, uguale per tutte le modifiche")
    ap.add_argument("--m6a-fail-fraction", type=float, default=0.25,
                     help="soglia FISSA di fail_fraction per m6A (per questa modifica corrisponde comunque a un percentile alto della sua distribuzione)")
    ap.add_argument("--adaptive-fail-percentile", type=float, default=90.0,
                     help="percentile di fail_fraction (calcolato SEPARATAMENTE per ciascuna modifica) usato come soglia per m5C/Pseudouridina/Inosina")
    ap.add_argument("--recover-by-coordinate", action="store_true", default=True,
                     help="tenta di recuperare i siti intergenica/gene_non_in_gtf con un lookup per sole coordinate")
    ap.add_argument("--no-recover-by-coordinate", dest="recover_by_coordinate", action="store_false")
    ap.add_argument("--metagene-mods", type=str, nargs="+", default=["m5C", "m6A", "Inosina", "Pseudouridina"],
                     help="modifiche per cui generare il metagene plot")
    ap.add_argument("--metagene-conditions", type=str, nargs="+", default=["WT", "STORM"],
                     help="condizioni da confrontare nel metagene plot (default: WT STORM; passa piu' valori o ometti per tutte)")
    ap.add_argument("--metagene-conditions-2", type=str, nargs="+", default=["Amr", "C412", "E11"],
                     help="secondo set di condizioni da confrontare nel metagene plot (default: Amr/C412=slow/E11=fast); passa una lista vuota per disattivare")
    ap.add_argument("--metagene-mods-2", type=str, nargs="+", default=["m6A"],
                     help="modifiche per cui generare il metagene plot secondario")
    ap.add_argument("--test-window", type=float, nargs=2, default=[1.9, 2.1],
                     help="finestra [x_min x_max] sul metagene (2=stop codon) per il test statistico")
    ap.add_argument("--test-cond-a", type=str, default="WT")
    ap.add_argument("--test-cond-b", type=str, default="STORM")
    ap.add_argument("--paired-mods", type=str, nargs="+", default=["m5C", "m6A"],
                     help="modifiche per cui fare il confronto sito-per-sito (paired)")
    ap.add_argument("--paired-cond-a", type=str, default="WT")
    ap.add_argument("--paired-cond-b", type=str, default="STORM")
    ap.add_argument("--enrich-mods", type=str, nargs="+", default=["m5C"],
                     help="modifiche per cui isolare/arricchire il sottoinsieme responsivo (delta forte)")
    ap.add_argument("--enrich-delta-thresholds", type=float, nargs="+", default=[-30.0, -20.0],
                     help="soglie di delta (cond_B - cond_A, una o piu') sotto le quali un sito e' 'responsivo'")
    ap.add_argument("--tf-top-n", type=int, default=20,
                     help="quanti TF piu' responsivi cercare nel database dei bersagli (TFT)")
    ap.add_argument("--genes-of-interest", type=str, nargs="+",
                     default=["NSUN2", "NSUN4", "NSUN5", "NSUN6", "NSUN7", "ALKBH1"],
                     help="gene symbol da controllare puntualmente (es. enzimi scrittori/rimuoventi)")
    ap.add_argument("--genes-of-interest-mods", type=str, nargs="+", default=["m6A"],
                     help="modifiche per cui generare il report sui geni di interesse")
    ap.add_argument("--go-enrichment", action="store_true", default=True,
                     help="esegui l'arricchimento GO (BP/MF/CC) sui geni responsivi del delta enrichment")
    ap.add_argument("--no-go-enrichment", dest="go_enrichment", action="store_false")
    ap.add_argument("--go-dir", type=Path, default=Path(__file__).resolve().parent / "go_data",
                     help="cartella con (o dove scaricare) i file .gmt GO BP/MF/CC")
    ap.add_argument("--go-sig-threshold", type=float, default=0.1,
                     help="soglia di q-value sotto cui un termine GO viene stampato/evidenziato come significativo")
    ap.add_argument("--gene-family-name", type=str, default="ZNF/ZFP (zinc-finger C2H2)",
                     help="nome della famiglia genica da testare per arricchimento (regex su gene_symbol)")
    ap.add_argument("--gene-family-pattern", type=str, default=r"^(ZNF|ZFP|ZKSCAN|ZBTB|ZSCAN)",
                     help="pattern regex (case-insensitive) sul gene_symbol per definire la famiglia genica")
    args = ap.parse_args()

    # le tabelle di annotazione, se non specificate esplicitamente, si
    # cercano di default nella stessa cartella dei file di input (segue
    # sempre --input-dir, anche se quest'ultimo viene sovrascritto)
    if args.transcript_spans is None:
        args.transcript_spans = args.input_dir / "mane_transcript_spans.tsv"
    if args.region_features is None:
        args.region_features = args.input_dir / "mane_region_features.tsv"
    if args.all_gene_spans is None:
        args.all_gene_spans = args.input_dir / "all_gene_spans.tsv"

    args.outdir.mkdir(parents=True, exist_ok=True)

    print("== 1. Caricamento dati ==")
    df = load_all(args.input_dir)

    print("== 1bis. Filtro di qualita' ==")
    df = filter_modkit_sites(df, min_coverage=args.min_coverage,
                              m6a_fail_fraction=args.m6a_fail_fraction,
                              adaptive_percentile=args.adaptive_fail_percentile)

    print("== 2. Classificazione per % di modifica ==")
    df = classify_mod_level(df)
    df = add_rna_class(df)

    have_annotation = (not args.skip_annotation) and args.transcript_spans.exists() and args.region_features.exists()
    if have_annotation:
        print("== 3. Annotazione regione genica (esone/introne/UTR) + gene symbol ==")
        df, annot = annotate_regions(df, args.transcript_spans, args.region_features)
        s_mapping = summary_mapping_status(df)
        s_mapping.to_csv(args.outdir / "09_stato_mappatura_trascritto.csv", index=False)
        plot_mapping_status(s_mapping, args.outdir)
        if args.recover_by_coordinate and args.all_gene_spans.exists():
            print("== 3bis. Tentativo di recupero per coordinate (intergenica / gene_non_in_gtf) ==")
            try:
                coord_index = LightCoordinateIndex(args.all_gene_spans)
                recover_unmapped_by_coordinate(df, coord_index, args.outdir)
            except ImportError:
                print("  ATTENZIONE: pacchetto 'intervaltree' non installato -> salto il recupero per "
                      "coordinate (il resto della pipeline continua normalmente). Per abilitarlo: "
                      "pip install intervaltree")
    else:
        print("== 3. GTF non fornito: salto l'annotazione di regione ==")
        annot = None

    print("== 4. Tabelle di riepilogo ==")
    s_sample = summary_by_sample_mod(df)
    s_cond = summary_by_condition_mod(df)
    s_level = summary_by_level(df)
    s_genetype = summary_by_gene_type(df)
    s_region = summary_by_region(df)
    s_gene = per_gene_status(df)

    s_sample.to_csv(args.outdir / "01_siti_per_campione_e_modifica.csv", index=False)
    s_cond.to_csv(args.outdir / "02_siti_per_condizione_e_modifica.csv", index=False)
    s_level.to_csv(args.outdir / "03_distribuzione_livello_modifica.csv", index=False)
    s_genetype.to_csv(args.outdir / "04_siti_per_gene_type.csv", index=False)
    if not s_region.empty:
        s_region.to_csv(args.outdir / "05_siti_per_regione_genica.csv", index=False)
    s_gene.to_csv(args.outdir / "06_stato_modifica_per_gene.csv", index=False)

    s_rna = summary_by_rna_class(df)
    s_rna.to_csv(args.outdir / "08_classificazione_rna.csv", index=False)
    plot_rna_class(s_rna, args.outdir)

    if annot is not None:
        gene_map = pd.DataFrame(
            sorted(annot.gene_symbol.items()), columns=["gene_id", "gene_symbol"]
        )
        gene_map.to_csv(args.outdir / "07_mapping_gene_id_symbol.csv", index=False)
        print(f"  Mapping gene_id -> gene_symbol salvato ({len(gene_map):,} geni, da {args.transcript_spans.name})")

    print("== 5. Grafici ==")
    plot_counts_by_condition(s_cond, args.outdir)
    plot_level_distribution(s_level, args.outdir)
    if not s_region.empty:
        plot_region_distribution(s_region, args.outdir)

    if annot is not None:
        print("== 6. Metagene plot ==")
        for mod in args.metagene_mods:
            # limitiamo il dataframe alla sola modifica corrente: riduce
            # sensibilmente la memoria necessaria e velocizza il calcolo
            mg_df = df[df["mod_name"] == mod].copy()
            if mg_df.empty:
                print(f"  [metagene] nessun dato per {mod}, salto")
                continue
            mg_df = compute_metagene_positions(mg_df, annot)
            for mode in ("density", "counts", "pct_mod"):
                plot_metagene(mg_df, mod, args.outdir, conditions=args.metagene_conditions, mode=mode)
            del mg_df

        if args.metagene_conditions_2:
            print(f"== 6bis. Metagene plot secondario ({','.join(args.metagene_conditions_2)}) ==")
            for mod in args.metagene_mods_2:
                mg_df = df[df["mod_name"] == mod].copy()
                if mg_df.empty:
                    print(f"  [metagene 2] nessun dato per {mod}, salto")
                    continue
                mg_df = compute_metagene_positions(mg_df, annot)
                label2 = "_".join(args.metagene_conditions_2)
                for mode in ("density", "counts", "pct_mod"):
                    plot_metagene(mg_df, mod, args.outdir, conditions=args.metagene_conditions_2, mode=mode,
                                  subset_label=label2)
                del mg_df

        print("== 7. Test statistico sulla finestra del metagene ==")
        for mod in args.metagene_mods:
            mg_df = df[df["mod_name"] == mod].copy()
            if mg_df.empty:
                continue
            mg_df = compute_metagene_positions(mg_df, annot)
            test_metagene_window(
                mg_df, mod, args.outdir,
                x_min=args.test_window[0], x_max=args.test_window[1],
                cond_a=args.test_cond_a, cond_b=args.test_cond_b,
            )
            del mg_df

        print("== 8. Confronto sito-per-sito (paired) ==")
        paired_results = {}
        for mod in args.paired_mods:
            paired_results[mod] = paired_site_comparison(df, mod, args.paired_cond_a, args.paired_cond_b, args.outdir)

        print("== 9. Arricchimento del sottoinsieme responsivo ==")
        go_ready = False
        pathway_ready = False
        if args.go_enrichment:
            go_ready = ensure_gmts(args.go_dir, GO_GMT_URLS)
            pathway_ready = ensure_gmts(args.go_dir, PATHWAY_GMT_URLS)

        for mod in args.enrich_mods:
            paired = paired_results.get(mod)
            if paired is None:
                paired = paired_site_comparison(df, mod, args.paired_cond_a, args.paired_cond_b, args.outdir)
            if paired is None or paired.empty:
                continue

            for delta_threshold in args.enrich_delta_thresholds:
                print(f"-- {mod}, soglia responsivita' delta<{delta_threshold} --")
                responsive, enrich_df, gene_stats = enrich_responsive_subset(
                    paired, mod, args.outdir, delta_threshold=delta_threshold
                )
                label = f"{mod}_delta_lt_{delta_threshold}"

                foreground_genes = []
                tf_gene_set = set()
                if args.go_enrichment and (go_ready or pathway_ready) and not gene_stats.empty:
                    gene_col = "gene_symbol" if "gene_symbol" in gene_stats.columns else "gene_id"
                    background_genes = gene_stats[gene_col].dropna().unique().tolist()
                    foreground_genes = gene_stats.loc[
                        gene_stats["n_siti_responsivi"] >= 1, gene_col
                    ].dropna().unique().tolist()
                    print(f"  [GO/pathway] {mod} (delta<{delta_threshold}): {len(foreground_genes):,} geni "
                          f"responsivi (foreground) su {len(background_genes):,} geni testati (background)")
                    if len(foreground_genes) >= 3:
                        if go_ready:
                            go_res = run_go_enrichment_bp_mf_cc(
                                foreground_genes, background_genes, args.go_dir, args.outdir,
                                label=label, sig_threshold=args.go_sig_threshold,
                            )
                            # geni annotati come "DNA binding transcription factor activity" (GO:MF),
                            # usati sotto per isolare i TF (non tutti i geni responsivi lo sono!)
                            mf = go_res.get("GO_MF")
                            if mf is not None and not mf.empty:
                                tf_row = mf[mf["termine_GO"] == "GO_DNA_BINDING_TRANSCRIPTION_FACTOR_ACTIVITY"]
                                if not tf_row.empty:
                                    tf_gene_set = set(tf_row.iloc[0]["geni"].split(","))
                        if pathway_ready:
                            run_pathway_enrichment(
                                foreground_genes, background_genes, args.go_dir, args.outdir,
                                label=label, sig_threshold=args.go_sig_threshold,
                            )
                    else:
                        print(f"  [GO/pathway] {mod}: troppo pochi geni responsivi ({len(foreground_genes)}) per un arricchimento sensato")

                print(f"  -- Arricchimento famiglia genica ({args.gene_family_name}) --")
                family_enrichment(
                    gene_stats, args.gene_family_name, args.gene_family_pattern, args.outdir,
                    mod_name=mod, delta_threshold=delta_threshold,
                )

                if annot is not None:
                    print(f"  -- Metagene e regione ristretti a {args.gene_family_name} --")
                    mg_df = df[df["mod_name"] == mod].copy()
                    mg_df = compute_metagene_positions(mg_df, annot)
                    mg_df_family = mg_df[mg_df["gene_symbol"].astype(str).str.match(args.gene_family_pattern, case=False)]
                    for m in ("density", "counts", "pct_mod"):
                        plot_metagene(mg_df_family, mod, args.outdir, conditions=[args.paired_cond_a, args.paired_cond_b],
                                      mode=m, subset_label=args.gene_family_name)
                    region_distribution_subset(
                        df, mod, args.gene_family_pattern, args.outdir,
                        subset_label=args.gene_family_name, conditions=[args.paired_cond_a, args.paired_cond_b],
                    )
                    del mg_df, mg_df_family

                # dettaglio sui TF piu' responsivi (SOLO quelli annotati come
                # fattori di trascrizione via GO:MF, non tutti i geni responsivi):
                # quali geni regolano (database TFT)
                if pathway_ready and tf_gene_set and not gene_stats.empty:
                    gene_col = "gene_symbol" if "gene_symbol" in gene_stats.columns else "gene_id"
                    top_tfs = (
                        gene_stats[gene_stats[gene_col].isin(tf_gene_set)]
                        .sort_values(["n_siti_responsivi", "frazione_responsiva"], ascending=False)
                        [gene_col].head(args.tf_top_n).tolist()
                    )
                    print(f"  -- Bersagli noti (TFT) dei {len(top_tfs)} TF piu' responsivi (tra i {len(tf_gene_set)} TF nel foreground) --")
                    tf_target_lookup(top_tfs, args.go_dir / "c3.tft.v7.0.symbols.gmt", args.outdir, label=label)


        print("== 10. Report su geni di interesse ==")
        for mod in args.genes_of_interest_mods:
            report_genes_of_interest(df, args.genes_of_interest, mod, args.outdir)

    # dataset completo (per analisi successive / caricamento rapido)
    df.to_parquet(args.outdir / "00_dataset_completo.parquet", index=False)

    print(f"\nFatto. Risultati in: {args.outdir}")


if __name__ == "__main__":
    main()
