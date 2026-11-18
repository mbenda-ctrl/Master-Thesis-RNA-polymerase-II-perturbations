# Transcriptional and epitranscriptional consequences of RNA Polymerase II perturbations

Code from my Master's thesis internship. The project integrates **Pol II ChIP-seq**, **RNA-seq**,
**nascent RNA** (Nanopore, metabolic labelling) and **m6A** calls (Nanopore direct RNA, modkit) in HEK293T cells to
study how different perturbations of RNA Polymerase II affect the transcription cycle, the transcriptome
and the m6A landscape.

> **In breve (IT)** – Codice del tirocinio magistrale: ChIP-seq di RNA Pol II (QC, allineamento, peak calling,
> DiffBind, annotazione, metagene), espressione differenziale RNA-seq e RNA nascente, analisi dei siti m6A da
> modkit e integrazione dei tre livelli (Venn, heatmap integrate, arricchimento GO).

## Experimental design

| Experiment | Control | Perturbations |
|---|---|---|
| Pol II elongation rate | **Amr** (α-amanitin-resistant Pol II) | **Fast** (fast-elongating mutant), **Slow** (slow-elongating mutant) |
| Treatments | **UT** (untreated) | **Storm** (METTL3 catalytic inhibitor), **NaCl** (hyperosmotic stress) |
| Hyperosmotic stress, public ChIP-seq ([GSE152062](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152062)) | **UT\*** | **KCl** (used as the Pol II reference for NaCl) |

Two biological replicates per condition. ChIP-seq: anti-Pol II IP + Input, Illumina paired-end, 61 bp reads.

## Repository structure

```
01_chipseq/                 Pol II ChIP-seq: QC -> alignment -> bigWig -> peaks -> DiffBind -> annotation
  metagene/                 metagene profiles (TSS / TES / scaled gene body), per RNA-seq category, differential
  public_dataset_KCl_UT/    processing of the public UT*/KCl dataset (GSE152062)
02_rnaseq/                  DESeq2 differential expression, DE lists, volcano plots, GO enrichment
03_nascent_rna/             per-gene aggregation of RNAkinet calls, DESeq2 on modified (nascent) reads, GO
04_modkit_m6a/              modkit site analysis, site filtering, differential m6A per gene
05_integration/             RNA-seq x nascent x Pol II x m6A: Venn/Fisher, integrated heatmaps, GO on clusters
config/                     DiffBind sample sheet and contrasts, SRA runs of the public dataset
environment/                conda environment (ChIP-seq tools), R and Python dependencies
```

Every script has its parameters at the top (or as command-line arguments) and a usage line in the header.
Default paths are relative to the repository root (`data/` for inputs, `results/` for outputs).

## 1. ChIP-seq (`01_chipseq/`)

| Step | Script | What it does |
|---|---|---|
| 1 | `01_fastqc_multiqc.sh` | FastQC on R1 and R2 of every sample; MultiQC reports (all, IP only, Input only) with R1/R2 merged per sample (`multiqc_config.yaml`). Optional read-only fastp reports (`RUN_FASTP=1`: `--disable_quality_filtering --disable_length_filtering --disable_adapter_trimming`) |
| 2 | `02_bowtie2_index_hg38.sh` | Bowtie2 index of GRCh38 (no-alt analysis set) |
| 3 | `03_align_bowtie2_pe.sh` | Bowtie2 2.4.5 paired-end alignment (default end-to-end, `--sensitive`), `samtools view` / `sort` / `index` (samtools 1.19.2). No duplicate removal |
| 4 | `04_bam_qc.sh` | `samtools flagstat`, multimapping estimates (MAPQ < 30, MAPQ 0, XS tag), MultiQC |
| 5 | `05_bigwig_cpm.sh` | deepTools `bamCoverage`: CPM, bin 10 bp, `--extendReads 200` |
| 6 | `06_bigwig_IPminusInput.sh` | `bigwigCompare --operation subtract` (IP − Input, CPM) |
| 7 | `07_macs3_peak_calling.sh` | MACS3 3.0.2, `-f BAMPE -g hs`; narrow peaks q < 0.01, broad peaks `--broad --broad-cutoff 0.1` |
| 8 | `08a_peaks_to_genes.py`, `08b_make_gene_bed.py` | One transcript per gene (MANE Select, GENCODE v48 basic, protein-coding); a peak is assigned to a gene if it falls in the transcription unit (TSS − 3 kb … TES + 300 bp). Consensus set: genes with a peak in ≥ 2 samples, unit ≥ 1 kb (13,542 genes; ≥ 5 kb for metagenes) |
| 9 | `09_annotate_peaks_ChIPseeker.R` | ChIPseeker `annotatePeak` on narrow and broad peaks, promoter = TSS ± 1 kb |
| 10 | `10_diffbind.R` | DiffBind 3.12.0: summit-centred 401 bp consensus (peaks in ≥ 2 samples), DESeq2 with library-size normalisation, report of all sites (`th = 1`) |
| 11 | `11_diffbind_volcano.py` | Volcano plots of the DiffBind reports |

**Metagenes (`01_chipseq/metagene/`)**, deepTools `computeMatrix` on IP − Input bigWigs, MANE genes ≥ 5 kb,
`--missingDataAsZero`:

- `01_run_metagene_profiles.sh` (+ `01a`): reference-point matrices on TSS and TES (± 3 kb, 10 bp bins) and
  scale-regions on the gene body (body scaled to 35 kb, ± 3 kb, 50 bp bins); per-sample mean profile in CPM
  and as log2(profile + 1).
- `02_make_bed_categorie.py`, `03_run_metagene_categorie.sh` (+ `03a`–`03c`): Pol II profiles of genes
  grouped by RNA-seq category (UP / NO CHANGE / DOWN; NO CHANGE = |log2FC| < 0.3 and baseMean ≥ 50).
  Differential profiles: log2((treated + 0.1) / (control + 0.1)) per gene and bin, 5-bin moving average,
  95 % bootstrap CI (500 resamplings of genes).
- `04_run_differential_metagene.sh` (+ `04a`): differential log2FC metagenes (scaled 35 kb, ± 3 kb, 100 bp bins)
  on the genes with differentially bound DiffBind sites.

## 2. RNA-seq differential expression (`02_rnaseq/`)

- `01_deseq2_rnaseq.R`: DESeq2 on the raw gene count matrix. One model per contrast built on its four samples
  (`~ condition`), no count pre-filter (DESeq2 independent filtering), median-of-ratios size factors, Wald test,
  Benjamini–Hochberg. LFC shrinkage with `lfcShrink(type = "normal")`; MLE log2FC also reported. Outputs: full
  tables, significant genes, volcano and MA plots, QC (PCA, sample distances, dispersion).
  Contrasts: Fast vs Amr, Slow vs Amr, Storm vs UT, NaCl vs UT.
- `02_crea_liste_rnaseq.py`: UP / DOWN lists used downstream: padj < 0.05 (NaCl vs UT: padj < 0.10) and
  |log2FC MLE| ≥ 0.58 (≥ 1.5-fold).
- `03_volcano_rnaseq.R`: volcano plots labelling the top 10 UP and top 10 DOWN genes.
- `04_go_rnaseq_updown.R`: GO enrichment of the UP / DOWN lists (universe = genes tested in the contrast).

## 3. Nascent RNA (`03_nascent_rna/`)

Nanopore reads classified as labelled (nascent) or not with **RNAkinet** (the read-level classification was
performed by a colleague and is not part of this repository).

- `01_aggrega_rnakinet_per_gene.py`: per gene and sample, total and modified (labelled) reads.
- `02_deseq2_nascent.R`: a single DESeq2 model on the modified-read counts of all twelve samples (`~ condition`),
  size factors computed on the modified reads; per contrast, genes with ≥ 5 modified reads in each of the four
  samples, BH correction on those genes only. Differentially synthesised genes: padj < 0.05 (NaCl 0.10) and
  |shrunken log2FC| ≥ 0.58.
- `03_go_nascent_updown.R`: GO enrichment of the nascent UP / DOWN lists.

## 4. m6A from modkit (`04_modkit_m6a/`)

Input: modkit per-site tables (bedMethyl, annotated with gene id/name/type). Base calling and `modkit` calling
were performed upstream and are not included.

- `01_modkit_analysis.py`: site counts per modification (m6A, m5C, pseudouridine, inosine), sample and condition,
  modification-level classes, distribution over gene regions (5′UTR, CDS, 3′UTR, introns, ncRNA exons) on
  MANE transcripts, per-gene modification status, metagene plots.
- `02_filtra_siti.py`: per-sample filter, coverage ≥ 98 reads and modified fraction ≥ 20 %.
- `03_diff_m6A.py`: differential m6A per protein-coding gene (treated vs control). Stoichiometry tested with a
  quasi-binomial GLM; differentially methylated genes = |log2FC stoichiometry| > 0.58 and BH FDR < 0.05.
  Change in the number of sites = difference ≥ 2 sites with consistent replicates.

## 5. Integration (`05_integration/`)

- `01_integrazione_rnaseq_nascent.R`: Venn diagrams and 2 × 2 contingency matrices between RNA-seq and nascent
  UP / DOWN genes.
- `02_chip_integrazione.R`: DiffBind regions assigned to MANE genes (gene body ± 1 kb); Venn diagrams of Pol II,
  RNA-seq and nascent genes with Fisher's exact test. ChIP thresholds: FDR < 0.05 (Slow vs Amr);
  FDR < 0.01 and |Fold| > 1.5 (KCl vs UT\*, paired with NaCl vs UT); all regions for Fast vs Amr and Storm vs UT.
- `heatmaps/`: integrated gene-level heatmaps with baseline columns (control level) and log2FC columns
  (Pol II, RNA-seq, nascent). Pol II per gene = weighted mean (weight 2^Conc) of the DiffBind Fold of the regions
  whose centre falls in the gene ± 1 kb. Colour scale saturated at ± 2, grey = not significant, white = missing;
  rows ordered by Ward clustering on the log2FC (missing = 0, clipped at ± 3), k = 4 clusters.
  `heatmap_guidata_rna.py` (genes DE in RNA-seq and/or nascent RNA), `heatmap_guidata_chip.py` (genes with
  significant Pol II changes), `heatmap_guidata_stech.py` (genes with changed m6A stoichiometry); optional m6A columns.
- `03_m6a_polII_venn.R`: overlap between differentially methylated genes and Pol II DiffBind genes.
- `go_enrichment/go_cluster_clusterProfiler.R`: GO Biological Process enrichment of the heatmap clusters with
  clusterProfiler `compareCluster(fun = "enrichGO")`, org.Hs.eg.db, term size 10–500, BH padj < 0.05,
  redundancy removed with `simplify` (Wang similarity ≥ 0.7).

## Software

FastQC, fastp 0.22.0, MultiQC, Bowtie2 2.4.5, samtools 1.19.2, deepTools 3.5.6, MACS3 3.0.2, R 4.3 with
DESeq2 1.42.0, DiffBind 3.12.0, ChIPseeker, clusterProfiler, org.Hs.eg.db; Python 3 with numpy, pandas, scipy,
statsmodels, matplotlib. See `environment/`.

```bash
conda env create -f environment/chipseq.yml      # ChIP-seq tools
pip install -r environment/requirements.txt      # Python packages for downstream scripts
Rscript environment/install_R_packages.R         # R / Bioconductor packages
```

## Data

Raw and processed data are not stored in this repository.

- Pol II ChIP-seq coverage tracks (bigWig, CPM): [Zenodo 17602869](https://zenodo.org/records/17602869)
- MACS3 peaks (bigBed): [Zenodo 17749415](https://zenodo.org/records/17749415)
- Public UT\*/KCl Pol II ChIP-seq: [GEO GSE152062](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152062)

## Author

Matteo Benda – Master's thesis internship, Università degli Studi di Milano-Bicocca.
