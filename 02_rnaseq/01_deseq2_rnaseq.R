#!/usr/bin/env Rscript
# =============================================================================
#  Analisi differenziale DESeq2 - HEK RNA-seq (conte grezze)
#
#  Confronti:
#     Fast  vs Amr      Slow vs Amr        (controllo = Amr)
#     Storm vs WT       NaCl vs WT         (controllo = WT)
#
#  Uso:  Rscript 01_deseq2_rnaseq.R [counts.tsv] [cartella_output]
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
# scrive un CSV usando il nome originale della colonna ID (gene_id / transcript_id)
write_res <- function(x, path) {
  names(x)[names(x) == "gene_id"] <- ID_COL
  write.csv(x, path, row.names = FALSE)
}

# ----------------------------------------------------------------------------
# 1. CONFIGURAZIONE  (modifica solo qui)
# ----------------------------------------------------------------------------
args        <- commandArgs(trailingOnly = TRUE)
COUNTS_FILE <- if (length(args) >= 1) args[1] else "data/expressed_genes_counts_matrix.tsv"
OUT_DIR     <- if (length(args) >= 2) args[2] else "results/rnaseq/RNAseq"

PADJ_CUT <- 0.05   # soglia FDR
LFC_CUT  <- 0.58   # |log2FC| per la lista dei "significativi" (0.58 ≈ fold change 1.5x)
USE_PREFILTER <- FALSE  # TRUE = applica il filtro sotto; FALSE = tieni tutti i geni
MIN_COUNT <- 10    # filtro: almeno MIN_COUNT conte...
MIN_SAMPLES <- 2   # ...in almeno MIN_SAMPLES campioni (= dimensione del gruppo più piccolo)
# Independent filtering di DESeq2: esclude dal calcolo dell'FDR (padj = NA)
# i geni con conte medie basse. TRUE = default DESeq2 (consigliato).
INDEP_FILTER <- TRUE

# Filtro per confronto: tiene solo i geni con almeno CONTRAST_MIN_COUNT conte
# in TUTTI i campioni del confronto (es. Fast_r1, Fast_r2, Amr_r1, Amr_r2 >= 10).
# Il padj (BH) viene ricalcolato solo su quei geni; sostituisce l'independent filtering.
PER_CONTRAST_FILTER <- FALSE
CONTRAST_MIN_COUNT  <- 10

# Correzione per test multipli usata per definire i significativi:
#   "BH"     = Benjamini-Hochberg (default DESeq2)
#   "qvalue" = Storey q-value (più permissivo: stima la quota di geni non DE)
#   "IHW"    = Independent Hypothesis Weighting, pesato sul baseMean
# La colonna padj_BH contiene sempre il valore BH standard di DESeq2.
PADJ_METHOD <- "BH"

# Nome nel file -> nome della condizione (E11 = Fast, C412 = Slow; WT = UT)
# I campioni vengono rinominati <Condizione>_r1 / _r2 (es. Fast_r1, UT_r2).
COND_MAP <- c(Amr   = "Amr",
              E11   = "Fast",
              C412  = "Slow",
              NACL  = "NaCl",
              STORM = "Storm",
              WT    = "UT")

# Confronti: c(trattato, controllo)
CONTRASTS <- list(
  Fast_vs_Amr = c("Fast",  "Amr"),
  Slow_vs_Amr = c("Slow",  "Amr"),
  Storm_vs_UT = c("Storm", "UT"),
  NaCl_vs_UT  = c("NaCl",  "UT")
)

# FALSE = un unico modello con tutti i 12 campioni (consigliato con 2 repliche:
#         la dispersione è stimata su più campioni).
# TRUE  = un modello separato per ogni confronto, solo con i 4 campioni coinvolti.
FIT_SEPARATE <- TRUE

# ----------------------------------------------------------------------------
# 2. LETTURA DELLA MATRICE
#    Funziona con geni (gene_id, ENSG) e trascritti (transcript_id, ENST).
#    La riga di intestazione può essere in cima o in fondo: la si cerca.
#    Le righe di riepilogo di htseq-count (__ambiguous, __no_feature...) si tolgono.
# ----------------------------------------------------------------------------
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
lines <- readLines(COUNTS_FILE)
hdr_i <- grep("^(gene_id|transcript_id)\t", lines)
if (length(hdr_i) != 1) stop("Riga di intestazione 'gene_id'/'transcript_id' non trovata (o duplicata).")
header <- strsplit(lines[hdr_i], "\t")[[1]]
ID_COL <- header[1]                       # "gene_id" oppure "transcript_id"
header[1] <- "gene_id"                    # nome interno
raw <- read.delim(text = lines[-hdr_i], header = FALSE, col.names = header,
                  check.names = FALSE, stringsAsFactors = FALSE)
special <- grepl("^__", raw$gene_id)
if (any(special)) {
  message("Righe di riepilogo htseq-count rimosse: ", paste(raw$gene_id[special], collapse = ", "))
  raw <- raw[!special, ]
}
rownames(raw) <- raw$gene_id
counts <- as.matrix(raw[, -1])
storage.mode(counts) <- "numeric"
if (any(abs(counts - round(counts)) > 1e-6))
  message("Attenzione: conte non intere (es. RSEM/Salmon) -> arrotondate.")
counts <- round(counts)
storage.mode(counts) <- "integer"
message(sprintf("Matrice: %d %s x %d campioni", nrow(counts),
                ifelse(ID_COL == "transcript_id", "trascritti", "geni"), ncol(counts)))

# ----------------------------------------------------------------------------
# 3. METADATI
# ----------------------------------------------------------------------------
orig_names <- colnames(counts)

parse_meta <- function(nm) {
  key  <- sub("^HEK_(.+)_r[0-9]+$", "\\1", nm)
  cond <- unname(COND_MAP[key])
  rep  <- sub(".*_(r[0-9]+)$", "\\1", nm)
  data.frame(sample    = paste0(cond, "_", rep),
             condition = cond,
             replicate = rep, stringsAsFactors = FALSE)
}
coldata <- parse_meta(colnames(counts))
coldata$original_label <- orig_names
if (anyNA(coldata$condition)) stop("Campioni senza condizione: controlla COND_MAP.")
if (anyDuplicated(coldata$sample)) stop("Nomi campione duplicati: controlla COND_MAP.")
colnames(counts)  <- coldata$sample
rownames(coldata) <- coldata$sample

# ordina per condizione
ord <- order(match(coldata$condition, unique(COND_MAP)), coldata$replicate)
counts  <- counts[, ord]
coldata <- coldata[ord, ]
coldata$condition <- factor(coldata$condition, levels = unique(COND_MAP))
write.csv(coldata, file.path(OUT_DIR, "sample_sheet_used.csv"), row.names = FALSE)
print(coldata)

# ----------------------------------------------------------------------------
# 4. ANNOTAZIONE (ENSEMBL -> SYMBOL), se org.Hs.eg.db è installato
# ----------------------------------------------------------------------------
ens <- sub("\\..*$", "", rownames(counts))
annot <- data.frame(gene_id = rownames(counts), ensembl = ens,
                    symbol = NA_character_, description = NA_character_)
ANN_KEY <- if (all(grepl("^ENST", ens))) "ENSEMBLTRANS" else "ENSEMBL"
if (has_pkg("org.Hs.eg.db") && has_pkg("AnnotationDbi")) {
  suppressMessages({
    annot$symbol <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, ens, "SYMBOL",
                                          ANN_KEY, multiVals = "first")
    annot$description <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, ens, "GENENAME",
                                               ANN_KEY, multiVals = "first")
  })
} else message("org.Hs.eg.db non installato: niente simboli genici.")

# Trascritti: org.Hs.eg.db copre pochi ENST. Se è installato EnsDb.Hsapiens.v86
# (BiocManager::install("EnsDb.Hsapiens.v86")) si usa quello, che dà anche
# il gene di appartenenza (ENSG) e il tipo di trascritto.
if (ANN_KEY == "ENSEMBLTRANS") {
  annot$parent_gene <- NA_character_
  annot$tx_biotype  <- NA_character_
  if (has_pkg("EnsDb.Hsapiens.v86") && has_pkg("ensembldb")) {
    edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    tx  <- suppressWarnings(suppressMessages(
      ensembldb::select(edb, keys = unique(ens), keytype = "TXID",
                        columns = c("TXID", "GENEID", "GENENAME", "TXBIOTYPE"))))
    m <- match(ens, tx$TXID)
    annot$parent_gene <- tx$GENEID[m]
    annot$tx_biotype  <- tx$TXBIOTYPE[m]
    annot$symbol <- ifelse(is.na(tx$GENENAME[m]) | tx$GENENAME[m] == "",
                           annot$symbol, tx$GENENAME[m])
    message(sprintf("EnsDb v86: annotati %d / %d trascritti", sum(!is.na(m)), length(ens)))
  } else message("Per annotare bene i trascritti installa EnsDb.Hsapiens.v86.")
}
rownames(annot) <- annot$gene_id
label_of <- function(ids) ifelse(is.na(annot[ids, "symbol"]), annot[ids, "ensembl"],
                                 annot[ids, "symbol"])

# ----------------------------------------------------------------------------
# 5. FUNZIONI
# ----------------------------------------------------------------------------
make_dds <- function(cts, cd) {
  dds <- DESeqDataSetFromMatrix(cts, cd, design = ~ condition)
  keep <- if (USE_PREFILTER) rowSums(counts(dds) >= MIN_COUNT) >= MIN_SAMPLES
          else rowSums(counts(dds)) > 0   # toglie solo i geni con zero conte ovunque
  dds <- dds[keep, ]
  dds$condition <- droplevels(dds$condition)
  DESeq(dds, quiet = TRUE)
}

# shrinkage del LFC: apeglm > ashr > normal (in base a cosa è installato)
shrink_lfc <- function(dds, trt, ref, res) {
  if (has_pkg("apeglm")) {
    d2 <- dds
    d2$condition <- relevel(d2$condition, ref = ref)
    d2 <- nbinomWaldTest(d2, quiet = TRUE)
    coef <- paste0("condition_", trt, "_vs_", ref)
    s <- lfcShrink(d2, coef = coef, type = "apeglm", quiet = TRUE)
    attr(s, "method") <- "apeglm"
  } else if (has_pkg("ashr")) {
    s <- lfcShrink(dds, contrast = c("condition", trt, ref), res = res,
                   type = "ashr", quiet = TRUE)
    attr(s, "method") <- "ashr"
  } else {
    s <- lfcShrink(dds, contrast = c("condition", trt, ref), type = "normal", quiet = TRUE)
    attr(s, "method") <- "normal"
  }
  s
}

# p aggiustati secondo PADJ_METHOD; si parte dai geni che passano
# l'independent filtering di DESeq2 (padj BH non NA)
adjust_p <- function(res) {
  ok  <- !is.na(res$padj)
  out <- rep(NA_real_, nrow(res))
  if (PADJ_METHOD == "BH") return(res$padj)
  if (PADJ_METHOD == "qvalue") {
    if (!has_pkg("qvalue")) stop("Installa qvalue: BiocManager::install('qvalue')")
    out[ok] <- qvalue::qvalue(res$pvalue[ok])$qvalues
  } else if (PADJ_METHOD == "IHW") {
    if (!has_pkg("IHW")) stop("Installa IHW: BiocManager::install('IHW')")
    ok <- !is.na(res$pvalue)
    out[ok] <- IHW::adj_pvalues(IHW::ihw(res$pvalue[ok], res$baseMean[ok], alpha = PADJ_CUT))
  } else stop("PADJ_METHOD non valido")
  out
}

volcano_plot <- function(df, title) {
  df$status <- "NS"
  df$status[df$padj < PADJ_CUT & df$log2FoldChange >=  LFC_CUT] <- "Up"
  df$status[df$padj < PADJ_CUT & df$log2FoldChange <= -LFC_CUT] <- "Down"
  df <- df[!is.na(df$padj), ]
  df$mlog10p <- -log10(pmax(df$padj, 1e-300))
  top <- head(df[df$status != "NS", ][order(df$padj[df$status != "NS"]), ], 20)
  n_up <- sum(df$status == "Up"); n_dn <- sum(df$status == "Down")
  p <- ggplot(df, aes(log2FoldChange, mlog10p, colour = status)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_colour_manual(values = c(Down = "#2166AC", NS = "grey75", Up = "#B2182B")) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = 2, linewidth = 0.3) +
    geom_hline(yintercept = -log10(PADJ_CUT), linetype = 2, linewidth = 0.3) +
    labs(title = title, subtitle = sprintf("Up: %d   Down: %d   (padj %s < %.2g, |LFC| >= %g)",
                                           n_up, n_dn, PADJ_METHOD, PADJ_CUT, LFC_CUT),
         x = "log2 fold change (shrunk)", y = "-log10 padj", colour = NULL) +
    theme_bw()
  if (has_pkg("ggrepel") && nrow(top) > 0)
    p <- p + ggrepel::geom_text_repel(data = top, aes(label = label), size = 2.8,
                                      colour = "black", max.overlaps = 30)
  p
}

# ----------------------------------------------------------------------------
# 6. MODELLO GLOBALE + QC
# ----------------------------------------------------------------------------
dds_all <- make_dds(counts, coldata)
message(sprintf("Geni dopo il filtro: %d", nrow(dds_all)))
vsd <- vst(dds_all, blind = TRUE)

write_res(cbind(annot[rownames(dds_all), c("gene_id", "symbol")],
                counts(dds_all, normalized = TRUE)),
          file.path(OUT_DIR, "normalized_counts.csv"))
write_res(cbind(annot[rownames(vsd), c("gene_id", "symbol")], assay(vsd)),
          file.path(OUT_DIR, "vst_counts.csv"))

# PCA
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE, ntop = 500)
pv  <- round(100 * attr(pca, "percentVar"))
pca$sample <- colnames(vsd)
p_pca <- ggplot(pca, aes(PC1, PC2, colour = condition)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1: ", pv[1], "%"), y = paste0("PC2: ", pv[2], "%"),
       title = "PCA (VST, top 500 geni variabili)") +
  theme_bw()
if (has_pkg("ggrepel"))
  p_pca <- p_pca + ggrepel::geom_text_repel(aes(label = sample), size = 2.8, show.legend = FALSE)
ggsave(file.path(OUT_DIR, "QC_PCA.pdf"), p_pca, width = 7, height = 5.5)

# distanze tra campioni
sd_mat <- as.matrix(dist(t(assay(vsd))))
ann_col <- data.frame(condition = coldata[colnames(vsd), "condition"], row.names = colnames(vsd))
pheatmap(sd_mat, annotation_col = ann_col, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", main = "Distanze tra campioni (VST)",
         filename = file.path(OUT_DIR, "QC_sample_distances.pdf"), width = 8, height = 7)
pheatmap(cor(assay(vsd)), annotation_col = ann_col, display_numbers = TRUE,
         number_format = "%.3f", fontsize_number = 6, main = "Correlazione di Pearson (VST)",
         filename = file.path(OUT_DIR, "QC_sample_correlation.pdf"), width = 9, height = 8)

pdf(file.path(OUT_DIR, "QC_dispersion.pdf")); plotDispEsts(dds_all); invisible(dev.off())

# ----------------------------------------------------------------------------
# 7. CONFRONTI
# ----------------------------------------------------------------------------
summary_tab <- list()
for (nm in names(CONTRASTS)) {
  trt <- CONTRASTS[[nm]][1]; ref <- CONTRASTS[[nm]][2]
  message("== ", nm)

  if (FIT_SEPARATE) {
    sel <- coldata$condition %in% c(trt, ref)
    dds <- make_dds(counts[, sel], droplevels(coldata[sel, ]))
  } else dds <- dds_all

  smp_t <- colnames(dds)[dds$condition == trt]
  smp_r <- colnames(dds)[dds$condition == ref]

  res <- results(dds, contrast = c("condition", trt, ref), alpha = PADJ_CUT,
                 independentFiltering = INDEP_FILTER && !PER_CONTRAST_FILTER)
  if (PER_CONTRAST_FILTER) {
    keep_c <- rowSums(counts(dds)[, c(smp_r, smp_t)] >= CONTRAST_MIN_COUNT) ==
              length(c(smp_r, smp_t))
    res$padj <- NA
    res$padj[keep_c] <- p.adjust(res$pvalue[keep_c], method = "BH")
    message(sprintf("   geni con >= %d conte in tutti i campioni del confronto: %d",
                    CONTRAST_MIN_COUNT, sum(keep_c)))
  }
  shr <- shrink_lfc(dds, trt, ref, res)
  nc <- counts(dds, normalized = TRUE)

  df <- data.frame(
    gene_id        = rownames(res),
    symbol         = annot[rownames(res), "symbol"],
    description    = annot[rownames(res), "description"],
    baseMean       = res$baseMean,
    log2FoldChange = shr$log2FoldChange,      # LFC con shrinkage (per ranking/plot)
    lfcSE          = shr$lfcSE,
    log2FC_MLE     = res$log2FoldChange,      # LFC non "shrunk"
    stat           = res$stat,
    pvalue         = res$pvalue,
    padj_BH        = res$padj,
    padj           = adjust_p(res),           # metodo scelto in PADJ_METHOD
    mean_ctrl      = rowMeans(nc[, smp_r, drop = FALSE]),
    mean_trt       = rowMeans(nc[, smp_t, drop = FALSE]),
    nc[, c(smp_r, smp_t)],
    check.names = FALSE)
  names(df)[names(df) == "mean_ctrl"] <- paste0("mean_norm_", ref)
  names(df)[names(df) == "mean_trt"]  <- paste0("mean_norm_", trt)
  if ("parent_gene" %in% names(annot)) {           # solo per i trascritti
    extra <- data.frame(parent_gene = annot[df$gene_id, "parent_gene"],
                        tx_biotype  = annot[df$gene_id, "tx_biotype"])
    df <- data.frame(df[, 1:3], extra, df[, -(1:3)], check.names = FALSE)
  }
  df$label <- label_of(df$gene_id)
  df <- df[order(df$padj, -abs(df$log2FoldChange), na.last = TRUE), ]

  sig <- subset(df, !is.na(padj) & padj < PADJ_CUT & abs(log2FoldChange) >= LFC_CUT)

  out_df <- df[, setdiff(names(df), "label")]
  write_res(out_df, file.path(OUT_DIR, paste0(nm, "_all_genes.csv")))
  write_res(sig[, setdiff(names(sig), "label")],
            file.path(OUT_DIR, paste0(nm, "_significant.csv")))

  # MA plot
  pdf(file.path(OUT_DIR, paste0(nm, "_MAplot.pdf")), width = 7, height = 5)
  plotMA(shr, ylim = c(-6, 6), alpha = PADJ_CUT,
         main = sprintf("%s  (shrinkage: %s)", nm, attr(shr, "method")))
  invisible(dev.off())

  ggsave(file.path(OUT_DIR, paste0(nm, "_volcano.pdf")),
         volcano_plot(df, paste(trt, "vs", ref)), width = 7, height = 6)

  # heatmap top 50 DEG (z-score VST sui campioni del confronto)
  top_ids <- head(sig$gene_id, 50)
  if (length(top_ids) >= 2) {
    m <- assay(vsd)[top_ids, c(smp_r, smp_t)]
    m <- t(scale(t(m)))
    rownames(m) <- label_of(top_ids)
    pheatmap(m, annotation_col = ann_col[c(smp_r, smp_t), , drop = FALSE],
             cluster_cols = FALSE, fontsize_row = 6,
             main = sprintf("%s - top %d DEG (z-score VST)", nm, length(top_ids)),
             filename = file.path(OUT_DIR, paste0(nm, "_heatmap_top50.pdf")),
             width = 6, height = 9)
  }

  summary_tab[[nm]] <- data.frame(
    contrast   = nm, treated = trt, control = ref,
    genes_tested = sum(!is.na(res$padj)),
    padj_method  = PADJ_METHOD,
    padj_sig     = sum(df$padj < PADJ_CUT, na.rm = TRUE),
    up           = sum(sig$log2FoldChange > 0),
    down         = sum(sig$log2FoldChange < 0),
    shrinkage    = attr(shr, "method"))
}

summary_df <- do.call(rbind, summary_tab)
write.csv(summary_df, file.path(OUT_DIR, "DEG_summary.csv"), row.names = FALSE)
print(summary_df, row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("Fatto. Risultati in: ", normalizePath(OUT_DIR))
