#!/usr/bin/env Rscript
# =============================================================================
#  DESeq2 sull'RNA nascente (reads modificate 5EU, RNAkinet) - HEK
#
#  Impostazione:
#   - UN modello DESeq2 con tutti i 12 campioni (~ condition)
#   - size factor stimati UNA volta, su tutti i 12 campioni insieme
#   - per ogni confronto si tengono i geni con >= MIN_READS reads modificate
#     in CIASCUNO dei 4 campioni del confronto; il padj (BH) è calcolato
#     solo su quei geni
#
#  Uso:
#   Rscript 02_deseq2_nascent.R <cartella_csv> <cartella_output>
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

# -----------------------------------------------------------------------------
# PASSO 0. Parametri
# -----------------------------------------------------------------------------
args      <- commandArgs(trailingOnly = TRUE)
INPUT_DIR <- if (length(args) >= 1) args[1] else "input"
OUT_DIR   <- if (length(args) >= 2) args[2] else "DESeq2_nascent_global"

MIN_READS <- 5        # reads modificate minime in OGNI campione del confronto
PADJ_CUT  <- 0.05     # soglia FDR di default
# soglie FDR specifiche per confronto (sovrascrivono PADJ_CUT)
PADJ_BY_CONTRAST <- c(NaCl_vs_UT = 0.10)
LFC_CUT   <- 0.58     # |log2FC| per UP/DOWN (0.58 ~ fold change 1.5x)
LFC_TYPE  <- "shrunk" # "shrunk" (lfcShrink) oppure "MLE" (log2FC grezzo di DESeq2)

# Su quali conte si stimano i size factor (sempre su tutti i 12 campioni):
# "median_ratio"       = median-of-ratios sulle reads MODIFICATE (default DESeq2)
# "mod_libsize"        = totale delle reads modificate di ogni campione
# "total_median_ratio" = median-of-ratios sulle reads TOTALI (modificate + non)
SF_MODE <- "median_ratio"


COND_MAP  <- c(AMR = "Amr", E11 = "Fast", C412 = "Slow",
               NACL = "NaCl", STORM = "Storm", WT = "UT")
CONTRASTS <- list(Fast_vs_Amr = c("Fast", "Amr"),
                  Slow_vs_Amr = c("Slow", "Amr"),
                  Storm_vs_UT = c("Storm", "UT"),
                  NaCl_vs_UT  = c("NaCl", "UT"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# PASSO 1. Lettura dei CSV RNAkinet -> matrici gene x campione
#   Un file per campione. Il nome del campione (es. AMR_R1) si prende dal nome
#   del file. Il gene_id perde la versione (.16) e le righe dello stesso gene
#   vengono sommate, come nella pipeline Python.
# -----------------------------------------------------------------------------
files <- list.files(INPUT_DIR, pattern = "rnakinet.*\\.csv$", full.names = TRUE)
if (length(files) != 12) warning("Trovati ", length(files), " file invece di 12")

read_one <- function(f) {
  s <- sub(".*HEK_([A-Za-z0-9]+_R[0-9]).*", "\\1", basename(f))
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$gid <- sub("\\.[0-9]+$", "", d$gene_id)
  data.frame(sample = s,
             gid  = d$gid,
             name = d$gene_name,
             mod  = d$n_reads_modificate_5eu,
             tot  = d$n_reads)
}
long <- do.call(rbind, lapply(files, read_one))

# xtabs somma i valori per (gene, campione); i geni assenti in un campione = 0
mod <- unclass(xtabs(mod ~ gid + sample, data = long))
tot <- unclass(xtabs(tot ~ gid + sample, data = long))
mod <- round(mod); tot <- round(tot)
storage.mode(mod) <- "integer"; storage.mode(tot) <- "integer"
gene_name <- tapply(long$name, long$gid, `[`, 1)
message(sprintf("Passo 1: %d geni x %d campioni", nrow(mod), ncol(mod)))

# -----------------------------------------------------------------------------
# PASSO 2. Nomi finali dei campioni (<Condizione>_r1 / _r2)
# -----------------------------------------------------------------------------
orig <- colnames(mod)
lab  <- orig
cond <- unname(COND_MAP[sub("_R[0-9]$", "", lab)])
repl <- tolower(sub(".*_(R[0-9])$", "\\1", lab))
coldata <- data.frame(sample = paste0(cond, "_", repl), condition = cond,
                      replicate = repl, file_label = orig, stringsAsFactors = FALSE)
ord <- order(match(coldata$condition, unique(COND_MAP)), coldata$replicate)
coldata <- coldata[ord, ]
mod <- mod[, coldata$file_label]; tot <- tot[, coldata$file_label]
colnames(mod) <- colnames(tot) <- rownames(coldata) <- coldata$sample
coldata$condition <- factor(coldata$condition, levels = unique(COND_MAP))
message("Passo 2: ", paste(coldata$sample, collapse = " "))

# -----------------------------------------------------------------------------
# PASSO 3. Size factor GLOBALI (una sola volta, su tutti i 12 campioni)
#   median_ratio: per ogni gene presente (> 0) in tutti i 12 campioni si calcola
#   la media geometrica tra campioni; per ogni campione il size factor è la
#   MEDIANA dei rapporti conta / media geometrica. Si usa la matrice completa,
#   prima di qualsiasi filtro.
# -----------------------------------------------------------------------------
if (SF_MODE == "median_ratio") {
  sf <- estimateSizeFactorsForMatrix(mod)
} else if (SF_MODE == "total_median_ratio") {
  sf <- estimateSizeFactorsForMatrix(tot)
} else {
  lib <- colSums(mod)
  sf  <- lib / exp(mean(log(lib)))
}
coldata$reads_modified <- colSums(mod)
coldata$reads_total    <- colSums(tot)
coldata$pct_modified   <- round(100 * coldata$reads_modified / coldata$reads_total, 2)
coldata$size_factor    <- round(sf, 4)
message(sprintf("Passo 3: size factor (%s) calcolati su %d geni presenti in tutti i campioni",
                SF_MODE, sum(rowSums(mod > 0) == ncol(mod))))
print(coldata[, c("file_label", "condition", "reads_modified", "pct_modified", "size_factor")])

# -----------------------------------------------------------------------------
# PASSO 4. Filtro dei geni per confronto
#   Per ogni confronto: TRUE se il gene ha >= MIN_READS reads modificate in
#   ciascuno dei 4 campioni (2 trattati + 2 controlli).
#   Nel modello entrano i geni che passano il filtro in ALMENO un confronto.
# -----------------------------------------------------------------------------
keep_by_contrast <- sapply(CONTRASTS, function(cc) {
  s <- coldata$sample[coldata$condition %in% cc]
  rowSums(mod[, s] >= MIN_READS) == length(s)
})
keep_model <- rowSums(keep_by_contrast) > 0
message(sprintf("Passo 4: geni nel modello %d; per confronto: %s", sum(keep_model),
                paste(names(CONTRASTS), colSums(keep_by_contrast), sep = "=", collapse = ", ")))

# -----------------------------------------------------------------------------
# PASSO 5. Modello DESeq2 unico con i 12 campioni
#   I size factor del passo 3 vengono imposti a mano: DESeq() non li ricalcola
#   ("using pre-existing size factors"). Poi stima le dispersioni usando tutti
#   i 12 campioni e adatta il modello negative binomiale ~ condition.
# -----------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(mod[keep_model, ], coldata, design = ~ condition)
sizeFactors(dds) <- sf[colnames(dds)]
dds <- DESeq(dds, quiet = TRUE)
message("Passo 5: modello stimato; trend dispersione = ",
        attr(dispersionFunction(dds), "fitType"))

vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE, ntop = 500)
pv  <- round(100 * attr(pca, "percentVar")); pca$sample <- colnames(vsd)
p <- ggplot(pca, aes(PC1, PC2, colour = condition)) + geom_point(size = 3) +
  labs(x = paste0("PC1: ", pv[1], "%"), y = paste0("PC2: ", pv[2], "%"),
       title = "PCA reads modificate") + theme_bw()
if (has_pkg("ggrepel")) p <- p + ggrepel::geom_text_repel(aes(label = sample), size = 2.8, show.legend = FALSE)
ggsave(file.path(OUT_DIR, "QC_PCA.pdf"), p, width = 7, height = 5.5)
write.csv(coldata, file.path(OUT_DIR, "sample_sheet_size_factors.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# PASSO 6. Test per ogni confronto
#   a) results(): test di Wald trattato vs controllo sul modello unico.
#      independentFiltering = FALSE perché il filtro lo abbiamo deciso noi.
#   b) si tengono solo i geni del filtro di quel confronto e si ricalcola
#      il padj (Benjamini-Hochberg) su quei geni soltanto.
#   c) lfcShrink: log2FC "ristretto" verso 0 per i geni poco informativi.
#   d) UP / DOWN: padj < PADJ_CUT e |log2FC| >= LFC_CUT.
# -----------------------------------------------------------------------------
shrink <- function(trt, ref) {
  if (has_pkg("apeglm")) {
    d2 <- dds; d2$condition <- relevel(d2$condition, ref = ref)
    d2 <- nbinomWaldTest(d2, quiet = TRUE)
    s <- lfcShrink(d2, coef = paste0("condition_", trt, "_vs_", ref), type = "apeglm", quiet = TRUE)
    attr(s, "method") <- "apeglm"
  } else {
    s <- lfcShrink(dds, contrast = c("condition", trt, ref), type = "normal", quiet = TRUE)
    attr(s, "method") <- "normal"
  }
  s
}

summary_tab <- list()
nc <- counts(dds, normalized = TRUE)
for (nm in names(CONTRASTS)) {
  trt <- CONTRASTS[[nm]][1]; ref <- CONTRASTS[[nm]][2]
  padj_cut <- if (nm %in% names(PADJ_BY_CONTRAST)) PADJ_BY_CONTRAST[[nm]] else PADJ_CUT
  s_t <- coldata$sample[coldata$condition == trt]
  s_r <- coldata$sample[coldata$condition == ref]

  res <- results(dds, contrast = c("condition", trt, ref), independentFiltering = FALSE)  # (a)
  genes <- rownames(dds)[keep_by_contrast[rownames(dds), nm]]                             # (b)
  res <- res[genes, ]
  res$padj <- p.adjust(res$pvalue, method = "BH")
  shr <- shrink(trt, ref)[genes, ]                                                        # (c)

  lfc <- if (LFC_TYPE == "shrunk") shr$log2FoldChange else res$log2FoldChange
  status <- ifelse(res$padj < padj_cut & lfc >=  LFC_CUT, "UP",                           # (d)
            ifelse(res$padj < padj_cut & lfc <= -LFC_CUT, "DOWN", "NS"))
  status[is.na(status)] <- "NS"

  pct <- function(s) round(100 * rowSums(mod[genes, s]) / pmax(rowSums(tot[genes, s]), 1), 2)
  df <- data.frame(gene_id = genes, gene_name = gene_name[genes],
                   baseMean = res$baseMean,
                   log2FC_shrunk = shr$log2FoldChange, log2FC_MLE = res$log2FoldChange,
                   lfcSE = res$lfcSE, stat = res$stat, pvalue = res$pvalue, padj = res$padj,
                   status = status,
                   pct_mod_ctrl = pct(s_r), pct_mod_trt = pct(s_t),
                   round(nc[genes, c(s_r, s_t)], 2), check.names = FALSE)
  names(df)[names(df) == "pct_mod_ctrl"] <- paste0("pct_mod_", ref)
  names(df)[names(df) == "pct_mod_trt"]  <- paste0("pct_mod_", trt)
  df <- df[order(df$padj), ]
  write.csv(df, file.path(OUT_DIR, paste0(nm, "_all_genes.csv")), row.names = FALSE)
  write.csv(df[df$status != "NS", ], file.path(OUT_DIR, paste0(nm, "_significant.csv")), row.names = FALSE)

  v <- df; v$lfc <- if (LFC_TYPE == "shrunk") v$log2FC_shrunk else v$log2FC_MLE
  pv <- ggplot(v, aes(lfc, -log10(pmax(padj, 1e-300)), colour = status)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_colour_manual(values = c(DOWN = "#2166AC", NS = "grey75", UP = "#B2182B")) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = 2, linewidth = 0.3) +
    geom_hline(yintercept = -log10(padj_cut), linetype = 2, linewidth = 0.3) +
    labs(title = paste(trt, "vs", ref, "- RNA nascente"),
         subtitle = sprintf("UP %d  DOWN %d  (padj < %.2g)", sum(status == "UP"), sum(status == "DOWN"), padj_cut),
         x = paste0("log2FC (", LFC_TYPE, ")"), y = "-log10 padj", colour = NULL) + theme_bw()
  ggsave(file.path(OUT_DIR, paste0(nm, "_volcano.pdf")), pv, width = 7, height = 6)

  summary_tab[[nm]] <- data.frame(contrast = nm, genes_tested = length(genes),
                                  padj_cut = padj_cut,
                                  padj_sig = sum(res$padj < padj_cut, na.rm = TRUE),
                                  UP = sum(status == "UP"), DOWN = sum(status == "DOWN"),
                                  lfc = LFC_TYPE, shrinkage = attr(shr, "method"),
                                  size_factors = SF_MODE)
}

summary_df <- do.call(rbind, summary_tab)
write.csv(summary_df, file.path(OUT_DIR, "DEG_summary.csv"), row.names = FALSE)
print(summary_df, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
