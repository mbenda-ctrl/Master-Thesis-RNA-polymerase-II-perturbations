#!/usr/bin/env Rscript
# =============================================================================
#  Volcano plot RNA-seq con i geni più significativi etichettati (gene symbol)
#
#  - Punti: tutti i geni testati (tabella DESeq2 completa, <confronto>_all_genes.csv)
#  - UP / DOWN: geni delle liste RNA-seq usate nelle analisi (<lista>_UP.csv / _DOWN.csv);
#    per questi geni log2FC e padj sono quelli della lista
#  - Etichette: i TOP_N UP e TOP_N DOWN con padj più basso (a parità: |log2FC| più alto)
#
#  Uso: Rscript 03_volcano_rnaseq.R [cartella_DESeq2] [cartella_liste] [cartella_output]
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(ggrepel); library(patchwork) })

args    <- commandArgs(trailingOnly = TRUE)
ALL_DIR <- if (length(args) >= 1) args[1] else "RNAseq"
LST_DIR <- if (length(args) >= 2) args[2] else "rnaseq_lists"
OUT_DIR <- if (length(args) >= 3) args[3] else "volcano_RNAseq"
RNAKINET <- Sys.getenv("RNAKINET_DIR", "results/nascent/rnakinet_genes")   # simboli di riserva (stessa annotazione Ensembl)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TOP_N   <- 10      # geni etichettati per direzione
LFC_CUT <- 0.58
CONTRASTS <- list(   # nome = c(tabella completa, prefisso liste, soglia padj, titolo)
  Fast_vs_Amr = c("Fast_vs_Amr", "Fast_vs_Amr", 0.05, "Fast vs Amr"),
  Slow_vs_Amr = c("Slow_vs_Amr", "Slow_vs_Amr", 0.05, "Slow vs Amr"),
  Storm_vs_UT = c("Storm_vs_UT", "STORM_vs_WT", 0.05, "Storm vs UT"),
  NaCl_vs_UT  = c("NaCl_vs_UT",  "NACL_vs_WT",  0.10, "NaCl vs UT"))
COL <- c(UP = "#b2182b", DOWN = "#2166ac", NS = "#c8c8c8")

nv <- function(x) sub("\\.[0-9]+$", "", x)
rk_files <- list.files(RNAKINET, "rnakinet.*csv$", full.names = TRUE)
map <- if (length(rk_files)) do.call(rbind, lapply(rk_files, function(f) read.csv(f)[, c("gene_id", "gene_name")])) else
         data.frame(gene_id = character(), gene_name = character())
map$gene_id <- nv(map$gene_id); map <- map[!duplicated(map$gene_id), ]
sym_rk <- setNames(map$gene_name, map$gene_id)

plots <- list(); summ <- list()
for (nm in names(CONTRASTS)) {
  cc <- CONTRASTS[[nm]]; padj_cut <- as.numeric(cc[3])
  a <- read.csv(file.path(ALL_DIR, paste0(cc[1], "_all_genes.csv")), stringsAsFactors = FALSE)
  a <- data.frame(id = nv(a$gene_id), symbol = a$symbol, lfc = a$log2FC_MLE, padj = a$padj)
  rd <- function(cl) { d <- read.csv(file.path(LST_DIR, paste0(cc[2], "_", cl, ".csv")), row.names = 1)
                       data.frame(id = nv(rownames(d)), lfc_l = d$log2FoldChange, padj_l = d$padj, cls = cl) }
  L <- rbind(rd("UP"), rd("DOWN"))
  d <- merge(a, L, by = "id", all = TRUE)
  d$lfc  <- ifelse(is.na(d$lfc_l), d$lfc, d$lfc_l)          # valori della lista per i geni in lista
  d$padj <- ifelse(is.na(d$padj_l), d$padj, d$padj_l)
  d$cls[is.na(d$cls)] <- "NS"
  d <- d[!is.na(d$padj) & !is.na(d$lfc), ]
  d$symbol <- ifelse(is.na(d$symbol) | d$symbol == "", sym_rk[d$id], d$symbol)
  d$symbol <- ifelse(is.na(d$symbol) | d$symbol == "", d$id, d$symbol)
  ymin <- min(d$padj[d$padj > 0]); d$y <- -log10(pmax(d$padj, ymin))   # padj = 0 -> limite numerico
  d$cls <- factor(d$cls, levels = c("NS", "DOWN", "UP"))
  d <- d[order(d$cls), ]                                                # NS disegnati sotto

  top <- do.call(rbind, lapply(c("UP", "DOWN"), function(k) {
    x <- d[d$cls == k, ]; head(x[order(x$padj, -abs(x$lfc)), ], TOP_N) }))
  write.csv(top[, c("id", "symbol", "cls", "lfc", "padj")], file.path(OUT_DIR, paste0("etichettati_", nm, ".csv")), row.names = FALSE)

  n_up <- sum(d$cls == "UP"); n_dn <- sum(d$cls == "DOWN")
  p <- ggplot(d, aes(lfc, y, colour = cls)) +
    geom_point(size = 0.7, alpha = 0.7, stroke = 0) +
    scale_colour_manual(values = COL, breaks = c("UP", "DOWN"),
                        labels = c(paste0("UP (", n_up, ")"), paste0("DOWN (", n_dn, ")")), name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = 2, linewidth = 0.3, colour = "#666666") +
    geom_hline(yintercept = -log10(padj_cut), linetype = 2, linewidth = 0.3, colour = "#666666") +
    geom_text_repel(data = top, aes(label = symbol), colour = "#111111", size = 2.9, fontface = "italic",
                    max.overlaps = Inf, min.segment.length = 0, segment.size = 0.25, segment.colour = "#555555",
                    box.padding = 0.35, seed = 1) +
    coord_cartesian(xlim = c(-1, 1) * min(max(abs(d$lfc), na.rm = TRUE), 10)) +
    labs(title = cc[4], subtitle = sprintf("padj < %.2g, |log2FC| >= %.2f", padj_cut, LFC_CUT),
         x = "log2 fold change", y = expression(-log[10]~padj)) +
    theme_classic(base_size = 11) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, colour = "#444444"))
  ggsave(file.path(OUT_DIR, paste0("volcano_", nm, ".pdf")), p, width = 6, height = 6)
  ggsave(file.path(OUT_DIR, paste0("volcano_", nm, ".png")), p, width = 6, height = 6, dpi = 300)
  plots[[nm]] <- p
  summ[[nm]] <- data.frame(confronto = nm, geni = nrow(d), UP = n_up, DOWN = n_dn)
}
fig <- wrap_plots(plots, nrow = 2)
ggsave(file.path(OUT_DIR, "volcano_tutti.pdf"), fig, width = 12, height = 12)
ggsave(file.path(OUT_DIR, "volcano_tutti.png"), fig, width = 12, height = 12, dpi = 200, bg = "white")
print(do.call(rbind, summ), row.names = FALSE)
