#!/usr/bin/env Rscript
# =============================================================================
#  Integrazione RNA-seq x RNA nascente
#   - Venn dei geni significativi (UP + DOWN) per ogni confronto
#   - matrici 3x3 di conteggio RNA-seq (UP / NC / DOWN) x nascente (UP / NC / DOWN)
#
#  Uso: Rscript 01_integrazione_rnaseq_nascent.R <cartella_RNAseq> <cartella_nascent> <cartella_output>
#   <cartella_RNAseq>  = liste RNA-seq <confronto>_UP.csv / _DOWN.csv
#   <cartella_nascent> = output di 02_deseq2_nascent.R (<confronto>_all_genes.csv)
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

args    <- commandArgs(trailingOnly = TRUE)
RNA_DIR <- if (length(args) >= 1) args[1] else "RNAseq"
NAS_DIR <- if (length(args) >= 2) args[2] else "nascent"
OUT_DIR <- if (length(args) >= 3) args[3] else "integrazione"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CONTRASTS <- c("Fast_vs_Amr", "Slow_vs_Amr", "Storm_vs_UT", "NaCl_vs_UT")
COL_RNA   <- "#2a78d6"   # RNA-seq
COL_NAS   <- "#eb6834"   # RNA nascente

# Le classi UP/DOWN vengono dalle liste RNA-seq e dalla colonna "status" del nascente
# (padj < 0.05, NaCl < 0.10; |log2FC| >= 0.58).
no_ver <- function(x) sub("\\.[0-9]+$", "", x)

# RNA-seq: liste UP / DOWN già filtrate (<prefisso>_UP.csv, <prefisso>_DOWN.csv)
RNA_FILES <- c(Fast_vs_Amr = "Fast_vs_Amr", Slow_vs_Amr = "Slow_vs_Amr",
               Storm_vs_UT = "STORM_vs_WT", NaCl_vs_UT = "NACL_vs_WT")
load_rna <- function(nm) {
  rd <- function(cl) {
    d <- read.csv(file.path(RNA_DIR, paste0(RNA_FILES[[nm]], "_", cl, ".csv")), row.names = 1)
    data.frame(gene = no_ver(rownames(d)), lfc_rna = d$log2FoldChange,
               padj_rna = d$padj, class_rna = cl)
  }
  rbind(rd("UP"), rd("DOWN"))
}
load_nas <- function(nm) {
  d <- read.csv(file.path(NAS_DIR, paste0(nm, "_all_genes.csv")), stringsAsFactors = FALSE)
  d <- d[!is.na(d$padj), ]
  data.frame(gene = no_ver(d$gene_id), gene_name = d$gene_name,
             lfc_nas = d$log2FC_shrunk, padj_nas = d$padj,
             class_nas = d$status)
}

# --- Venn a due insiemi disegnato con ggplot ---------------------------------
circle <- function(x0, r = 1, n = 200) {
  t <- seq(0, 2 * pi, length.out = n); data.frame(x = x0 + r * cos(t), y = r * sin(t))
}
venn_plot <- function(n_rna, n_nas, n_both, title) {
  ggplot() +
    geom_polygon(data = circle(-0.55), aes(x, y), fill = COL_RNA, alpha = 0.25,
                 colour = COL_RNA, linewidth = 0.8) +
    geom_polygon(data = circle(0.55), aes(x, y), fill = COL_NAS, alpha = 0.25,
                 colour = COL_NAS, linewidth = 0.8) +
    annotate("text", x = -1.05, y = 0, label = n_rna - n_both, size = 4.2) +
    annotate("text", x = 0,     y = 0, label = n_both, size = 4.6, fontface = "bold") +
    annotate("text", x = 1.05,  y = 0, label = n_nas - n_both, size = 4.2) +
    annotate("text", x = -0.9, y = 1.22, label = paste0("RNA-seq (", n_rna, ")"),
             colour = "#0b0b0b", size = 3.4, fontface = "bold") +
    annotate("text", x = 0.9, y = 1.22, label = paste0("Nascente (", n_nas, ")"),
             colour = "#0b0b0b", size = 3.4, fontface = "bold") +
    coord_equal(xlim = c(-1.75, 1.75), ylim = c(-1.1, 1.4)) +
    labs(title = title) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
}

# --- matrice 2x2 --------------------------------------------------------------
matrix_plot <- function(tab, title) {
  df <- as.data.frame(as.table(tab)); names(df) <- c("rna", "nas", "n")
  df$rna <- factor(df$rna, levels = c("DOWN", "UP"))
  df$nas <- factor(df$nas, levels = c("UP", "DOWN"))
  df$pct <- 100 * df$n / ave(df$n, df$rna, FUN = sum)        # % sulla riga RNA-seq
  df$lab <- sprintf("%d\n(%.0f%%)", df$n, df$pct)
  df$ink <- ifelse(df$pct > 55, "white", "#0b0b0b")
  ggplot(df, aes(nas, rna, fill = pct)) +
    geom_tile(colour = "#fcfcfb", linewidth = 1.5) +
    geom_text(aes(label = lab, colour = ink), size = 3.3, lineheight = 0.9) +
    scale_colour_identity() +
    scale_fill_gradient(low = "#cde2fb", high = "#184f95", limits = c(0, 100),
                        name = "% della riga") +
    labs(title = title, x = "RNA nascente", y = "RNA-seq") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
}

venns <- list(); mats <- list(); summary <- list()
for (nm in CONTRASTS) {
  r <- load_rna(nm); s <- load_nas(nm)
  m <- merge(r, s, by = "gene", all = TRUE)
  m$tested_nascent <- !is.na(m$class_nas)
  m$class_nas[is.na(m$class_nas)] <- "non testato"
  m$class_rna[is.na(m$class_rna)] <- "NS"
  write.csv(m[order(m$class_rna, m$class_nas), ],
            file.path(OUT_DIR, paste0("integrazione_", nm, ".csv")), row.names = FALSE)

  # Venn: tutti i geni delle liste RNA-seq x tutti i significativi del nascente
  sig_r <- m$class_rna %in% c("UP", "DOWN")
  sig_n <- m$class_nas %in% c("UP", "DOWN")
  venns[[nm]] <- venn_plot(sum(sig_r), sum(sig_n), sum(sig_r & sig_n), gsub("_", " ", nm))

  # matrice 2x2: geni significativi (UP o DOWN) in ENTRAMBE le analisi
  k <- sig_r & sig_n
  tab <- table(factor(m$class_rna[k], levels = c("UP", "DOWN")),
               factor(m$class_nas[k], levels = c("UP", "DOWN")))
  write.csv(as.data.frame.matrix(tab), file.path(OUT_DIR, paste0("matrice_", nm, ".csv")))
  mats[[nm]] <- matrix_plot(tab, gsub("_", " ", nm))

  summary[[nm]] <- data.frame(
    confronto = nm, RNAseq_UP = sum(m$class_rna == "UP"), RNAseq_DOWN = sum(m$class_rna == "DOWN"),
    RNAseq_non_testati_nascente = sum(sig_r & !m$tested_nascent),
    nascente_UP = sum(m$class_nas == "UP"), nascente_DOWN = sum(m$class_nas == "DOWN"),
    in_comune = sum(k), UP_UP = tab["UP", "UP"], UP_DOWN = tab["UP", "DOWN"],
    DOWN_UP = tab["DOWN", "UP"], DOWN_DOWN = tab["DOWN", "DOWN"])
}

sm <- do.call(rbind, summary)
write.csv(sm, file.path(OUT_DIR, "riepilogo_integrazione.csv"), row.names = FALSE)
print(sm, row.names = FALSE)

fig_v <- wrap_plots(venns, nrow = 1) +
  plot_annotation(title = "Geni significativi (UP + DOWN): RNA-seq x RNA nascente",
                  subtitle = "|log2FC| >= 0.58; padj < 0.05 (NaCl < 0.10)")
ggsave(file.path(OUT_DIR, "venn_RNAseq_nascente.pdf"), fig_v, width = 16, height = 4.6)
ggsave(file.path(OUT_DIR, "venn_RNAseq_nascente.png"), fig_v, width = 16, height = 4.6, dpi = 200, bg = "white")

fig_m <- wrap_plots(mats, nrow = 1) +
  plot_annotation(title = "Geni significativi in entrambe le analisi: direzione RNA-seq x RNA nascente",
                  subtitle = "|log2FC| >= 0.58; padj < 0.05 (NaCl < 0.10). Percentuale sulla riga RNA-seq.")
ggsave(file.path(OUT_DIR, "matrici_RNAseq_nascente.pdf"), fig_m, width = 14, height = 3.8)
ggsave(file.path(OUT_DIR, "matrici_RNAseq_nascente.png"), fig_m, width = 14, height = 3.8, dpi = 200, bg = "white")
