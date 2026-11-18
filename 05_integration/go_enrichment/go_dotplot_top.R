#!/usr/bin/env Rscript
# Ridisegna il dotplot GO dai file *_ridotto.csv di go_cluster.R con un numero di termini a scelta (nessun ricalcolo)
# Uso: Rscript go_dotplot_top.R <cartella_GO> <tabella_heatmap.tsv> <TOP> [BP] [etichetta]
suppressPackageStartupMessages(library(ggplot2))
a <- commandArgs(trailingOnly = TRUE)
DIR <- a[1]; HEAT <- a[2]; TOP <- as.integer(a[3]); ONT <- if (length(a) >= 4) a[4] else "BP"
LAB <- if (length(a) >= 5) a[5] else sub("heatmap_(rna|chip)_(.*)\\.tsv", "\\2", basename(HEAT))
H <- read.delim(HEAT, row.names = 1); ks <- sort(unique(H$gruppo))
lab <- sapply(ks, function(k) sprintf("Gruppo %d\n(n = %d)", k, sum(H$gruppo == k)))
D <- do.call(rbind, lapply(ks, function(k) {
  r <- read.csv(file.path(DIR, sprintf("GO_%s_gruppo%d_ridotto.csv", ONT, k)))
  if (nrow(r)) cbind(head(r, TOP), gr = k) else NULL }))
D$gruppo_lab <- factor(lab[match(D$gr, ks)], levels = lab)
wrap <- function(x) vapply(strwrap(x, 60, simplify = FALSE), paste, character(1), collapse = "\n")
D$termine_w <- wrap(D$termine)
D$termine_w <- factor(D$termine_w, levels = unique(D$termine_w[order(D$gr, -D$pvalue)]))
p <- ggplot(D, aes(gruppo_lab, termine_w)) + scale_x_discrete(drop = FALSE) +
  geom_point(aes(size = GeneRatio, colour = -log10(padj))) +
  scale_colour_gradient(low = "#9ec5f4", high = "#104281", name = expression(-log[10]~padj)) +
  scale_size_continuous(range = c(1.8, 6.5), name = "GeneRatio") +
  labs(title = sprintf("Arricchimento GO %s per gruppo", ONT),
       subtitle = sprintf("%s - universo: genoma; primi %d termini non ridondanti per gruppo (padj < 0.05)", LAB, TOP), x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major = element_line(colour = "#ececec"), axis.text.y = element_text(size = 7.5),
        plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 8, colour = "#444444"))
h <- max(5, 0.26 * nlevels(D$termine_w) + 2)
f <- file.path(DIR, sprintf("GO_%s_dotplot_top%d", ONT, TOP))
ggsave(paste0(f, ".pdf"), p, width = 9.5, height = h, limitsize = FALSE)
ggsave(paste0(f, ".png"), p, width = 9.5, height = h, dpi = 250, limitsize = FALSE, bg = "white")
cat(table(D$gr), "termini per gruppo\n")
