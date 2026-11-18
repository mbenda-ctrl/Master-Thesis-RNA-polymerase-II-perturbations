#!/usr/bin/env Rscript
# =============================================================================
#  Venn a 4 insiemi: m6A (diff. metilati, meno siti, siti persi) x Pol II DiffBind
#
#  - m6A: fogli <confronto>_diff_metilati / _meno_siti / _persi del file Excel
#  - Pol II: geni del BED (corpo genico +/- 1 kb) con regioni DiffBind:
#       Slow vs Amr  FDR < 0.05
#       NaCl vs WT   ChIP KCl vs UT, FDR < 0.01 e |Fold| > 1.5
#       Fast, Storm  tutti i siti (nessuna regione significativa)
#  - Si lavora in ENSG (senza versione); i simboli del BED -> ENSG via file RNAkinet
#  - Test di Fisher: tra i geni del BED, i geni m6A sono più spesso Pol II?
#
#  Uso: Rscript 03_m6a_polII_venn.R   (dalla cartella final_sep)
# =============================================================================
suppressPackageStartupMessages({ library(GenomicRanges); library(ggplot2); library(patchwork); library(readxl) })

XLSX <- "m6a/geni_diff_metilati_e_perdita_siti.xlsx"
BED  <- "chip/geni_intersezione_collega_MANE.bed"
RNAKINET <- "../nascent/input"
OUT  <- "m6a_polII"; FLANK <- 1000
dir.create(OUT, showWarnings = FALSE)

SETUP <- list(
  Slow_vs_Amr = list(m6a = "Slow_vs_Amr", diffbind = "chip/Slow_vs_Amr_DiffBind_FULL.tsv", fdr = 0.05, fold = 0,
                     lab = "Pol II Slow vs Amr (FDR < 0.05)"),
  NaCl_vs_WT  = list(m6a = "NaCl_vs_WT", diffbind = "chip/KCl_vs_UT_DiffBind_all_sites.tsv", fdr = 0.01, fold = 1.5,
                     lab = "Pol II KCl vs UT (FDR < 0.01, |Fold| > 1.5)"),
  Fast_vs_Amr = list(m6a = "Fast_vs_Amr", diffbind = "chip/Fast_vs_Amr_DiffBind_FULL.tsv", fdr = 1, fold = 0,
                     lab = "Pol II Fast vs Amr (tutti i siti)"),
  Storm_vs_WT = list(m6a = "Storm_vs_WT", diffbind = "chip/Storm_vs_UT_DiffBind_FULL.tsv", fdr = 1, fold = 0,
                     lab = "Pol II Storm vs UT (tutti i siti)"))
COL <- c("#2a78d6", "#eb6834", "#1baf7a", "#6250d6")   # diff. metilati, meno siti, persi, Pol II

# ---- geni del BED e mappa simbolo -> ENSG ---------------------------------------
nv  <- function(x) sub("\\.[0-9]+$", "", x)
bed <- read.delim(BED, header = FALSE, col.names = c("chr", "start", "end", "symbol"))
genes <- resize(GRanges(bed$chr, IRanges(bed$start + 1, bed$end), symbol = bed$symbol),
                width = bed$end - bed$start + 2 * FLANK, fix = "center")
map <- do.call(rbind, lapply(list.files(RNAKINET, "rnakinet.*csv$", full.names = TRUE),
                             function(f) read.csv(f)[, c("gene_id", "gene_name")]))
map$gene_id <- nv(map$gene_id); map <- map[!duplicated(map$gene_id), ]
ensg_of <- split(map$gene_id, map$gene_name)
to_ensg <- function(sy) unique(unlist(ensg_of[intersect(sy, names(ensg_of))]))
bed_ensg <- to_ensg(bed$symbol)

polII_genes <- function(cfg) {
  d <- read.delim(cfg$diffbind)
  d <- d[d$FDR <= cfg$fdr & abs(d$Fold) >= cfg$fold, ]
  ov <- findOverlaps(genes, GRanges(d$seqnames, IRanges(d$start, d$end)))
  to_ensg(unique(genes$symbol[queryHits(ov)]))
}

# ---- Venn a 4 ellissi ------------------------------------------------------------
E <- data.frame(x = c(0.38, 0.47, 0.53, 0.62), y = c(0.45, 0.53, 0.53, 0.45),
                a = 0.33, b = 0.18, th = c(-40, -40, 40, 40) * pi / 180)
inside <- function(px, py, e) { dx <- px - e$x; dy <- py - e$y
  u <- dx * cos(e$th) + dy * sin(e$th); v <- -dx * sin(e$th) + dy * cos(e$th); (u / e$a)^2 + (v / e$b)^2 <= 1 }
grid <- expand.grid(x = seq(0, 1, 0.003), y = seq(0, 1, 0.003))
memb <- sapply(1:4, function(i) inside(grid$x, grid$y, E[i, ]))
code <- apply(memb, 1, function(r) paste(as.integer(r), collapse = ""))
pos  <- aggregate(cbind(x, y) ~ code, data = cbind(grid, code)[code != "0000", ], FUN = median)
ell  <- function(i, n = 300) { t <- seq(0, 2 * pi, length.out = n); e <- E[i, ]
  data.frame(x = e$x + e$a * cos(t) * cos(e$th) - e$b * sin(t) * sin(e$th),
             y = e$y + e$a * cos(t) * sin(e$th) + e$b * sin(t) * cos(e$th), id = i) }

venn4 <- function(sets, labs, title) {
  all <- unique(unlist(sets))
  m <- sapply(sets, function(s) all %in% s)
  cc <- apply(m, 1, function(r) paste(as.integer(r), collapse = ""))
  pos$n <- as.integer(table(factor(cc, levels = pos$code)))
  shapes <- do.call(rbind, lapply(1:4, ell)); shapes$id <- factor(shapes$id)
  lab_df <- data.frame(x = c(0.10, 0.30, 0.70, 0.90), y = c(0.80, 0.93, 0.93, 0.80),
                       lab = sprintf("%s\n(%d)", labs, lengths(sets)), hj = c(0, 0.5, 0.5, 1))
  ggplot() +
    geom_polygon(data = shapes, aes(x, y, group = id, fill = id, colour = id), alpha = 0.15, linewidth = 0.7) +
    scale_fill_manual(values = COL, guide = "none") + scale_colour_manual(values = COL, guide = "none") +
    geom_text(data = pos, aes(x, y, label = n), size = 3.2,
              fontface = ifelse(pos$code == "1111", "bold", "plain")) +
    geom_point(data = lab_df, aes(x = x + c(-0.025, -0.11, -0.11, 0.025), y = y + 0.012), colour = COL, size = 3.2) +
    geom_text(data = lab_df, aes(x, y, label = lab, hjust = hj), size = 3, fontface = "bold", colour = "#0b0b0b") +
    coord_equal(xlim = c(0, 1), ylim = c(0.1, 1)) + labs(title = title) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
}

# ---- per confronto ----------------------------------------------------------------
plots <- list(); summ <- list(); enr <- list()
for (nm in names(SETUP)) {
  cfg <- SETUP[[nm]]
  rd <- function(s) read_excel(XLSX, sheet = paste0(cfg$m6a, "_", s))
  dm <- rd("diff_metilati"); ms <- rd("meno_siti"); ps <- rd("persi")
  sets <- list(diff_metilati = unique(nv(dm$gene_id)), meno_siti = unique(nv(ms$gene_id)),
               persi = unique(nv(ps$gene_id)), PolII = polII_genes(cfg))
  plots[[nm]] <- venn4(sets, c("m6A diff. metilati", "m6A meno siti", "m6A siti persi", cfg$lab),
                       gsub("_", " ", nm))

  # tabella per gene
  ids <- unique(unlist(sets))
  sym <- setNames(c(dm$gene_name, ms$gene_name, ps$gene_name, map$gene_name),
                  c(nv(dm$gene_id), nv(ms$gene_id), nv(ps$gene_id), map$gene_id))
  t <- data.frame(gene_id = ids, gene_name = unname(sym[ids]),
                  sapply(sets, function(s) ids %in% s), in_BED = ids %in% bed_ensg)
  t$direzione_stech <- dm$direzione[match(ids, nv(dm$gene_id))]
  t <- t[order(-rowSums(t[, names(sets)])), ]
  write.csv(t, file.path(OUT, paste0("geni_", nm, "_m6A_PolII.csv")), row.names = FALSE)

  # test di Fisher sui geni del BED (solo qui un gene può essere Pol II)
  for (s in c("diff_metilati", "meno_siti", "persi")) {
    g <- intersect(sets[[s]], bed_ensg); o <- setdiff(bed_ensg, g)
    tab <- matrix(c(sum(g %in% sets$PolII), sum(!g %in% sets$PolII),
                    sum(o %in% sets$PolII), sum(!o %in% sets$PolII)), 2)
    ft <- fisher.test(tab)
    enr[[length(enr) + 1]] <- data.frame(confronto = nm, insieme_m6A = s, n_geni = length(sets[[s]]),
      n_nel_BED = length(g), pct_PolII = round(100 * mean(g %in% sets$PolII), 1),
      pct_PolII_altri_BED = round(100 * mean(o %in% sets$PolII), 1),
      odds_ratio = round(unname(ft$estimate), 2), p_value = signif(ft$p.value, 3))
  }
  summ[[nm]] <- data.frame(confronto = nm, diff_metilati = length(sets[[1]]), meno_siti = length(sets[[2]]),
    persi = length(sets[[3]]), PolII = length(sets[[4]]),
    m6A_qualsiasi_e_PolII = length(intersect(unique(unlist(sets[1:3])), sets$PolII)),
    tutti_e_quattro = length(Reduce(intersect, sets)))
}
S <- do.call(rbind, summ); F <- do.call(rbind, enr)
write.csv(S, file.path(OUT, "riepilogo_m6A_PolII.csv"), row.names = FALSE)
write.csv(F, file.path(OUT, "arricchimento_m6A_PolII_Fisher.csv"), row.names = FALSE)
print(S, row.names = FALSE); print(F, row.names = FALSE)

fig <- wrap_plots(plots, nrow = 2) +
  plot_annotation(title = "m6A (diff. metilati, meno siti, siti persi) x Pol II DiffBind",
                  subtitle = "Pol II: geni del BED (corpo genico +/- 1 kb) con regioni DiffBind; NaCl incrociato con la ChIP KCl")
ggsave(file.path(OUT, "venn_m6A_PolII.pdf"), fig, width = 13, height = 11)
ggsave(file.path(OUT, "venn_m6A_PolII.png"), fig, width = 13, height = 11, dpi = 200, bg = "white")

# ---- Venn a 2 cerchi separati: ogni lista m6A contro la Pol II ------------------
circ <- function(x0, r = 1, n = 200) { t <- seq(0, 2 * pi, length.out = n); data.frame(x = x0 + r * cos(t), y = r * sin(t)) }
venn2 <- function(A, B, labA, labB, colA, title) {
  nAB <- length(intersect(A, B))
  ggplot() +
    geom_polygon(data = circ(-0.55), aes(x, y), fill = colA, alpha = 0.2, colour = colA, linewidth = 0.8) +
    geom_polygon(data = circ(0.55), aes(x, y), fill = COL[4], alpha = 0.2, colour = COL[4], linewidth = 0.8) +
    annotate("text", x = -1.05, y = 0, label = length(A) - nAB, size = 4) +
    annotate("text", x = 0, y = 0, label = nAB, size = 4.4, fontface = "bold") +
    annotate("text", x = 1.05, y = 0, label = length(B) - nAB, size = 4) +
    annotate("text", x = -0.85, y = 1.28, label = sprintf("%s\n(%d)", labA, length(A)), size = 3.1, fontface = "bold", colour = "#0b0b0b") +
    annotate("text", x = 0.85, y = 1.28, label = sprintf("%s\n(%d)", labB, length(B)), size = 3.1, fontface = "bold", colour = "#0b0b0b") +
    coord_equal(xlim = c(-1.75, 1.75), ylim = c(-1.1, 1.55)) + labs(title = title) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10))
}
NOMI <- c(diff_metilati = "m6A diff. metilati", meno_siti = "m6A meno siti", persi = "m6A siti persi")
p2 <- list()
for (nm in names(SETUP)) {
  cfg <- SETUP[[nm]]
  rd <- function(s) read_excel(XLSX, sheet = paste0(cfg$m6a, "_", s))
  sets <- list(diff_metilati = unique(nv(rd("diff_metilati")$gene_id)), meno_siti = unique(nv(rd("meno_siti")$gene_id)),
               persi = unique(nv(rd("persi")$gene_id)))
  pol <- polII_genes(cfg)
  polab <- sub(" \\(", "\n(", sub("Pol II [A-Za-z]+ vs [A-Za-z]+ ", "Pol II ", cfg$lab))
  for (i in seq_along(sets)) {
    s <- names(sets)[i]
    fr <- F[F$confronto == nm & F$insieme_m6A == s, ]
    p2[[paste(nm, s)]] <- venn2(sets[[s]], pol, NOMI[[s]], "Pol II", COL[i],
      sprintf("%s - %s\nFisher (geni del BED): OR %.1f, p = %.1e", gsub("_", " ", nm), NOMI[[s]], fr$odds_ratio, fr$p_value))
  }
}
row_lab <- sapply(SETUP, function(x) x$lab)
fig2 <- wrap_plots(p2, ncol = 3) +
  plot_annotation(title = "m6A x Pol II DiffBind: un Venn per ogni lista m6A",
                  subtitle = paste("Pol II:", paste(sprintf("%s = %s", gsub("_", " ", names(row_lab)), sub("Pol II [A-Za-z]+ vs [A-Za-z]+ ", "", row_lab)), collapse = "; ")))
ggsave(file.path(OUT, "venn_m6A_PolII_separati.pdf"), fig2, width = 13, height = 17)
ggsave(file.path(OUT, "venn_m6A_PolII_separati.png"), fig2, width = 13, height = 17, dpi = 180, bg = "white")
for (nm in names(SETUP)) {                     # anche una figura per condizione
  f <- wrap_plots(p2[grep(paste0("^", nm, " "), names(p2))], ncol = 3) +
    plot_annotation(title = paste("m6A x Pol II -", gsub("_", " ", nm)), subtitle = SETUP[[nm]]$lab)
  ggsave(file.path(OUT, paste0("venn_m6A_PolII_", nm, ".png")), f, width = 13, height = 4.8, dpi = 220, bg = "white")
  ggsave(file.path(OUT, paste0("venn_m6A_PolII_", nm, ".pdf")), f, width = 13, height = 4.8)
}

# ---- solo i diff. metilati (4 condizioni, griglia 2x2) -----------------------------
p3 <- p2[grep(" diff_metilati$", names(p2))]
fig3 <- wrap_plots(p3, ncol = 2) +
  plot_annotation(title = "m6A diff. metilati x Pol II DiffBind",
                  subtitle = paste("Pol II:", paste(sprintf("%s = %s", gsub("_", " ", names(row_lab)), sub("Pol II [A-Za-z]+ vs [A-Za-z]+ ", "", row_lab)), collapse = "; ")))
ggsave(file.path(OUT, "venn_m6A_diffmetilati_PolII.pdf"), fig3, width = 10, height = 9)
ggsave(file.path(OUT, "venn_m6A_diffmetilati_PolII.png"), fig3, width = 10, height = 9, dpi = 220, bg = "white")
