#!/usr/bin/env Rscript
# =============================================================================
#  Arricchimento GO (BP) dei geni UP e DOWN dell'RNA-seq, per ogni confronto
#  - liste: rnaseq_lists/<lista>_UP.csv / _DOWN.csv (DESeq2, |log2FC| >= 0.58, padj < 0.05; NaCl padj < 0.1)
#  - universo: geni testati nel confronto con padj non NA (passati il filtro indipendente di DESeq2)
#  - ENSG -> Entrez con org.Hs.eg.db; annotazioni propagate (GOALL), termini 10-500 geni nell'universo
#  - ipergeometrico a una coda, BH per lista; riduzione ridondanza Jaccard >= 0.5 (come go_cluster.R)
#  Uso: Rscript 04_go_rnaseq_updown.R [cartella_output] [BP|MF|CC]
# =============================================================================
suppressPackageStartupMessages({ library(org.Hs.eg.db); library(GO.db); library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args) >= 1) args[1] else "results/rnaseq/GO_RNAseq_UP_DOWN"
ONT <- if (length(args) >= 2) args[2] else "BP"
MIN_SIZE <- 10; MAX_SIZE <- 500; PADJ_CUT <- 0.05; JACC_MAX <- 0.5; TOP <- 10
CONF <- list(Slow_vs_Amr = c("Slow_vs_Amr", "Slow_vs_Amr", "Slow vs Amr"),
             Fast_vs_Amr = c("Fast_vs_Amr", "Fast_vs_Amr", "Fast vs Amr"),
             Storm_vs_UT = c("STORM_vs_WT", "Storm_vs_UT", "Storm vs UT"),
             NaCl_vs_UT  = c("NACL_vs_WT",  "NaCl_vs_UT",  "NaCl vs UT"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
nv <- function(x) sub("\\.[0-9]+$", "", x)
term_all <- suppressMessages(AnnotationDbi::select(GO.db, keys = keys(GO.db), keytype = "GOID", columns = c("TERM", "ONTOLOGY")))
term_name <- setNames(term_all$TERM, term_all$GOID)

reduce <- function(res) {
  keep <- c(); kg <- list()
  for (i in seq_len(nrow(res))) {
    gi <- strsplit(res$geni[i], "/")[[1]]
    if (!any(vapply(kg, function(x) length(intersect(x, gi)) / length(union(x, gi)) >= JACC_MAX, logical(1)))) {
      keep <- c(keep, i); kg[[length(kg) + 1]] <- gi }
  }
  res[keep, ]
}
riep <- list(); tops <- list()
for (nm in names(CONF)) {
  lst <- CONF[[nm]][1]; tab <- CONF[[nm]][2]; tit <- CONF[[nm]][3]
  A <- read.csv(file.path(Sys.getenv("RNA_ALL_DIR", "results/rnaseq/RNAseq"), paste0(tab, "_all_genes.csv")))
  A <- A[!is.na(A$padj), ]; A$ens <- nv(A$gene_id)
  map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(A$ens), keytype = "ENSEMBL", columns = c("ENTREZID", "SYMBOL")))
  map <- map[!is.na(map$ENTREZID) & !duplicated(map$ENSEMBL), ]
  go <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(map$ENTREZID), keytype = "ENTREZID", columns = c("GOALL", "ONTOLOGYALL")))
  go <- unique(go[!is.na(go$GOALL) & go$ONTOLOGYALL == ONT, c("ENTREZID", "GOALL")])
  go$ENS <- map$ENSEMBL[match(go$ENTREZID, map$ENTREZID)]
  sets <- split(go$ENS, go$GOALL); sets <- sets[lengths(sets) >= MIN_SIZE & lengths(sets) <= MAX_SIZE]
  ann <- unique(go$ENS); N <- length(ann)
  sym <- setNames(ifelse(is.na(A$symbol) | A$symbol == "", A$ens, A$symbol), A$ens)
  for (dirz in c("UP", "DOWN")) {
    L <- read.csv(file.path(Sys.getenv("RNA_LISTS_DIR", "results/rnaseq/rnaseq_lists"), paste0(lst, "_", dirz, ".csv")))
    g <- intersect(nv(L[[1]]), ann); n <- length(g)
    res <- do.call(rbind, lapply(names(sets), function(id) {
      s <- sets[[id]]; hit <- g[g %in% s]; k <- length(hit); M <- length(s)
      if (k < 2) return(NULL)
      data.frame(GO = id, termine = term_name[[id]], geni_lista = k, dim_termine = M, GeneRatio = k / n, BgRatio = M / N,
                 fold = (k / n) / (M / N), pvalue = phyper(k - 1, M, N - M, n, lower.tail = FALSE),
                 geni = paste(sort(unname(sym[hit])), collapse = "/"))
    }))
    res$padj <- p.adjust(c(res$pvalue, rep(1, length(sets) - nrow(res))), "BH")[seq_len(nrow(res))]
    res <- res[order(res$pvalue), ]
    rs <- reduce(res[res$padj < PADJ_CUT, ])
    write.csv(res, file.path(OUT, sprintf("GO_%s_%s_%s_completo.csv", ONT, nm, dirz)), row.names = FALSE)
    write.csv(rs,  file.path(OUT, sprintf("GO_%s_%s_%s_ridotto.csv",  ONT, nm, dirz)), row.names = FALSE)
    cat(sprintf("%s %s: %d geni in lista (%d annotati); universo %d; termini padj<%.2f: %d (ridotti %d)\n",
                nm, dirz, nrow(L), n, N, PADJ_CUT, sum(res$padj < PADJ_CUT), nrow(rs)))
    riep[[length(riep) + 1]] <- data.frame(confronto = nm, direzione = dirz, geni_lista = nrow(L), geni_annotati = n,
                                           universo = N, termini_sig = sum(res$padj < PADJ_CUT), termini_ridotti = nrow(rs))
    if (nrow(rs)) tops[[paste(nm, dirz)]] <- cbind(head(rs, TOP), confronto = tit, direzione = dirz, n_lista = nrow(L))
  }
}
R <- do.call(rbind, riep); write.csv(R, file.path(OUT, "riepilogo_GO_RNAseq.csv"), row.names = FALSE)

# ---- dotplot: una figura per confronto (UP | DOWN) + una figura con tutti ----
wrap <- function(x) vapply(strwrap(x, 55, simplify = FALSE), paste, character(1), collapse = "\n")
plotta <- function(D, titolo, file) {
  D$col <- factor(sprintf("%s\n%s (n = %d)", D$confronto, D$direzione, D$n_lista),
                  levels = unique(sprintf("%s\n%s (n = %d)", D$confronto, D$direzione, D$n_lista)))
  D$termine_w <- wrap(D$termine)
  D$termine_w <- factor(D$termine_w, levels = unique(D$termine_w[order(as.integer(D$col), -D$pvalue)]))
  p <- ggplot(D, aes(col, termine_w)) + geom_point(aes(size = GeneRatio, colour = -log10(padj))) +
    scale_colour_gradient(low = "#9ec5f4", high = "#104281", name = expression(-log[10]~padj)) +
    scale_size_continuous(range = c(2, 7), name = "GeneRatio") +
    labs(title = titolo, subtitle = sprintf("geni UP e DOWN, top %d termini non ridondanti\nuniverso = geni testati in DESeq2 (padj non NA)", TOP),
         x = NULL, y = NULL) +
    theme_bw(base_size = 10) + theme(axis.text.y = element_text(size = 8), plot.title = element_text(face = "bold"),
                                     plot.subtitle = element_text(size = 8, colour = "#444444"), panel.grid.major = element_line(colour = "#ececec"))
  h <- max(4.5, 0.3 * nlevels(D$termine_w) + 2)
  ggsave(paste0(file, ".pdf"), p, width = 4 + 1.6 * nlevels(D$col), height = h, limitsize = FALSE)
  ggsave(paste0(file, ".png"), p, width = 4 + 1.6 * nlevels(D$col), height = h, dpi = 220, limitsize = FALSE, bg = "white")
}
for (nm in names(CONF)) {
  D <- do.call(rbind, tops[grep(paste0("^", nm, " "), names(tops))])
  if (!is.null(D)) plotta(D, sprintf("GO %s - RNA-seq %s", ONT, CONF[[nm]][3]),
                          file.path(OUT, sprintf("GO_%s_%s_dotplot", ONT, nm)))
}
print(R)
