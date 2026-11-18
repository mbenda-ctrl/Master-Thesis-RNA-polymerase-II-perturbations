#!/usr/bin/env Rscript
# =============================================================================
#  Integrazione ChIP (DiffBind) x RNA-seq x RNA nascente
#
#  1. BED: geni su cui si assegnano le regioni ChIP (geni_intersezione_MANE.bed, 01_chipseq/12b)
#  2. Gene "ChIP" = gene (corpo genico +/- FLANK) che si sovrappone ad almeno una
#     regione DiffBind con FDR <= fdr del confronto (1 = tutti i siti). Gain = Fold > 0, Loss = Fold < 0.
#  3. RNA-seq: liste UP/DOWN complete; nascente: colonna status, complete.
#     Si lavora in ENSG; i simboli del BED -> ENSG tramite i file RNAkinet.
#  4. Per ogni confronto: Venn a 3 insiemi per gli UP e per i DOWN
#     (RNA-seq, nascente, ChIP), più tabelle per gene.
#
#  Uso: Rscript 02_chip_integrazione.R [confronti separati da virgola] [cartella_output]
#       es. Rscript 02_chip_integrazione.R Slow_vs_Amr,Fast_vs_Amr,Storm_vs_UT results/integration/chip
#  Percorsi di input modificabili con variabili d'ambiente:
#       GENE_BED, DIFFBIND_DIR, RNAKINET_DIR, RNA_LISTS_DIR, NASCENT_DIR
# =============================================================================
suppressPackageStartupMessages({ library(GenomicRanges); library(ggplot2); library(patchwork) })

# ---- PASSO 0: parametri e file ----------------------------------------------
env <- function(k, d) { v <- Sys.getenv(k); if (nzchar(v)) v else d }
BED       <- env("GENE_BED",      "results/chipseq/gene_counts/geni_intersezione_MANE.bed")
DB_DIR    <- env("DIFFBIND_DIR",  "results/chipseq/diffbind")      # <confronto>_DiffBind_all_sites.tsv (01_chipseq/10)
RNAKINET  <- env("RNAKINET_DIR",  "results/nascent/rnakinet_genes") # per la mappa ENSG -> simbolo
RNA_DIR   <- env("RNA_LISTS_DIR", "results/rnaseq/rnaseq_lists")   # 02_rnaseq/02_crea_liste_rnaseq.py
NAS_DIR   <- env("NASCENT_DIR",   "results/nascent/DESeq2_nascent") # 03_nascent_rna/02_deseq2_nascent.R
OUT_DIR   <- "results/integration/chip_integrazione"
FLANK     <- 1000                       # bp aggiunte a ciascun lato del gene

# fdr = soglia FDR delle regioni DiffBind (1 = tutti i siti del consensus); fold = |Fold| minimo (opzionale)
SETUP <- list(
  Slow_vs_Amr = list(diffbind = file.path(DB_DIR, "Slow_vs_Amr_DiffBind_all_sites.tsv"), sep = "\t", fdr = 0.05,
                     rna = "Slow_vs_Amr", nas = "Slow_vs_Amr", chip_label = "ChIP Slow vs Amr"),
  # NaCl vs UT: ChIP del dataset pubblico KCl vs UT* (da integrare)
  NaCl_vs_UT  = list(diffbind = file.path(DB_DIR, "KCl_vs_UT_DiffBind_all_sites.tsv"), sep = "\t", fdr = 0.01, fold = 1.5,
                     rna = "NACL_vs_WT", nas = "NaCl_vs_UT", chip_label = "ChIP KCl vs UT"),
  Fast_vs_Amr = list(diffbind = file.path(DB_DIR, "Fast_vs_Amr_DiffBind_all_sites.tsv"), sep = "\t", fdr = 1,
                     rna = "Fast_vs_Amr", nas = "Fast_vs_Amr", chip_label = "ChIP Fast vs Amr (tutti i siti)"),
  Storm_vs_UT = list(diffbind = file.path(DB_DIR, "Storm_vs_UT_DiffBind_all_sites.tsv"), sep = "\t", fdr = 1,
                     rna = "STORM_vs_WT", nas = "Storm_vs_UT", chip_label = "ChIP Storm vs UT (tutti i siti)"))

# Confronti da eseguire e cartella di output (da riga di comando)
#   Rscript 02_chip_integrazione.R Fast_vs_Amr chip_integrazione_Fast
args <- commandArgs(trailingOnly = TRUE)
RUN  <- if (length(args) >= 1) strsplit(args[1], ",")[[1]] else c("Slow_vs_Amr", "Fast_vs_Amr", "Storm_vs_UT")
if (length(args) >= 2) OUT_DIR <- args[2]
SETUP <- SETUP[RUN]
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
COL <- c(RNA = "#2a78d6", NAS = "#eb6834", CHIP = "#1baf7a")

# ---- PASSO 1: universo dei geni (BED) ---------------------------------------
bed <- read.delim(BED, header = FALSE, col.names = c("chr", "start", "end", "symbol"))
genes <- GRanges(bed$chr, IRanges(bed$start + 1, bed$end), symbol = bed$symbol)  # BED 0-based
genes_ext <- resize(genes, width(genes) + 2 * FLANK, fix = "center")
universe <- unique(bed$symbol)

# ---- PASSO 2: mappa simbolo <-> ENSG (dai file RNAkinet, stessa annotazione) --
map <- do.call(rbind, lapply(list.files(RNAKINET, "rnakinet.*csv$", full.names = TRUE),
                             function(f) read.csv(f)[, c("gene_id", "gene_name")]))
map$gene_id <- sub("\\.[0-9]+$", "", map$gene_id)
map <- map[!duplicated(map$gene_id), ]
sym_of  <- setNames(map$gene_name, map$gene_id)
ensg_of <- tapply(map$gene_id, map$gene_name, function(x) x, simplify = FALSE)  # simbolo -> ENSG
ensg_from_symbols <- function(sy) unique(unlist(ensg_of[intersect(sy, names(ensg_of))]))

# ---- funzioni ----------------------------------------------------------------
chip_genes <- function(cfg) {
  d  <- read.table(cfg$diffbind, header = TRUE, sep = cfg$sep, stringsAsFactors = FALSE)
  fold_min <- if (is.null(cfg$fold)) 0 else cfg$fold      # |Fold| minimo (solo se indicato)
  d  <- d[d$FDR <= cfg$fdr & abs(d$Fold) >= fold_min, ]
  pk <- GRanges(d$seqnames, IRanges(d$start, d$end), fold = d$Fold)
  ov <- findOverlaps(genes_ext, pk)
  hit <- data.frame(symbol = genes$symbol[queryHits(ov)], fold = pk$fold[subjectHits(ov)])
  agg <- do.call(rbind, lapply(split(hit$fold, hit$symbol), function(f)
    data.frame(n_regioni = length(f), n_gain = sum(f > 0), n_loss = sum(f < 0),
               fold_medio = mean(f))))
  agg$symbol <- rownames(agg)
  agg$chip <- ifelse(agg$n_gain > 0 & agg$n_loss > 0, "misto",
               ifelse(agg$n_gain > 0, "gain", "loss"))
  list(table = agg, n_regions = nrow(d), n_regions_in_genes = length(unique(subjectHits(ov))))
}
# Liste RNA-seq e nascente COMPLETE (le stesse dei Venn RNA-seq x nascente), in ENSG
expr_sets <- function(cfg) {
  up <- read.csv(file.path(RNA_DIR, paste0(cfg$rna, "_UP.csv")), row.names = 1)
  dn <- read.csv(file.path(RNA_DIR, paste0(cfg$rna, "_DOWN.csv")), row.names = 1)
  na <- read.csv(file.path(NAS_DIR, paste0(cfg$nas, "_all_genes.csv")))
  nv <- function(x) unique(sub("\\.[0-9]+$", "", x))
  list(rna_up = nv(rownames(up)), rna_down = nv(rownames(dn)),
       nas_up = nv(na$gene_id[na$status == "UP"]), nas_down = nv(na$gene_id[na$status == "DOWN"]))
}

# Venn a 3 insiemi (cerchi, conteggi per ognuna delle 7 regioni)
circ <- function(x0, y0, r = 1, n = 200) {
  t <- seq(0, 2 * pi, length.out = n); data.frame(x = x0 + r * cos(t), y = y0 + r * sin(t))
}
venn3 <- function(A, B, C, labs, title) {
  r <- function(a, b, c) sum(a & b & c)
  all <- unique(c(A, B, C)); a <- all %in% A; b <- all %in% B; c <- all %in% C
  n <- c(A = r(a, !b, !c), B = r(!a, b, !c), C = r(!a, !b, c), AB = r(a, b, !c),
         AC = r(a, !b, c), BC = r(!a, b, c), ABC = r(a, b, c))
  pos <- data.frame(k = names(n), n = n,
                    x = c(-0.9, 0.9, 0, 0, -0.5, 0.5, 0),
                    y = c(0.45, 0.45, -1.0, 0.65, -0.35, -0.35, 0.05))
  ggplot() +
    geom_polygon(data = circ(-0.55, 0.3), aes(x, y), fill = COL["RNA"], alpha = 0.2, colour = COL["RNA"], linewidth = 0.8) +
    geom_polygon(data = circ(0.55, 0.3),  aes(x, y), fill = COL["NAS"], alpha = 0.2, colour = COL["NAS"], linewidth = 0.8) +
    geom_polygon(data = circ(0, -0.6),    aes(x, y), fill = COL["CHIP"], alpha = 0.2, colour = COL["CHIP"], linewidth = 0.8) +
    geom_text(data = pos, aes(x, y, label = n), size = 3.8,
              fontface = ifelse(pos$k == "ABC", "bold", "plain")) +
    annotate("text", x = -1.25, y = 1.5, label = paste0(labs[1], "\n(", length(A), ")"), size = 3.2, fontface = "bold", colour = "#0b0b0b") +
    annotate("text", x = 1.25,  y = 1.5, label = paste0(labs[2], "\n(", length(B), ")"), size = 3.2, fontface = "bold", colour = "#0b0b0b") +
    annotate("text", x = 0,     y = -1.85, label = paste0(labs[3], " (", length(C), ")"), size = 3.2, fontface = "bold", colour = "#0b0b0b") +
    coord_equal(xlim = c(-1.9, 1.9), ylim = c(-2.0, 1.85)) +
    labs(title = title) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
}

# ---- PASSO 3: per ogni confronto ---------------------------------------------
plots <- list(); summ <- list()
for (nm in names(SETUP)) {
  cfg <- SETUP[[nm]]
  ch <- chip_genes(cfg); s <- expr_sets(cfg)
  ch$table$ensg <- sapply(ch$table$symbol, function(x) paste(ensg_from_symbols(x), collapse = ";"))
  chip_all <- ensg_from_symbols(ch$table$symbol)          # geni ChIP in ENSG

  # Venn unico: significativi UP + DOWN insieme
  rna_sig <- union(s$rna_up, s$rna_down); nas_sig <- union(s$nas_up, s$nas_down)
  plots[[nm]] <- venn3(rna_sig, nas_sig, chip_all,
    c("RNA-seq UP+DOWN", "Nascente UP+DOWN", cfg$chip_label), gsub("_", " ", nm))
  writeLines(paste(Reduce(intersect, list(rna_sig, nas_sig, chip_all)),
                   sym_of[Reduce(intersect, list(rna_sig, nas_sig, chip_all))], sep = "\t"),
             file.path(OUT_DIR, paste0("tripla_", nm, "_UP+DOWN.txt")))

  # tabella per gene: tutti i geni presenti in almeno uno dei tre insiemi
  ids <- unique(c(unlist(s), chip_all))
  t <- data.frame(gene_id = ids, symbol = unname(sym_of[ids]))
  t$RNAseq   <- ifelse(ids %in% s$rna_up, "UP", ifelse(ids %in% s$rna_down, "DOWN", ""))
  t$nascente <- ifelse(ids %in% s$nas_up, "UP", ifelse(ids %in% s$nas_down, "DOWN", ""))
  t$in_BED   <- ids %in% ensg_from_symbols(universe)
  ct <- ch$table[, c("symbol", "n_regioni", "n_gain", "n_loss", "fold_medio", "chip")]
  names(ct)[1] <- "symbol_BED"
  t$symbol_BED <- ifelse(t$in_BED, t$symbol, NA)
  t <- merge(t, ct, by = "symbol_BED", all.x = TRUE)
  t$chip[is.na(t$chip)] <- ""
  t <- t[order(t$chip == "", t$RNAseq == "", t$nascente == ""), c("gene_id", "symbol", "RNAseq",
         "nascente", "in_BED", "chip", "n_regioni", "n_gain", "n_loss", "fold_medio")]
  write.csv(t, file.path(OUT_DIR, paste0("geni_", nm, "_ChIP_RNAseq_nascente.csv")), row.names = FALSE)

  # liste delle intersezioni triple
  for (dir in c("UP", "DOWN")) {
    r <- s[[paste0("rna_", tolower(dir))]]; n <- s[[paste0("nas_", tolower(dir))]]
    tri <- Reduce(intersect, list(r, n, chip_all))
    writeLines(paste(tri, sym_of[tri], sep = "\t"),
               file.path(OUT_DIR, paste0("tripla_", nm, "_", dir, ".txt")))
  }
  summ[[nm]] <- data.frame(
    confronto = nm, chip = cfg$chip_label, regioni_FDR = ch$n_regions,
    regioni_su_geni_BED = ch$n_regions_in_genes, geni_ChIP = length(chip_all),
    geni_gain = sum(ch$table$chip == "gain"), geni_loss = sum(ch$table$chip == "loss"),
    geni_misti = sum(ch$table$chip == "misto"),
    geni_ChIP_ENSG = length(chip_all),
    RNAseq_UP = length(s$rna_up), RNAseq_DOWN = length(s$rna_down),
    nasc_UP = length(s$nas_up), nasc_DOWN = length(s$nas_down),
    tripla_UP = length(Reduce(intersect, list(s$rna_up, s$nas_up, chip_all))),
    tripla_DOWN = length(Reduce(intersect, list(s$rna_down, s$nas_down, chip_all))),
    tripla_UP_DOWN = length(Reduce(intersect, list(rna_sig, nas_sig, chip_all))))
}
sm <- do.call(rbind, summ)
write.csv(sm, file.path(OUT_DIR, "riepilogo_ChIP.csv"), row.names = FALSE)
print(sm, row.names = FALSE)

fig <- wrap_plots(plots, nrow = 1) +
  plot_annotation(title = "RNA-seq x RNA nascente x ChIP (DiffBind): geni significativi UP + DOWN",
                  subtitle = sprintf("RNA-seq e nascente: stesse liste dei Venn RNA-seq x nascente.\nChIP: geni del BED con regioni DiffBind (%s) sul corpo genico +/- %d kb",
                                     paste(sapply(SETUP, function(x) if (x$fdr >= 1) "tutti i siti" else paste0("FDR < ", x$fdr, if (!is.null(x$fold)) paste0(", |Fold| > ", x$fold) else "")), collapse = " / "), FLANK / 1000))
W <- max(9, 1.5 + 4.75 * length(SETUP))
ggsave(file.path(OUT_DIR, "venn_ChIP_RNAseq_nascente.pdf"), fig, width = W, height = 5.6)
ggsave(file.path(OUT_DIR, "venn_ChIP_RNAseq_nascente.png"), fig, width = W, height = 5.6, dpi = 200, bg = "white")

# ---- PASSO 4: arricchimento (Fisher) ----------------------------------------
# Solo i geni del BED possono essere "ChIP": il test si fa sui geni del BED.
# I geni UP/DOWN sono occupati in modo differenziale più (o meno) spesso degli altri?
bed_ensg <- ensg_from_symbols(universe)
enr <- list()
for (nm in names(SETUP)) {
  t <- read.csv(file.path(OUT_DIR, paste0("geni_", nm, "_ChIP_RNAseq_nascente.csv")))
  chip <- t$gene_id[t$chip != ""]
  for (src in c("RNAseq", "nascente")) for (d in c("UP", "DOWN")) {
    g <- intersect(t$gene_id[t[[src]] == d], bed_ensg); oth <- setdiff(bed_ensg, g)
    tab <- matrix(c(sum(g %in% chip), sum(!g %in% chip), sum(oth %in% chip), sum(!oth %in% chip)), 2)
    ft <- fisher.test(tab)
    enr[[length(enr) + 1]] <- data.frame(confronto = nm, dati = src, classe = d,
      n_geni_totali = sum(t[[src]] == d), n_geni_nel_BED = length(g),
      pct_con_ChIP = round(100 * mean(g %in% chip), 1), pct_altri_geni_BED = round(100 * mean(oth %in% chip), 1),
      odds_ratio = round(unname(ft$estimate), 2), p_value = signif(ft$p.value, 3))
  }
}
enr <- do.call(rbind, enr)
write.csv(enr, file.path(OUT_DIR, "arricchimento_ChIP_Fisher.csv"), row.names = FALSE)
print(enr, row.names = FALSE)
