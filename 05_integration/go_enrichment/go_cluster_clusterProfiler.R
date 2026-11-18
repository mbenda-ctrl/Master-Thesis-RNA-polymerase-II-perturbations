#!/usr/bin/env Rscript
# =============================================================================
#  Arricchimento GO dei cluster di una heatmap con clusterProfiler
#
#  - gruppi: colonna "gruppo" della tabella della heatmap (gene symbol come nome di riga)
#  - clusterProfiler::compareCluster(fun = "enrichGO"): test ipergeometrico, BH,
#    termini con 10-500 geni, annotazioni org.Hs.eg.db
#  - riduzione ridondanza: clusterProfiler::simplify (similarità semantica di Wang >= 0.7,
#    tiene il termine con p.adjust più basso)
#  - universo:
#      genoma = default di enrichGO (tutti i geni con annotazione GO nell'OrgDb)
#      bed    = geni del BED (conteggi_gene.tab) testati in RNA-seq e nel nascent
#
#  Uso: Rscript go_cluster_clusterProfiler.R <tabella_heatmap.tsv> <cartella_output> [BP|MF|CC] [genoma|bed] [n_termini_plot]
#  Esempio: Rscript go_cluster_clusterProfiler.R heatmap/heatmap_rna/heatmap_rna_Slow_vs_Amr.tsv GO_CP_Slow_genoma BP genoma 25
# =============================================================================
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot); library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
HEAT <- args[1]; OUT <- args[2]
ONT  <- if (length(args) >= 3) args[3] else "BP"
UNIV <- if (length(args) >= 4) args[4] else "genoma"
NSHOW <- if (length(args) >= 5) as.integer(args[5]) else 10
CONFR <- sub("heatmap_(rna|chip)_(.*)\\.tsv", "\\2", basename(HEAT))
# file usati solo per l'universo "bed"
BED <- Sys.getenv("GENE_COUNTS", "results/chipseq/gene_counts/conteggi_gene.tab")   # 01_chipseq/12b
RNA_ALL <- file.path(Sys.getenv("RNA_ALL_DIR", "results/rnaseq/RNAseq"), sprintf("%s_all_genes.csv", CONFR))
NAS_ALL <- file.path(Sys.getenv("NASCENT_DIR", "results/nascent/DESeq2_nascent"), sprintf("%s_all_genes.csv", CONFR))
MAP_DIR <- Sys.getenv("RNAKINET_DIR", "results/nascent/rnakinet_genes")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

H <- read.delim(HEAT, row.names = 1)
cl <- split(rownames(H), paste0("Cluster ", H$gruppo))
cl <- cl[order(as.integer(sub("Cluster ", "", names(cl))))]
names(cl) <- sprintf("%s\n(n = %d)", names(cl), lengths(cl))

universe <- NULL
if (UNIV == "bed") {
  nv  <- function(x) sub("\\.[0-9]+$", "", x)
  map <- do.call(rbind, lapply(list.files(MAP_DIR, "rnakinet.*csv$", full.names = TRUE),
                               function(f) read.csv(f)[, c("gene_id", "gene_name")]))
  map$gene_id <- nv(map$gene_id)
  bed <- rownames(read.delim(BED, row.names = 1))
  rna_ids <- nv(read.csv(RNA_ALL)$gene_id)
  nas <- read.csv(NAS_ALL); nas_ids <- nv(nas$gene_id[!is.na(nas$log2FC_MLE)])
  tested <- unique(map$gene_name[map$gene_id %in% intersect(rna_ids, nas_ids)])
  universe <- union(intersect(bed, tested), rownames(H))
  cat(sprintf("Universo BED: %d geni\n", length(universe)))
} else cat("Universo: genoma (default enrichGO)\n")

cc <- compareCluster(geneClusters = cl, fun = "enrichGO", OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                     ont = ONT, universe = universe, pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500, readable = FALSE)
if (is.null(cc) || nrow(as.data.frame(cc)) == 0) stop("Nessun termine significativo in nessun cluster")
ccs <- simplify(cc, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")

tab  <- as.data.frame(cc);  tab$Cluster  <- sub("\n.*", "", tab$Cluster)
tabs <- as.data.frame(ccs); tabs$Cluster <- sub("\n.*", "", tabs$Cluster)
write.csv(tab,  file.path(OUT, sprintf("GO_%s_clusterProfiler_completo.csv", ONT)), row.names = FALSE)
write.csv(tabs, file.path(OUT, sprintf("GO_%s_clusterProfiler_simplify.csv", ONT)), row.names = FALSE)
print(table(tab$Cluster)); print(table(tabs$Cluster))

p <- dotplot(ccs, showCategory = NSHOW, font.size = 9, label_format = 60) +
  ggtitle(sprintf("GO %s - %s (clusterProfiler)", ONT, CONFR),
          subtitle = sprintf("universo: %s; termini 10-500 geni, BH padj < 0.05; simplify (Wang >= 0.7); primi %d per cluster",
                             if (UNIV == "bed") "geni del BED testati in RNA-seq e nascent" else "genoma (org.Hs.eg.db)", NSHOW)) +
  theme(plot.subtitle = element_text(size = 8))
h <- max(5, 0.24 * length(unique(p$data$Description)) + 2.5)
ggsave(file.path(OUT, sprintf("GO_%s_clusterProfiler_dotplot.pdf", ONT)), p, width = 10, height = h, limitsize = FALSE)
ggsave(file.path(OUT, sprintf("GO_%s_clusterProfiler_dotplot.png", ONT)), p, width = 10, height = h, dpi = 250, bg = "white", limitsize = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
