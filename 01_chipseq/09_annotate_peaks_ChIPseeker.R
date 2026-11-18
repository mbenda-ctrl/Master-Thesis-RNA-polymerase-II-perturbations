#!/usr/bin/env Rscript
# =============================================================================
#  PASSO 9 - Annotazione genomica dei picchi narrow e broad con ChIPseeker
#
#  annotatePeak con regione promotore = TSS +/- 1 kb (tssRegion = c(-1000, 1000))
#  e annotazione TxDb.Hsapiens.UCSC.hg38.knownGene (simboli da org.Hs.eg.db).
#  Per ogni tipo di picco (narrow, broad):
#    - tabella annotata per campione  (<campione>_annotated.tsv)
#    - riepilogo delle categorie       (annotation_summary_<tipo>.tsv)
#    - barplot della distribuzione     (annotation_bar_<tipo>.pdf)
#  Separazione promoter / gene body: i picchi narrow annotati "Promoter" descrivono
#  la regione del promotore, i picchi broad non-Promoter il corpo genico.
#
#  Uso: Rscript 09_annotate_peaks_ChIPseeker.R [cartella_picchi] [cartella_output]
#       (cartella_picchi contiene le sottocartelle narrow/ e broad/ del passo 7)
# =============================================================================
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(ggplot2)
})

args     <- commandArgs(trailingOnly = TRUE)
PEAK_DIR <- if (length(args) >= 1) args[1] else "results/chipseq/peaks"
OUT_DIR  <- if (length(args) >= 2) args[2] else "results/chipseq/annotation"
TSS      <- c(-1000, 1000)      # regione promotore: TSS +/- 1 kb
txdb     <- TxDb.Hsapiens.UCSC.hg38.knownGene

for (tipo in c("narrow", "broad")) {
  ext   <- if (tipo == "narrow") "_peaks.narrowPeak$" else "_peaks.broadPeak$"
  files <- list.files(file.path(PEAK_DIR, tipo), pattern = ext, full.names = TRUE)
  if (!length(files)) { message("nessun picco ", tipo, " in ", PEAK_DIR); next }
  names(files) <- sub(paste0("_", tipo, ext), "", basename(files))
  od <- file.path(OUT_DIR, tipo); dir.create(od, recursive = TRUE, showWarnings = FALSE)

  ann <- lapply(files, function(f)
    annotatePeak(readPeakFile(f), tssRegion = TSS, TxDb = txdb,
                 annoDb = "org.Hs.eg.db", verbose = FALSE))

  for (s in names(ann))
    write.table(as.data.frame(ann[[s]]), file.path(od, paste0(s, "_annotated.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)

  # riepilogo: % di picchi per categoria (Promoter, 5'UTR, 3'UTR, esoni, introni, intergenico...)
  summ <- do.call(rbind, lapply(names(ann), function(s) {
    st <- ann[[s]]@annoStat; data.frame(sample = s, Feature = st$Feature, Frequency = st$Frequency)
  }))
  write.table(summ, file.path(OUT_DIR, paste0("annotation_summary_", tipo, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  pdf(file.path(OUT_DIR, paste0("annotation_bar_", tipo, ".pdf")), width = 9, height = 0.45 * length(ann) + 2)
  print(plotAnnoBar(ann) + ggtitle(paste("Annotazione picchi", tipo)))
  print(plotDistToTSS(ann) + ggtitle(paste("Distanza dal TSS -", tipo)))
  invisible(dev.off())
  message(tipo, ": annotati ", length(ann), " campioni")
}
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
