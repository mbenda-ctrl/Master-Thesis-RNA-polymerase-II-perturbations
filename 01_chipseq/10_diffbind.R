#!/usr/bin/env Rscript
# =============================================================================
#  PASSO 10 - Occupazione differenziale di Pol II con DiffBind (3.12.0)
#
#  Per ogni confronto (config/diffbind_contrasts.tsv):
#    - oggetto DBA con i 2 replicati del trattato e i 2 del controllo
#      (BAM del passo 3 + picchi narrow del passo 7; per il dataset pubblico i .xls di MACS)
#    - dba.count con le impostazioni di default: picchi ricentrati sul summit e portati a
#      401 bp (+/- 200 bp), consensus = picchi presenti in almeno 2 campioni
#    - dba.analyze con DESeq2 (normalizzazione di default sulla dimensione della libreria)
#    - report di TUTTI i siti del consensus (th = 1), con conte normalizzate per campione
#  Colonne utili del report: Conc (log2 concentrazione media), Fold (log2FC, > 0 = piu' segnale
#  nel trattato), p.value, FDR (Benjamini-Hochberg).
#  Soglie usate a valle: FDR < 0.05 (Slow vs Amr); FDR < 0.01 e |Fold| > 1.5 (KCl vs UT).
#
#  Uso: Rscript 10_diffbind.R [samples.tsv] [contrasts.tsv] [cartella_BAM] [cartella_picchi] [output]
# =============================================================================
suppressPackageStartupMessages(library(DiffBind))

args      <- commandArgs(trailingOnly = TRUE)
SAMPLES   <- if (length(args) >= 1) args[1] else "config/diffbind_samples.tsv"
CONTRASTS <- if (length(args) >= 2) args[2] else "config/diffbind_contrasts.tsv"
BAM_DIR   <- if (length(args) >= 3) args[3] else "results/chipseq/bam"
PEAK_DIR  <- if (length(args) >= 4) args[4] else "results/chipseq/peaks/narrow"
OUT_DIR   <- if (length(args) >= 5) args[5] else "results/chipseq/diffbind"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ss <- read.delim(SAMPLES, stringsAsFactors = FALSE)
cc <- read.delim(CONTRASTS, stringsAsFactors = FALSE)

for (i in seq_len(nrow(cc))) {
  nm <- cc$contrast[i]; trt <- cc$treated[i]; ctl <- cc$control[i]
  s  <- ss[ss$Condition %in% c(trt, ctl), ]
  sheet <- data.frame(SampleID   = s$SampleID,
                      Condition  = s$Condition,
                      Replicate  = s$Replicate,
                      bamReads   = file.path(BAM_DIR, s$bam),
                      Peaks      = file.path(PEAK_DIR, s$peaks),
                      PeakCaller = cc$peak_format[i],          # "narrow" (narrowPeak) o "macs" (.xls)
                      stringsAsFactors = FALSE)
  miss <- c(sheet$bamReads, sheet$Peaks)[!file.exists(c(sheet$bamReads, sheet$Peaks))]
  if (length(miss)) { message(nm, ": file mancanti, salto\n  ", paste(miss, collapse = "\n  ")); next }
  message("== ", nm, " (", nrow(sheet), " campioni)")

  dba_obj <- dba(sampleSheet = sheet)
  dba_cnt <- dba.count(dba_obj, bParallel = FALSE)        # bParallel = FALSE: meno memoria
  dba_con <- dba.contrast(dba_cnt, contrast = c("Condition", trt, ctl))
  dba_res <- dba.analyze(dba_con, method = DBA_DESEQ2)

  rep_all <- as.data.frame(dba.report(dba_res, method = DBA_DESEQ2, th = 1,
                                      bCounts = TRUE, bCalled = TRUE))
  write.table(rep_all, file.path(OUT_DIR, paste0(nm, "_DiffBind_all_sites.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  saveRDS(dba_res, file.path(OUT_DIR, paste0(nm, "_dba.rds")))
  message(sprintf("   siti testati %d | FDR < 0.05: %d", nrow(rep_all), sum(rep_all$FDR < 0.05)))

  pdf(file.path(OUT_DIR, paste0(nm, "_QC.pdf")), width = 7, height = 6)
  dba.plotPCA(dba_cnt, label = DBA_CONDITION)
  dba.plotMA(dba_res)
  invisible(dev.off())
}
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
