# R / Bioconductor packages used in the analysis (R 4.3, Bioconductor 3.18)
# Rscript environment/install_R_packages.R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran <- c("ggplot2", "ggrepel", "patchwork", "pheatmap", "readxl")
bioc <- c("DESeq2",                       # 1.42.0
          "DiffBind",                     # 3.12.0
          "ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene",
          "clusterProfiler", "enrichplot", "org.Hs.eg.db", "GO.db", "AnnotationDbi",
          "GenomicRanges",
          "apeglm", "ashr")               # optional: LFC shrinkage (the analysis used type = "normal")

install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE)
