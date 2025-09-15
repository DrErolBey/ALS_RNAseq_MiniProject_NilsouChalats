# 00_setup.R ??? package install
# Author: Nilsou Chalats
# R: 4.3.x | Bioc: 3.18

# code/00_setup.R  (install to user library)
options(repos=c(CRAN="https://cloud.r-project.org"), timeout=1200)

# Windows user lib
userlib <- file.path(Sys.getenv("USERPROFILE"), "R", "win-library",
                     paste0(R.version$major,".",R.version$minor))
dir.create(userlib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(userlib, .libPaths()))

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")

# Bioconductor
BiocManager::install(c(
  "DESeq2","apeglm","AnnotationDbi","org.Hs.eg.db","clusterProfiler","GO.db","GOSemSim","DOSE"
), ask=FALSE, update=FALSE)

# CRAN
install.packages(c("readr","dplyr","ggplot2"), dependencies=TRUE)

# smoke test
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
  library(DESeq2); library(apeglm)
  library(AnnotationDbi); library(org.Hs.eg.db)
  library(clusterProfiler)
})
cat("Setup OK\n")

