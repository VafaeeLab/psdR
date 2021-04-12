install.packages("devtools")
install.packages("roxygen2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleCellExperiment")
BiocManager::install("Linnorm")

install.packages("caTools")
install.packages("Seurat")

