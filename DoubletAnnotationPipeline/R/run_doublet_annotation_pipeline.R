library(Seurat)
library(Signac)
library(dplyr)
library(tidyr)
library(textclean)
library(ggplot2)
library(MASS)
require(stringr)
library(mclust)

source("R/getMarkerPeak_method.R")
source("R/annotateDoublets_method.R")

pbmc.atac <- readRDS(file = "data/example_seurat_pbmc1.Rds")

Idents(pbmc.atac) <- "grouped_clusters"

doublets <- read.table("../NewDoubletIds/PBMC1_doublets.txt") %>% t() %>% as.vector()

marker_peaks <- getMarkerPeaks(pbmc.atac, doublets = doublets)

doublet_annotations <- annotateDoublets(obj = pbmc.atac, marker_peaks = marker_peaks, doublets = doublets)
