## PCA Plots ##
### Cole Nawrocki ###

rm(list = ls())

# Check that we are in the correct environment
.libPaths()
# [1] "/Users/ccn22/micromamba/envs/geomx-env/lib/R/library"

# Packages
library(dplyr)
library(stringr)
library(ggplot2)
library(ggbiplot)
library(magrittr)
library(patchwork)
library(parallel)
library(edgeR)
options(mc.cores = 10)

## Data
meta <- readxl::read_excel("meta_cleaned_v5.xlsx") |> as.data.frame()
rownames(meta) <- meta$aoi_id
norm <- read.csv(file = "q3norm.csv", row.names = 1)
cts <- read.csv(file = "counts.csv", row.names = 1)
meta$AOI_code %<>% as.factor()
meta$Patient_number %<>% as.factor()

lmem <- function(gene, idx) {
  d <- dplyr::mutate(meta[idx,], expr=(norm[gene, idx] |> unlist()))
  testout <- lme4::lmer(formula = expr ~ 1 + (1 + 1|Patient_number), data = d, REML = F)
  return(deviance(testout))
}

lmod <- function(gene, idx) {
  d <- dplyr::mutate(meta[idx,], expr=(norm[gene, idx] |> unlist()))
  testout <- lm(formula = expr ~ 1 + Patient_number, data = d)
  return(deviance(testout))
}

### IMMUNE ---------------------------------------------------------------------
deviances <- parallel::mclapply(rownames(norm), lmem, idx = (meta$AOI_code == "immune"), mc.cores = 10)
names(deviances) <- rownames(norm)
hvgs <- sort(deviances |> unlist(), decreasing = T)
hvgs <- hvgs[hvgs >= quantile(hvgs, 0.9)] |> names()
mat <- norm[hvgs, (meta$AOI_code == "immune")] |> apply(MARGIN = 1, FUN = scale)
rownames(mat) <- colnames(norm[hvgs, (meta$AOI_code == "immune")])
pca <- prcomp(mat, center = F, scale. = F)
p1 <- ggbiplot::ggbiplot(pca, var.axes = F, groups = meta[(meta$AOI_code == "immune"),]$infiltration_type) + ggtitle("infiltration_type")
p2 <- ggbiplot::ggbiplot(pca, var.axes = F, groups = meta[(meta$AOI_code == "immune"),]$Patient_number) + ggtitle("Patient_number")
p3 <- ggbiplot::ggbiplot(pca, var.axes = F, groups = meta[(meta$AOI_code == "immune"),]$ROI_number) + ggtitle("ROI_number")
#pdf(file = "immune_pca.pdf", width = 12, height = 4)
p1|p2|p3
#dev.off()

### FIBROBLAST -----------------------------------------------------------------
deviances <- parallel::mclapply(rownames(norm), lmem, idx = (meta$AOI_code == "fibroblast"), mc.cores = 10)
names(deviances) <- rownames(norm)
hvgs <- sort(deviances |> unlist(), decreasing = T)
hvgs <- hvgs[hvgs >= quantile(hvgs, 0.9)] |> names()
mat <- norm[hvgs, (meta$AOI_code == "fibroblast")] |> apply(MARGIN = 1, FUN = scale)
rownames(mat) <- colnames(norm[hvgs, (meta$AOI_code == "fibroblast")])
pca <- prcomp(mat, center = F, scale. = F)
p1 <- ggbiplot::ggbiplot(pca, var.axes = F, groups = meta[(meta$AOI_code == "fibroblast"),]$infiltration_type) + ggtitle("infiltration_type")
p2 <- ggbiplot::ggbiplot(pca, var.axes = F, groups = meta[(meta$AOI_code == "fibroblast"),]$Patient_number) + ggtitle("Patient_number")
p3 <- ggbiplot::ggbiplot(pca, var.axes = F, groups = meta[(meta$AOI_code == "fibroblast"),]$ROI_number) + ggtitle("ROI_number")
#pdf(file = "fibroblast_pca.pdf", width = 12, height = 4)
p1|p2|p3
#dev.off()

### TUMOR ----------------------------------------------------------------------
deviances <- parallel::mclapply(rownames(norm), lmem, idx = (meta$AOI_code == "tumor"), mc.cores = 10)
names(deviances) <- rownames(norm)
hvgs <- sort(deviances |> unlist(), decreasing = T)
hvgs <- hvgs[hvgs >= quantile(hvgs, 0.9)] |> names()
mat <- norm[hvgs, (meta$AOI_code == "tumor")] |> apply(MARGIN = 1, FUN = scale)
rownames(mat) <- colnames(norm[hvgs, (meta$AOI_code == "tumor")])
pca <- prcomp(mat, center = F, scale. = F)
p1 <- ggbiplot::ggbiplot(pca, var.axes = F, groups = meta[(meta$AOI_code == "tumor"),]$infiltration_type) + ggtitle("infiltration_type")
p2 <- ggbiplot::ggbiplot(pca, var.axes = F, groups = meta[(meta$AOI_code == "tumor"),]$Patient_number) + ggtitle("Patient_number")
p3 <- ggbiplot::ggbiplot(pca, var.axes = F, groups = meta[(meta$AOI_code == "tumor"),]$ROI_number) + ggtitle("ROI_number")
#pdf(file = "tumor_pca.pdf", width = 12, height = 4)
p1|p2|p3
#dev.off()

### NOTES ----------------------------------------------------------------------
# I sent an email to the people at NanoString about a package that improves normalization. 
# Hopefully, I can use that to improve the results shown above. 

