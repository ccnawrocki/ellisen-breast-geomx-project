### J9/J10 Presentation Analysis ###

rm(list = ls())
.rs.restartR(clean = T)
.libPaths()

cts <- read.csv("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/analysis/ellisen-breast/counts.csv", row.names = 1) |> as.matrix()
meta <- openxlsx::read.xlsx("~/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/analysis/ellisen-breast/meta_cleaned_v5.xlsx")

all(meta$aoi_id == colnames(cts)) # TRUE

library(edgeR)
library(limma)
library(dplyr)
library(magrittr)
library(dplyr)

# Background subtraction... seems to help sometimes, but not doing here.
# bgsub <- sweep(x = cts, MARGIN = 2, STATS = cts["NegProbe-WTX",], FUN = "-") |> apply(MARGIN = 2, FUN = pmax, 0)
# keep <- rowMeans(bgsub) >= 0 & rowMeans(bgsub) <= 45 & rownames(bgsub) != "NegProbe-WTX"
# cts_filt <- bgsub[keep,]

keep <- rowMeans(cts) >= -Inf & rowMeans(cts) <= 25 & rownames(cts) != "NegProbe-WTX" # 12093 genes kept
cts_filt <- cts[keep,]

# Size factors via TMM which is typical for sequencing data
y <- DGEList(counts = cts_filt) |> edgeR::calcNormFactors(method = "upperquartile") # 14478

# voom and model fitting
# vlmseq <- voomLmFit(counts = y, normalize.method = "none",
#                     block = meta$Patient_number,
#                     adaptive.span = T, plot = T, save.plot = T, keep.EList = T)
# saveRDS(vlmseq, file = "full_vlmseq.RDS")
vlmseq <- readRDS("full_vlmseq.RDS")

# HVG selection
n_hvg <- 1500
resids <- (vlmseq$voom.xy$y[sort(vlmseq$voom.xy$x, decreasing = F) |> names()] - vlmseq$voom.line$y)
hvg_genes <- sort(resids, decreasing = T)[1:n_hvg] |> names()

# After: 
pdf(file = "voom_plot.pdf", width = 6, height = 6)
plot(vlmseq$voom.xy$x, vlmseq$voom.xy$y, 
     pch = 16, cex = 0.5, col = ifelse(test = names(vlmseq$voom.xy$x) %in% hvg_genes, yes = "red", no = "black"),
     xlab = "log2( count size + 0.5 )",
     ylab = "Sqrt( standard deviation )",
     main = "voom: Mean-variance trend")
lines(vlmseq$voom.line, col = "red", lwd = 2)
dev.off()

# Unsupervised clustering
library(ComplexHeatmap)
png(file = "unsup_clustering.png", width = 6, height = 8, res = 200, units = "in")
mat <- vlmseq$EList$E[hvg_genes,] |> apply(MARGIN = 1, FUN = scale) |> t()
Heatmap(matrix = mat, show_row_names = F, name = "scaled\nlogCPM", 
        top_annotation = HeatmapAnnotation(type = meta$AOI_code, col = list("type" = c("tumor" = "green4", "fibroblast" = "dodgerblue", "immune" = "red")), show_annotation_name = F), 
        col = circlize::colorRamp2(colors = viridis::viridis(n = 101, option = "C"), 
                                   breaks = seq(quantile(mat, 0.015), quantile(mat, 0.995), length.out = 101)))
dev.off()

png(file = "unsup_clustering_batch_corrected.png", width = 6, height = 8, res = 200, units = "in")
bc <- limma::removeBatchEffect(x = vlmseq$EList$E, batch = meta$Patient_number)
mat <- bc[hvg_genes,] |> apply(MARGIN = 1, FUN = scale) |> t()
Heatmap(matrix = mat, show_row_names = F, name = "scaled\nlogCPM", 
        top_annotation = HeatmapAnnotation(type = meta$AOI_code, col = list("type" = c("tumor" = "green4", "fibroblast" = "dodgerblue", "immune" = "red")), show_annotation_name = F), 
        col = circlize::colorRamp2(colors = viridis::viridis(n = 101, option = "C"), 
                                   breaks = seq(quantile(mat, 0.015), quantile(mat, 0.995), length.out = 101)))
dev.off()

# Segment markers... using presto for speed
markers <- presto::wilcoxauc(X = vlmseq$EList$E, y = meta$AOI_code)
markers_summary <- filter(markers, padj < 0.05) |> group_by(group) |> slice_max(order_by = logFC, n = 20) |> mutate(group = factor(group, levels = c("immune", "fibroblast", "tumor"))) |> arrange(desc(group))
mat <- bc[markers_summary$feature,] |> apply(MARGIN = 1, FUN = scale) |> t()

png(file = "segment_markers.png", width = 8, height = 8, res = 300, units = "in")
Heatmap(matrix = mat, column_split = meta$AOI_code, name = "scaled\nlogCPM", cluster_rows = F, row_names_gp = grid::gpar(fontsize = 7), 
        top_annotation = HeatmapAnnotation(type = meta$AOI_code, col = list("type" = c("tumor" = "green4", "fibroblast" = "dodgerblue", "immune" = "red")), show_annotation_name = F), 
        col = circlize::colorRamp2(colors = viridis::viridis(n = 101, option = "C"), breaks = seq(quantile(mat, 0.01), quantile(mat, 0.99), length.out = 101)))
dev.off()

# PCA
# pcs <- irlba::prcomp_irlba(x = vlmseq$EList$E[hvg_genes,] |> apply(MARGIN = 1, FUN = scale), n = 10, center = F, scale. = F)
pcs <- irlba::prcomp_irlba(x = bc[hvg_genes,] |> apply(MARGIN = 1, FUN = scale), n = 10, center = F, scale. = F)
plot(pcs)

pdf("pca.pdf", width = 6, height = 6)
plot(pcs$x[,1], pcs$x[,2], 
     pch = 16, cex = 2.5,
     col = case_when(meta$AOI_code == "tumor" ~ "green4", meta$AOI_code == "fibroblast" ~ "dodgerblue", meta$AOI_code == "immune" ~ "red"),
     xlab = paste("PC1", " (", round(pcs$sdev[1]**2/pcs$totalvar,3)*100, "%)", sep = ""), ylab = paste("PC1", " (", round(pcs$sdev[2]**2/pcs$totalvar, 3)*100, "%)", sep = "")
     )
dev.off()

# UMAP
pdf("umap.pdf", width = 6, height = 6)
um <- uwot::umap(X = pcs$x, n_neighbors = 20, spread = 3, metric = "cosine")
plot(um, 
     pch = 16, cex = 1, 
     col = case_when(meta$AOI_code == "tumor" ~ "green4", meta$AOI_code == "fibroblast" ~ "dodgerblue", meta$AOI_code == "immune" ~ "red"), 
     xlab = "UMAP1", ylab = "UMAP2",
     )
dev.off()

# DE analysis
tumoridx <- meta$AOI_code == "tumor"
tumor_cts <- cts[,tumoridx]
tumor_meta <- meta[tumoridx,]
tumor_meta[,c("Patient_number", "ROI_number", "AOI_number", "infiltration_type")] |> 
  arrange(infiltration_type) |> View()

table(tumor_meta$Patient_number, tumor_meta$infiltration_category)
#   desert excluded inflamed
# 1      1       14        0 --> take the excluded
# 2      4        1        4 --> take the inflamed
# 3      5        1        6 --> take the inflamed
# 4      0        0        9 --> take the inflamed
# 5      0        2       10 --> take the inflamed
# 6      1        9        2 --> take the excluded
# 7      0       10        0 --> take the excluded
# 8      0        5        3 --> take the excluded

# We will compare inflamed vs. excluded, subsetting to make the modeling easier
subtumormeta <- tumor_meta[
  (tumor_meta$Patient_number == 1 & tumor_meta$infiltration_category == "excluded") |
  (tumor_meta$Patient_number == 2 & tumor_meta$infiltration_category == "inflamed") |
  (tumor_meta$Patient_number == 3 & tumor_meta$infiltration_category == "inflamed") |
  (tumor_meta$Patient_number == 4 & tumor_meta$infiltration_category == "inflamed") |
  (tumor_meta$Patient_number == 5 & tumor_meta$infiltration_category == "inflamed") |
  (tumor_meta$Patient_number == 6 & tumor_meta$infiltration_category == "excluded") |
  (tumor_meta$Patient_number == 7 & tumor_meta$infiltration_category == "excluded") |
  (tumor_meta$Patient_number == 8 & tumor_meta$infiltration_category == "excluded")
,]
table(subtumormeta$Patient_number, subtumormeta$infiltration_category)
#   excluded inflamed
# 1       14        0
# 2        0        4
# 3        0        6
# 4        0        9
# 5        0       10
# 6        9        0
# 7       10        0
# 8        5        0

# Modeling
subtumorcts <- tumor_cts[,subtumormeta$aoi_id]
keep <- rowMeans(subtumorcts) >= -Inf & rowMeans(subtumorcts) <= 50 & rownames(subtumorcts) != "NegProbe-WTX" # 14107 genes kept
cts_filt <- subtumorcts[keep,]
y <- DGEList(counts = cts_filt) |> calcNormFactors(method = "upperquartile")

subtumormeta$infiltration_category %<>% as.factor()
mm <- model.matrix(~0+infiltration_category, data = subtumormeta)
colnames(mm) <- levels(subtumormeta$infiltration_category)
v <- voomLmFit(counts = y, normalize.method = "none", sample.weights = T, 
               design = mm, 
               block = subtumormeta$Patient_number,
               plot = T, save.plot = , keep.EList = T, adaptive.span = T)

# Testing
contrast.matrix <- makeContrasts(inflamed - excluded, levels = mm)
out <- v |> contrasts.fit(contrasts = contrast.matrix) |> eBayes() |> topTable(number = Inf, coef = 1)

toplot <- out[out$adj.P.Val < 0.05 & abs(out$logFC) > 1,] |> arrange(logFC) |> rownames()
rownames(subtumormeta) <- subtumormeta$aoi_id
plotidx <- arrange(subtumormeta, infiltration_category, Patient_number) |> rownames()
mat <- v$EList$E[toplot, plotidx] |> apply(MARGIN = 1, FUN = scale) |> t()

png(filename = "inflamed_vs_excluded_heatmap.png", width = 6, height = 8, res = 300, units = "in")
Heatmap(matrix = mat, show_row_names = T, name = "scaled\nlogCPM", cluster_rows = F, cluster_columns = T,
        column_split = subtumormeta[plotidx,]$infiltration_category,
        row_names_gp = grid::gpar(fontsize = 6),
        top_annotation = HeatmapAnnotation(group = subtumormeta[plotidx,]$infiltration_category,
                                           patient = subtumormeta[plotidx,]$Patient_number,
                                           show_annotation_name = F,
                                           col = list("group" = c("excluded" = "turquoise", "inflamed" = "violet"),
                                                      "patient" = ggprism::prism_color_pal()(8) |> setNames(1:8))), 
        col = circlize::colorRamp2(colors = viridis::viridis(n = 101, option = "C"), 
                                   breaks = seq(quantile(mat, 0.015), quantile(mat, 0.995), length.out = 101)))
dev.off()

library(ggplot2)
out$target <- rownames(out)
ggplot() + 
  # scattermore::geom_scattermore(data = out[out$adj.P.Val >= 0.05,], mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "grey", pointsize = 2) + 
  geom_point(data = out[out$adj.P.Val >= 0.05,], mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "grey", shape = ".") +
  geom_point(data = out[out$adj.P.Val < 0.05 & out$logFC > 0,], mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "violet", size = 2, shape = 16) + 
  geom_point(data = out[out$adj.P.Val < 0.05 & out$logFC < 0,], mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "turquoise", size = 2, shape = 16) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggrepel::geom_text_repel(data = out[out$adj.P.Val < 0.1,], mapping = aes(x = logFC, y = -log10(adj.P.Val), label = target), min.segment.length = 0, box.padding = 0.25, max.overlaps = 15) +
  ggthemes::theme_par() + 
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))
ggsave(filename = "inflamed_vs_excluded_volcano.pdf", width = 10, height = 10)

mean(out$logFC > 0) # 0.3552965
mean(out$adj.P.Val < 0.05) # 0.0103531

library(msigdbr)
hallmarkgenesets <- msigdbr(db_species = "HS", collection = "H") |> 
  dplyr::select(gs_name, gene_symbol)

prlist <- arrange(out, desc(t)) |> pull(t)
names(prlist) <- arrange(out, desc(t)) |> pull(target)

library(clusterProfiler)
set.seed(2001)
hallmark <- GSEA(geneList = prlist, TERM2GENE = hallmarkgenesets, eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
hallmark@result$ID %<>% factor(levels = (arrange(hallmark@result, desc(NES)) |> pull(ID)))

dat <- hallmark@result |> mutate(updown = (sign(NES) > 0)) |> group_by(updown) |> top_n(n = -10, wt = p.adjust)
ggplot(data = dat) + 
  geom_bar(mapping = aes(y = ID, x = NES, fill = updown), stat = "identity") +
  geom_text(data = dat |> filter(NES < 0),
            mapping = aes(x = 0.1, y = ID, label = ID |> gsub(x = _, pattern = "HALLMARK_", replacement = "")),
            size = 4.5, hjust = 0) +
  geom_text(data = dat |> filter(NES > 0),
            mapping = aes(x = -0.1, y = ID, label = ID |> gsub(x = _, pattern = "HALLMARK_", replacement = "")),
            size = 4.5, hjust = 1) +
  geom_vline(xintercept = 0) + 
  scale_fill_manual(values = c("turquoise", "violet")) +
  ggthemes::theme_par() +
  theme(axis.text.y.left = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), 
        axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20), plot.title = element_text(size = 20)) +
  Seurat::NoLegend() +
  labs(title = "Hallmark Gene Sets")
ggsave(filename = "inflamed_vs_excluded_pathways.pdf", width = 10, height = 10)

gseaplot(x = hallmark, geneSetID = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", color.line = "turquoise") +
  ggthemes::theme_par()
gseaplot(x = hallmark, geneSetID = "HALLMARK_KRAS_SIGNALING_DN", color.line = "turquoise") + 
  ggthemes::theme_par()
ggsave(filename = "inflamed_vs_excluded_HALLMARK_KRAS_SIGNALING_DN_GSEA_plot.pdf", width = 6, height = 7)
gseaplot(x = hallmark, geneSetID = "HALLMARK_MYC_TARGETS_V1", color.line = "violet") + 
  ggthemes::theme_par()
ggsave(filename = "inflamed_vs_excluded_HALLMARK_MYC_TARGETS_V1_GSEA_plot.pdf", width = 6, height = 7)
gseaplot(x = hallmark, geneSetID = "HALLMARK_OXIDATIVE_PHOSPHORYLATION", color.line = "violet") + 
  ggthemes::theme_par()

# # Compare fibroblast of hot tumors vs. fibroblast of cold tumors
# fibidx <- meta$AOI_code == "fibroblast"
# fib_cts <- cts[,fibidx]
# fib_meta <- meta[fibidx,]
# fib_meta[,c("Patient_number", "ROI_number", "AOI_number", "infiltration_type")] |> 
#   arrange(infiltration_type) |> View()
# 
# table(fib_meta$Patient_number, fib_meta$infiltration_category)
# #   desert excluded inflamed
# # 1      0        7        0 --> excluded
# # 2      0        1        3 --> inflamed
# # 3      3        0        0
# # 4      0        0        1 --> inflamed
# # 5      0        0        4 --> inflamed
# # 6      1        2        1
# # 7      0        3        0 --> excluded
# # 8      0        1        1
# 
# subfibmeta <- fib_meta[
#   (fib_meta$Patient_number == 1 & fib_meta$infiltration_category == "excluded") |
#   (fib_meta$Patient_number == 2 & fib_meta$infiltration_category == "inflamed") |
#   (fib_meta$Patient_number == 4 & fib_meta$infiltration_category == "inflamed") |
#   (fib_meta$Patient_number == 5 & fib_meta$infiltration_category == "inflamed") |
#   (fib_meta$Patient_number == 7 & fib_meta$infiltration_category == "excluded")
#   ,]
# table(subfibmeta$Patient_number, subfibmeta$infiltration_category)
# 
# subfibcts <- fib_cts[,subfibmeta$aoi_id]
# keep <- rowMeans(subfibcts) >= 1 & rowMeans(subfibcts) <= 45 & rownames(subfibcts) != "NegProbe-WTX" # 6657 genes kept
# cts_filt <- subfibcts[keep,]
# y <- DGEList(counts = cts_filt) |> calcNormFactors(method = "upperquartile")
# 
# subfibmeta$infiltration_category %<>% as.factor()
# mm <- model.matrix(~0+infiltration_category, data = subfibmeta)
# colnames(mm) <- levels(subfibmeta$infiltration_category)
# v <- voomLmFit(counts = y, normalize.method = "quantile", sample.weights = T, 
#                design = mm, 
#                block = subfibmeta$Patient_number,
#                plot = T, save.plot = , keep.EList = T, adaptive.span = T)
# 
# contrast.matrix <- makeContrasts(inflamed - excluded, levels = mm)
# out <- v |> contrasts.fit(contrasts = contrast.matrix) |> eBayes() |> topTable(number = Inf, coef = 1)
# 
# toplot <- out[out$adj.P.Val < 0.1,] |> arrange(logFC) |> rownames()
# rownames(subfibmeta) <- subfibmeta$aoi_id
# plotidx <- arrange(subfibmeta, infiltration_category, Patient_number) |> rownames()
# mat <- v$EList$E[toplot, plotidx] |> apply(MARGIN = 1, FUN = scale) |> t()
# Heatmap(matrix = mat, show_row_names = T, name = "scaled\nlogCPM", cluster_rows = F, cluster_columns = T,
#         column_split = subfibmeta[plotidx,]$infiltration_category,
#         top_annotation = HeatmapAnnotation(group = subfibmeta[plotidx,]$infiltration_category,
#                                            patient = subfibmeta[plotidx,]$Patient_number,
#                                            col = list("group" = c("excluded" = "darkblue", "inflamed" = "turquoise"), 
#                                                       "patient" = ggprism::prism_color_pal()(5) |> setNames(c(1, 2, 4, 5, 7)))), 
#         col = circlize::colorRamp2(colors = viridis::viridis(n = 101, option = "C"), 
#                                    breaks = seq(quantile(mat, 0.015), quantile(mat, 0.995), length.out = 101)))
# 
# out$target <- rownames(out)
# ggplot() + 
#   scattermore::geom_scattermore(data = out[out$adj.P.Val >= 0.05,], mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "grey", pointsize = 2) + 
#   geom_point(data = out[out$adj.P.Val < 0.05 & out$logFC > 0,], mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "turquoise", size = 2, shape = 16) + 
#   geom_point(data = out[out$adj.P.Val < 0.05 & out$logFC < 0,], mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "darkblue", size = 2, shape = 16) +
#   geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   ggrepel::geom_text_repel(data = out[out$adj.P.Val < 0.05,], mapping = aes(x = logFC, y = -log10(adj.P.Val), label = target), min.segment.length = 0, box.padding = 0.25, max.overlaps = 15) +
#   ggthemes::theme_par()

# Compare tumor border vs. tumor core
tumor_meta$position <- ifelse(test = grepl(pattern = "c", x = tumor_meta$infiltration_type), yes = "core", no = "border")
table(tumor_meta$Patient_number, tumor_meta$position)
#   border core
# 1      9    6
# 2      9    0 --> omit
# 3      7    5
# 4      6    3
# 5      9    3
# 6     10    2
# 7      8    2
# 8      6    2

# We will compare inflamed vs. excluded, subsetting to make the modeling easier
subtumormeta <- tumor_meta[
  tumor_meta$Patient_number != 2
  ,]
subtumorcts <- tumor_cts[,subtumormeta$aoi_id]
keep <- rowMeans(subtumorcts) >= -Inf & rowMeans(subtumorcts) <= 50 & rownames(subtumorcts) != "NegProbe-WTX" # 14592 genes kept
cts_filt <- subtumorcts[keep,]
y <- DGEList(counts = cts_filt) |> calcNormFactors(method = "TMM")

subtumormeta$position %<>% as.factor()
subtumormeta$Patient_number %<>% as.factor()
mm <- model.matrix(~position+Patient_number, data = subtumormeta)
v <- voomLmFit(counts = y, normalize.method = "none", sample.weights = T, 
               design = mm,
               plot = T, save.plot = , keep.EList = T, adaptive.span = T)
contrast.matrix <- makeContrasts(positioncore, levels = mm)
out <- v |> contrasts.fit(contrasts = contrast.matrix) |> eBayes() |> topTable(number = Inf, coef = 1)
out$target <- rownames(out)
ggplot() + 
  scattermore::geom_scattermore(data = out[out$adj.P.Val >= 0.05,], mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "grey", pointsize = 2) + 
  geom_point(data = out[out$adj.P.Val < 0.05 & out$logFC > 0,], mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "darkgreen", size = 3, shape = 16) + 
  geom_point(data = out[out$adj.P.Val < 0.05 & out$logFC < 0,], mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "lightgreen", size = 3, shape = 16) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggrepel::geom_text_repel(data = out[out$adj.P.Val < 0.05,], mapping = aes(x = logFC, y = -log10(adj.P.Val), label = target), min.segment.length = 0, box.padding = 0.25, max.overlaps = 15, size = 6) +
  ggthemes::theme_par() +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20))
ggsave(filename = "core_vs_border_volcano.pdf", width = 10, height = 10)

#!#!# SHOW FIBROBLAST DECONVOLUTION RESULTS #!#!#

