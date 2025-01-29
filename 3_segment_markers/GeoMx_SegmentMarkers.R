###### GeoMx Data Segment Markers Script ######
## Cole Nawrocki ##

# Check that we are in the correct environment
.libPaths()
# [1] "/Users/ccn22/micromamba/envs/geomx-env/lib/R/library"

# Packages
library(dplyr)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(limma)
library(edgeR)
library(ggbiplot)
library(magrittr)
library(patchwork)
library(lme4)
library(lmerTest)
library(parallel)
options(mc.cores = 8)

# Data: we will use the q3 normalized expression.
meta <- read.csv("meta.csv", row.names = 1)
norm <- read.csv("q3norm.csv", row.names = 1)

segment_cols <- c("tumor" = "green4", "immune" = "red", "fibroblast" = "dodgerblue")
run_cols <- c("R1"="yellow", "R2"="orange")
slide_cols <- grDevices::colors()[grep('(gr(a|e)y)|(white)', grDevices::colors(), invert = T)] |> sample(size = 8, replace = F)
names(slide_cols) <- unique(meta$slide.name)

pdf("Segment_Markers.pdf", width = 10, height = 8)

### Segment Markers ###
# Fitting the model and performing DE
mm <- model.matrix(~segment, data = meta)
contrast_list <- list()
for (seg in unique(meta$segment)) {
  contrast_list[[seg]] <- mm[meta$segment == seg,] |> colMeans()
}

lmem_de <- function(gene) {
  d <- dplyr::mutate(meta, expr=(norm[gene,] |> unlist()))
  modelout <- lmerTest::lmer(formula = expr ~ segment + (1+segment|slide.name), data = d)
  out <- list()
  for (seg in names(contrast_list)) {
    seg_oi <- contrast_list[[seg]]
    other <- Reduce(f = "+", 
                    x = contrast_list[names(contrast_list) != seg]) /
      (length(contrast_list[names(contrast_list) != seg]))
    seg_oi_vs_other <- seg_oi - other
    out[[seg]] <- contest(model = modelout, L = seg_oi_vs_other, joint = F, rhs = 0)
  }
  out <- bind_rows(out, .id = "group")
  out$gene <- gene
  return(out)
}

de <- parallel::mclapply(rownames(norm), lmem_de, mc.cores = 8)
dedf <- bind_rows(de)
dedf <- group_by(dedf, group) |> mutate(fdr = p.adjust(method = "fdr", `Pr(>|t|)`))

# Overall volcano
ggplot() + 
  geom_point(data = dedf, mapping = aes(x = Estimate, y = -log10(fdr))) + 
  geom_hline(yintercept = -log10(0.05)) + 
  geom_vline(xintercept = c(-1,1))

# Heatmap
toplot <- dedf |> 
  filter(fdr < 0.05) |> 
  group_by(group) |> 
  top_n(n = 25, wt = Estimate) |> 
  arrange(group, desc(Estimate))

# Heatmap
idx <- meta |> arrange(segment) |> rownames()
mat <- norm[toplot$gene, idx] |> apply(1, scale) |> t()
ha <- HeatmapAnnotation(segment = meta[idx, ]$segment, 
                        slide =  meta[idx, ]$slide.name, 
                        run = meta[idx, ]$Run, 
                        col = list("segment" = segment_cols, 
                                   "slide" = slide_cols, 
                                   "run" = run_cols))
ha2 <- rowAnnotation(marker_for = toplot$group, 
                     col = list("marker_for" = segment_cols))
Heatmap(matrix = mat, 
        top_annotation = ha, 
        right_annotation = ha2,
        row_names_gp = gpar(fontsize = 5), 
        cluster_rows = F, 
        cluster_columns = F, 
        name = "Scaled Q3\nNormalized\nExpression", 
        column_split = meta[idx, ]$segment
)

# Volcanos
# Tumor
tumor_res <- dedf |> filter(group == "tumor")
tumor_res$contrast <- "neither"
tumor_res$contrast <- ifelse((tumor_res$Estimate > 1) & (tumor_res$fdr < 0.05), yes = "tumor", no = tumor_res$contrast)
tumor_res$contrast <- ifelse((tumor_res$Estimate < -1) & (tumor_res$fdr < 0.05), yes = "other", no = tumor_res$contrast)

ggplot() + 
  geom_point(data = tumor_res, mapping = aes(x = Estimate, y = -log10(fdr), color = contrast), alpha = 0.5) + 
  ggprism::theme_prism() + 
  ggrepel::geom_text_repel(data = tumor_res |> filter(abs(Estimate) > 1 & fdr < 0.05), 
                  mapping = aes(x = Estimate, y = -log10(fdr), label = gene, color = contrast), size = 3) + 
  scale_color_manual(values = c("other"="purple", "tumor"="green4")) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dotted") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dotted")

# Fibroblasts
fibroblast_res <- dedf |> filter(group == "fibroblast")
fibroblast_res$contrast <- "neither"
fibroblast_res$contrast <- ifelse((fibroblast_res$Estimate > 1) & (fibroblast_res$fdr < 0.05), yes = "fibroblast", no = fibroblast_res$contrast)
fibroblast_res$contrast <- ifelse((fibroblast_res$Estimate < -1) & (fibroblast_res$fdr < 0.05), yes = "other", no = fibroblast_res$contrast)

ggplot() + 
  geom_point(data = fibroblast_res, mapping = aes(x = Estimate, y = -log10(fdr), color = contrast), alpha = 0.5) + 
  ggprism::theme_prism() + 
  ggrepel::geom_text_repel(data = fibroblast_res |> filter(abs(Estimate) > 1 & fdr < 0.05), 
                           mapping = aes(x = Estimate, y = -log10(fdr), label = gene, color = contrast), size = 3) + 
  scale_color_manual(values = c("other"="brown", "fibroblast"="dodgerblue")) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dotted") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dotted")

# Immune
immune_res <- dedf |> filter(group == "immune")
immune_res$contrast <- "neither"
immune_res$contrast <- ifelse((immune_res$Estimate > 1) & (immune_res$fdr < 0.05), yes = "immune", no = immune_res$contrast)
immune_res$contrast <- ifelse((immune_res$Estimate < -1) & (immune_res$fdr < 0.05), yes = "other", no = immune_res$contrast)

ggplot() + 
  geom_point(data = immune_res, mapping = aes(x = Estimate, y = -log10(fdr), color = contrast), alpha = 0.5) + 
  ggprism::theme_prism() + 
  ggrepel::geom_text_repel(data = immune_res |> filter(abs(Estimate) > 1 & fdr < 0.05), 
                           mapping = aes(x = Estimate, y = -log10(fdr), label = gene, color = contrast), size = 3) + 
  scale_color_manual(values = c("other"="turquoise", "immune"="red")) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dotted") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dotted")

dev.off()
