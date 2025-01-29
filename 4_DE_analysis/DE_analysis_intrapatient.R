## DE Analysis Intrapatient ##
### Cole Nawrocki ###

# Clearing environment and making sure we are in the correct micromamba environment.
rm(list = ls())
.libPaths()

# [1] "/Users/ccn22/micromamba/envs/geomx-env/lib/R/library"

# Packages
library(limma)
library(parallel)
library(multcomp)
library(dplyr)
library(stringr)
library(magrittr)
library(ggplot2)
library(ComplexHeatmap)
options(mc.cores = 10)
set.seed(2001)

# Volcano plot function
plot_volcano <- function(.deres) {
  .deres$significance <- "NS"
  .deres$significance[.deres$pvalue < 0.05] <- "p < 0.05"
  .deres$significance[.deres$fdr < 0.05] <- "fdr < 0.05"
  .deres$significance <- factor(.deres$significance, levels = c("NS", "p < 0.05", "fdr < 0.05"))
  vp <- ggplot(data = .deres) + 
    geom_hline(yintercept = -log10(0.05), linewidth = 0.2) +
    geom_vline(xintercept = 0, linewidth = 0.2) +
    scattermore::geom_scattermore(mapping = aes(x = coefficient, y = -log10(pvalue), color = significance), alpha = 0.75, stroke = 0.1, pointsize = 3) + 
    geom_text(data = data.frame(x = c(min(.deres$coefficient)+0.1, max(.deres$coefficient)-0.1), y = c(0, 0), group = str_split(.deres$contrast[1], pattern = " / ", simplify = T)[1,c(2,1)]), mapping = aes(x = x, y = y, label = group)) +
    ggrepel::geom_text_repel(data = filter(.deres, fdr < 0.05), mapping = aes(x = coefficient, y = -log10(pvalue), label = target), size = 2, max.overlaps = 30, segment.size = 0.2) +
    scale_color_manual(values = c("black", "orange", "red")) + 
    ggthemes::theme_few() + 
    theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) + 
    labs(title = .deres$contrast[1]) + 
    guides(color = guide_legend(override.aes = list(size = 5)))
  return(vp)
}

# Data
meta <- openxlsx::read.xlsx("meta_cleaned_v5.xlsx")
rownames(meta) <- meta$aoi_id
counts <- read.csv(file = "counts.csv", row.names = 1)

meta <- arrange(meta, Patient_number)
counts <- counts[,rownames(meta)]
meta_list <- meta |> split(meta$Patient_number)
cts_list <- t(counts) |> as.data.frame() |> split(meta$Patient_number)
cts_list <- lapply(cts_list, t)

run_cols <- c("1"="grey", "2"="black")
pt_cols <- ggprism::ggprism_data$colour_palettes$colors[1:8]
names(pt_cols) <- 1:8

# Cohort 1 of Comparisons ------------------------------------------------------
# tumor intra-patient 2: id vs ii --> limma
tmpcts <- cts_list[[2]]
tmpmeta <- meta_list[[2]]
tmpmeta$group <- paste(tmpmeta$AOI_code, tmpmeta$infiltration_category, sep = "_")
mm <- model.matrix(~1+group, data = tmpmeta)
tumor_ii <- mm[tmpmeta$group == "tumor_inflamed",] |> colMeans()
tumor_id <- mm[tmpmeta$group == "tumor_desert",] |> colMeans()

tmpcts <- edgeR::DGEList(counts = tmpcts)
tmpcts <- edgeR::calcNormFactors(object = tmpcts)

vtmp <- voom(counts = tmpcts, design = mm)
fit <- lmFit(object = vtmp, design = mm)
fit <- contrasts.fit(fit = fit, contrasts = (tumor_id - tumor_ii))
res <- eBayes(fit = fit) |> topTable(number = Inf)

res$target <- rownames(res)
res$contrast <- "tumor id / ii"

write.csv(res, "patient2_tumor_id_vs_ii_limma.csv")

res <- rename(res, fdr = adj.P.Val, coefficient = logFC, pvalue = P.Value)
v1 <- plot_volcano(res)

res <- res |> filter(fdr < 0.05)
tmpmeta$group <- case_when(tmpmeta$infiltration_category == "excluded" ~ "ie", 
                           tmpmeta$infiltration_category == "inflamed" ~ "ii", 
                           tmpmeta$infiltration_category == "desert" ~ "id")
if (dim(res)[1] > 0) {
  cs <- str_split(pattern = " / |\\+", simplify = T, string = res$contrast)[1,] |> gsub(pattern = "tumor ", replacement = "")
  infil_cols <- c("id" = "red", "ii" = "blue", "ie" = "green")[cs]
  idx <- tmpmeta[(tmpmeta$group %in% cs) & (tmpmeta$AOI_code == "tumor"),] |> arrange(group) |> pull(aoi_id)
  toplot <- res |> arrange(desc(coefficient)) %$% target
  
  ha <- HeatmapAnnotation(group = tmpmeta[idx,]$group, 
                          pt = tmpmeta[idx,]$Patient_number,
                          run = tmpmeta[idx,]$Run_number,
                          col = list("group"= infil_cols, 
                                     "pt" = pt_cols, 
                                     "run" = run_cols))
  mat <- vtmp[toplot, idx, drop = F] |> apply(1, scale) |> t()
  hm1 <- Heatmap(matrix = mat, 
                top_annotation = ha, 
                name = "Scaled\nlog2(CPM + 1)", 
                cluster_columns = F, cluster_rows = F,
                width = ncol(mat)*unit(3, "mm"), 
                height = nrow(mat)*unit(3, "mm"),
                rect_gp = gpar(col = "white", lwd = 0.1), 
                row_names_gp = gpar(fontsize = 5)
  )
}

dev.off()
# tumor intra-patient 3: id vs ii --> limma
tmpcts <- cts_list[[3]]
tmpmeta <- meta_list[[3]]
tmpmeta$group <- paste(tmpmeta$AOI_code, tmpmeta$infiltration_category, sep = "_")
mm <- model.matrix(~1+group, data = tmpmeta)
tumor_ii <- mm[tmpmeta$group == "tumor_inflamed",] |> colMeans()
tumor_id <- mm[tmpmeta$group == "tumor_desert",] |> colMeans()

tmpcts <- edgeR::DGEList(counts = tmpcts)
tmpcts <- edgeR::calcNormFactors(object = tmpcts)

vtmp <- voom(counts = tmpcts, design = mm)
fit <- lmFit(object = vtmp, design = mm)
fit <- contrasts.fit(fit = fit, contrasts = (tumor_id - tumor_ii))
res <- eBayes(fit = fit) |> topTable(number = Inf)

res$target <- rownames(res)
res$contrast <- "tumor id / ii"

write.csv(res, "patient3_tumor_id_vs_ii_limma.csv")

res <- rename(res, fdr = adj.P.Val, coefficient = logFC, pvalue = P.Value)
v2 <- plot_volcano(res)

res <- res |> filter(fdr < 0.05)
tmpmeta$group <- case_when(tmpmeta$infiltration_category == "excluded" ~ "ie", 
                           tmpmeta$infiltration_category == "inflamed" ~ "ii", 
                           tmpmeta$infiltration_category == "desert" ~ "id")
if (dim(res)[1] > 0) {
  cs <- str_split(pattern = " / |\\+", simplify = T, string = res$contrast)[1,] |> gsub(pattern = "tumor ", replacement = "")
  infil_cols <- c("id" = "red", "ii" = "blue", "ie" = "green")[cs]
  idx <- tmpmeta[(tmpmeta$group %in% cs) & (tmpmeta$AOI_code == "tumor"),] |> arrange(group) |> pull(aoi_id)
  toplot <- res |> arrange(desc(coefficient)) %$% target
  
  ha <- HeatmapAnnotation(group = tmpmeta[idx,]$group, 
                          pt = tmpmeta[idx,]$Patient_number,
                          run = tmpmeta[idx,]$Run_number,
                          col = list("group"= infil_cols, 
                                     "pt" = pt_cols, 
                                     "run" = run_cols))
  mat <- vtmp[toplot, idx, drop = F] |> apply(1, scale) |> t()
  hm2 <- Heatmap(matrix = mat, 
                top_annotation = ha, 
                name = "Scaled\nlog2(CPM + 1)", 
                cluster_columns = F, cluster_rows = F,
                width = ncol(mat)*unit(3, "mm"), 
                height = nrow(mat)*unit(3, "mm"),
                rect_gp = gpar(col = "white", lwd = 0.1), 
                row_names_gp = gpar(fontsize = 5)
  )
}

# tumor intra-patient 6: ie vs ii --> limma
tmpcts <- cts_list[[6]]
tmpmeta <- meta_list[[6]]
tmpmeta$group <- paste(tmpmeta$AOI_code, tmpmeta$infiltration_category, sep = "_")
mm <- model.matrix(~1+group, data = tmpmeta)
tumor_ii <- mm[tmpmeta$group == "tumor_inflamed",] |> colMeans()
tumor_ie <- mm[tmpmeta$group == "tumor_excluded",] |> colMeans()

tmpcts <- edgeR::DGEList(counts = tmpcts)
tmpcts <- edgeR::calcNormFactors(object = tmpcts)

vtmp <- voom(counts = tmpcts, design = mm)
fit <- lmFit(object = vtmp, design = mm)
fit <- contrasts.fit(fit = fit, contrasts = (tumor_ie - tumor_ii))
res <- eBayes(fit = fit) |> topTable(number = Inf)

res$target <- rownames(res)
res$contrast <- "tumor ie / ii"

write.csv(res, "patient6_tumor_ie_vs_ii_limma.csv")

res <- rename(res, fdr = adj.P.Val, coefficient = logFC, pvalue = P.Value)
v3 <- plot_volcano(res)

res <- res |> filter(fdr < 0.05)
tmpmeta$group <- case_when(tmpmeta$infiltration_category == "excluded" ~ "ie", 
                           tmpmeta$infiltration_category == "inflamed" ~ "ii", 
                           tmpmeta$infiltration_category == "desert" ~ "id")
if (dim(res)[1] > 0) {
  cs <- str_split(pattern = " / |\\+", simplify = T, string = res$contrast)[1,] |> gsub(pattern = "tumor ", replacement = "")
  infil_cols <- c("id" = "red", "ii" = "blue", "ie" = "green")[cs]
  idx <- tmpmeta[(tmpmeta$group %in% cs) & (tmpmeta$AOI_code == "tumor"),] |> arrange(group) |> pull(aoi_id)
  toplot <- res |> arrange(desc(coefficient)) %$% target
  
  ha <- HeatmapAnnotation(group = tmpmeta[idx,]$group, 
                          pt = tmpmeta[idx,]$Patient_number,
                          run = tmpmeta[idx,]$Run_number,
                          col = list("group"= infil_cols, 
                                     "pt" = pt_cols, 
                                     "run" = run_cols))
  mat <- vtmp[toplot, idx, drop = F] |> apply(1, scale) |> t()
  hm3 <- Heatmap(matrix = mat, 
                top_annotation = ha, 
                name = "Scaled\nlog2(CPM + 1)", 
                cluster_columns = F, cluster_rows = F,
                width = ncol(mat)*unit(3, "mm"), 
                height = nrow(mat)*unit(3, "mm"),
                rect_gp = gpar(col = "white", lwd = 0.1), 
                row_names_gp = gpar(fontsize = 5)
  )
}

# tumor intra-patient 8: ie vs ii --> limma
tmpcts <- cts_list[[8]]
tmpmeta <- meta_list[[8]]
tmpmeta$group <- paste(tmpmeta$AOI_code, tmpmeta$infiltration_category, sep = "_")
mm <- model.matrix(~1+group, data = tmpmeta)
tumor_ii <- mm[tmpmeta$group == "tumor_inflamed",] |> colMeans()
tumor_ie <- mm[tmpmeta$group == "tumor_excluded",] |> colMeans()

tmpcts <- edgeR::DGEList(counts = tmpcts)
tmpcts <- edgeR::calcNormFactors(object = tmpcts)

vtmp <- voom(counts = tmpcts, design = mm)
fit <- lmFit(object = vtmp, design = mm)
fit <- contrasts.fit(fit = fit, contrasts = (tumor_ie - tumor_ii))
res <- eBayes(fit = fit) |> topTable(number = Inf)

res$target <- rownames(res)
res$contrast <- "tumor ie / ii"

write.csv(res, "patient8_tumor_ie_vs_ii_limma.csv")

res <- rename(res, fdr = adj.P.Val, coefficient = logFC, pvalue = P.Value)
v4 <- plot_volcano(res)

res <- res |> filter(fdr < 0.05)
tmpmeta$group <- case_when(tmpmeta$infiltration_category == "excluded" ~ "ie", 
                           tmpmeta$infiltration_category == "inflamed" ~ "ii", 
                           tmpmeta$infiltration_category == "desert" ~ "id")
if (dim(res)[1] > 0) {
  cs <- str_split(pattern = " / |\\+", simplify = T, string = res$contrast)[1,] |> gsub(pattern = "tumor ", replacement = "")
  infil_cols <- c("id" = "red", "ii" = "blue", "ie" = "green")[cs]
  idx <- tmpmeta[(tmpmeta$group %in% cs) & (tmpmeta$AOI_code == "tumor"),] |> arrange(group) |> pull(aoi_id)
  toplot <- res |> arrange(desc(coefficient)) %$% target
  
  ha <- HeatmapAnnotation(group = tmpmeta[idx,]$group, 
                          pt = tmpmeta[idx,]$Patient_number,
                          run = tmpmeta[idx,]$Run_number,
                          col = list("group"= infil_cols, 
                                     "pt" = pt_cols, 
                                     "run" = run_cols))
  mat <- vtmp[toplot, idx, drop = F] |> apply(1, scale) |> t()
  hm4 <- Heatmap(matrix = mat, 
                top_annotation = ha, 
                name = "Scaled\nlog2(CPM + 1)", 
                cluster_columns = F, cluster_rows = F,
                width = ncol(mat)*unit(3, "mm"), 
                height = nrow(mat)*unit(3, "mm"),
                rect_gp = gpar(col = "white", lwd = 0.1), 
                row_names_gp = gpar(fontsize = 5)
  )
}

# Cohort 2 of Comparisons ------------------------------------------------------
# tumor intra-patient 1: ieb vs iec --> limma
tmpcts <- cts_list[[1]]
tmpmeta <- meta_list[[1]]
tmpmeta$group <- paste(tmpmeta$AOI_code, tmpmeta$infiltration_type, sep = "_")
mm <- model.matrix(~1+group, data = tmpmeta)
tumor_ieb <- mm[tmpmeta$group == "tumor_ieb",] |> colMeans()
tumor_iec <- mm[tmpmeta$group == "tumor_iec",] |> colMeans()

tmpcts <- edgeR::DGEList(counts = tmpcts)
tmpcts <- edgeR::calcNormFactors(object = tmpcts)

vtmp <- voom(counts = tmpcts, design = mm)
fit <- lmFit(object = vtmp, design = mm)
fit <- contrasts.fit(fit = fit, contrasts = (tumor_ieb - tumor_iec))
res <- eBayes(fit = fit) |> topTable(number = Inf)

res$target <- rownames(res)
res$contrast <- "tumor ieb / iec"

write.csv(res, "patient1_tumor_ieb_vs_iec_limma.csv")

res <- rename(res, fdr = adj.P.Val, coefficient = logFC, pvalue = P.Value)
v5 <- plot_volcano(res)

res <- res |> filter(fdr < 0.05)
tmpmeta$group <- tmpmeta$infiltration_type
if (dim(res)[1] > 0) {
  cs <- str_split(pattern = " / |\\+", simplify = T, string = res$contrast)[1,] |> gsub(pattern = "tumor ", replacement = "")
  infil_cols <- c("idb" = "red", "iib" = "blue", "ieb" = "green", "iec"="orange", "iic"="magenta")[cs]
  idx <- tmpmeta[(tmpmeta$group %in% cs) & (tmpmeta$AOI_code == "tumor"),] |> arrange(group) |> pull(aoi_id)
  toplot <- res |> arrange(desc(coefficient)) %$% target
  
  ha <- HeatmapAnnotation(group = tmpmeta[idx,]$group, 
                          pt = tmpmeta[idx,]$Patient_number,
                          run = tmpmeta[idx,]$Run_number,
                          col = list("group"= infil_cols, 
                                     "pt" = pt_cols, 
                                     "run" = run_cols))
  mat <- vtmp[toplot, idx, drop = F] |> apply(1, scale) |> t()
  hm5 <- Heatmap(matrix = mat, 
                top_annotation = ha, 
                name = "Scaled\nlog2(CPM + 1)", 
                cluster_columns = F, cluster_rows = F,
                width = ncol(mat)*unit(3, "mm"), 
                height = nrow(mat)*unit(3, "mm"),
                rect_gp = gpar(col = "white", lwd = 0.1), 
                row_names_gp = gpar(fontsize = 5)
  )
}

# tumor intra-patient 4: iib vs iic --> limma
tmpcts <- cts_list[[4]]
tmpmeta <- meta_list[[4]]
tmpmeta$group <- paste(tmpmeta$AOI_code, tmpmeta$infiltration_type, sep = "_")
mm <- model.matrix(~1+group, data = tmpmeta)
tumor_iib <- mm[tmpmeta$group == "tumor_iib",] |> colMeans()
tumor_iic <- mm[tmpmeta$group == "tumor_iic",] |> colMeans()

tmpcts <- edgeR::DGEList(counts = tmpcts)
tmpcts <- edgeR::calcNormFactors(object = tmpcts)

vtmp <- voom(counts = tmpcts, design = mm)
fit <- lmFit(object = vtmp, design = mm)
fit <- contrasts.fit(fit = fit, contrasts = (tumor_iib - tumor_iic))
res <- eBayes(fit = fit) |> topTable(number = Inf)

res$target <- rownames(res)
res$contrast <- "tumor iib / iic"

write.csv(res, "patient4_tumor_iib_vs_iic_limma.csv")

res <- rename(res, fdr = adj.P.Val, coefficient = logFC, pvalue = P.Value)
v6 <- plot_volcano(res)

res <- res |> filter(fdr < 0.05)
tmpmeta$group <- tmpmeta$infiltration_type
if (dim(res)[1] > 0) {
  cs <- str_split(pattern = " / |\\+", simplify = T, string = res$contrast)[1,] |> gsub(pattern = "tumor ", replacement = "")
  infil_cols <- c("idb" = "red", "iib" = "blue", "ieb" = "green", "iec"="orange", "iic"="magenta")[cs]
  idx <- tmpmeta[(tmpmeta$group %in% cs) & (tmpmeta$AOI_code == "tumor"),] |> arrange(group) |> pull(aoi_id)
  toplot <- res |> arrange(desc(coefficient)) %$% target
  
  ha <- HeatmapAnnotation(group = tmpmeta[idx,]$group, 
                          pt = tmpmeta[idx,]$Patient_number,
                          run = tmpmeta[idx,]$Run_number,
                          col = list("group"= infil_cols, 
                                     "pt" = pt_cols, 
                                     "run" = run_cols))
  mat <- vtmp[toplot, idx, drop = F] |> apply(1, scale) |> t()
  hm6 <- Heatmap(matrix = mat, 
                top_annotation = ha, 
                name = "Scaled\nlog2(CPM + 1)", 
                cluster_columns = F, cluster_rows = F,
                width = ncol(mat)*unit(3, "mm"), 
                height = nrow(mat)*unit(3, "mm"),
                rect_gp = gpar(col = "white", lwd = 0.1), 
                row_names_gp = gpar(fontsize = 5)
  )
}

# tumor intra-patient 5: iib vs iic --> limma
tmpcts <- cts_list[[5]]
tmpmeta <- meta_list[[5]]
tmpmeta$group <- paste(tmpmeta$AOI_code, tmpmeta$infiltration_type, sep = "_")
mm <- model.matrix(~1+group, data = tmpmeta)
tumor_iib <- mm[tmpmeta$group == "tumor_iib",] |> colMeans()
tumor_iic <- mm[tmpmeta$group == "tumor_iic",] |> colMeans()

tmpcts <- edgeR::DGEList(counts = tmpcts)
tmpcts <- edgeR::calcNormFactors(object = tmpcts)

vtmp <- voom(counts = tmpcts, design = mm)
fit <- lmFit(object = vtmp, design = mm)
fit <- contrasts.fit(fit = fit, contrasts = (tumor_iib - tumor_iic))
res <- eBayes(fit = fit) |> topTable(number = Inf)

res$target <- rownames(res)
res$contrast <- "tumor iib / iic"

write.csv(res, "patient5_tumor_iib_vs_iic_limma.csv")

res <- rename(res, fdr = adj.P.Val, coefficient = logFC, pvalue = P.Value)
v7 <- plot_volcano(res)

res <- res |> filter(fdr < 0.05)
tmpmeta$group <- tmpmeta$infiltration_type
if (dim(res)[1] > 0) {
  cs <- str_split(pattern = " / |\\+", simplify = T, string = res$contrast)[1,] |> gsub(pattern = "tumor ", replacement = "")
  infil_cols <- c("idb" = "red", "iib" = "blue", "ieb" = "green", "iec"="orange", "iic"="magenta")[cs]
  idx <- tmpmeta[(tmpmeta$group %in% cs) & (tmpmeta$AOI_code == "tumor"),] |> arrange(group) |> pull(aoi_id)
  toplot <- res |> arrange(desc(coefficient)) %$% target
  
  ha <- HeatmapAnnotation(group = tmpmeta[idx,]$group, 
                          pt = tmpmeta[idx,]$Patient_number,
                          run = tmpmeta[idx,]$Run_number,
                          col = list("group"= infil_cols, 
                                     "pt" = pt_cols, 
                                     "run" = run_cols))
  mat <- vtmp[toplot, idx, drop = F] |> apply(1, scale) |> t()
  hm7 <- Heatmap(matrix = mat, 
                top_annotation = ha, 
                name = "Scaled\nlog2(CPM + 1)", 
                cluster_columns = F, cluster_rows = F,
                width = ncol(mat)*unit(3, "mm"), 
                height = nrow(mat)*unit(3, "mm"),
                rect_gp = gpar(col = "white", lwd = 0.1), 
                row_names_gp = gpar(fontsize = 5)
  )
}

# tumor intra-patient 6: ieb vs iec --> limma
tmpcts <- cts_list[[6]]
tmpmeta <- meta_list[[6]]
tmpmeta$group <- paste(tmpmeta$AOI_code, tmpmeta$infiltration_type, sep = "_")
mm <- model.matrix(~1+group, data = tmpmeta)
tumor_ieb <- mm[tmpmeta$group == "tumor_ieb",] |> colMeans()
tumor_iec <- mm[tmpmeta$group == "tumor_iec",] |> colMeans()

tmpcts <- edgeR::DGEList(counts = tmpcts)
tmpcts <- edgeR::calcNormFactors(object = tmpcts)

vtmp <- voom(counts = tmpcts, design = mm)
fit <- lmFit(object = vtmp, design = mm)
fit <- contrasts.fit(fit = fit, contrasts = (tumor_ieb - tumor_iec))
res <- eBayes(fit = fit) |> topTable(number = Inf)

res$target <- rownames(res)
res$contrast <- "tumor ieb / iec"

write.csv(res, "patient6_tumor_ieb_vs_iec_limma.csv")

res <- rename(res, fdr = adj.P.Val, coefficient = logFC, pvalue = P.Value)
v8 <- plot_volcano(res)

res <- res |> filter(fdr < 0.05)
tmpmeta$group <- tmpmeta$infiltration_type
if (dim(res)[1] > 0) {
  cs <- str_split(pattern = " / |\\+", simplify = T, string = res$contrast)[1,] |> gsub(pattern = "tumor ", replacement = "")
  infil_cols <- c("idb" = "red", "iib" = "blue", "ieb" = "green", "iec"="orange", "iic"="magenta")[cs]
  idx <- tmpmeta[(tmpmeta$group %in% cs) & (tmpmeta$AOI_code == "tumor"),] |> arrange(group) |> pull(aoi_id)
  toplot <- res |> arrange(desc(coefficient)) %$% target
  
  ha <- HeatmapAnnotation(group = tmpmeta[idx,]$group, 
                          pt = tmpmeta[idx,]$Patient_number,
                          run = tmpmeta[idx,]$Run_number,
                          col = list("group"= infil_cols, 
                                     "pt" = pt_cols, 
                                     "run" = run_cols))
  mat <- vtmp[toplot, idx, drop = F] |> apply(1, scale) |> t()
  hm8 <- Heatmap(matrix = mat, 
                top_annotation = ha, 
                name = "Scaled\nlog2(CPM + 1)", 
                cluster_columns = F, cluster_rows = F,
                width = ncol(mat)*unit(3, "mm"), 
                height = nrow(mat)*unit(3, "mm"),
                rect_gp = gpar(col = "white", lwd = 0.1), 
                row_names_gp = gpar(fontsize = 5)
  )
}

# tumor intra-patient 7: ieb vs iec --> limma
tmpcts <- cts_list[[7]]
tmpmeta <- meta_list[[7]]
tmpmeta$group <- paste(tmpmeta$AOI_code, tmpmeta$infiltration_type, sep = "_")
mm <- model.matrix(~1+group, data = tmpmeta)
tumor_ieb <- mm[tmpmeta$group == "tumor_ieb",] |> colMeans()
tumor_iec <- mm[tmpmeta$group == "tumor_iec",] |> colMeans()

tmpcts <- edgeR::DGEList(counts = tmpcts)
tmpcts <- edgeR::calcNormFactors(object = tmpcts)

vtmp <- voom(counts = tmpcts, design = mm)
fit <- lmFit(object = vtmp, design = mm)
fit <- contrasts.fit(fit = fit, contrasts = (tumor_ieb - tumor_iec))
res <- eBayes(fit = fit) |> topTable(number = Inf)

res$target <- rownames(res)
res$contrast <- "tumor ieb / iec"

write.csv(res, "patient7_tumor_ieb_vs_iec_limma.csv")

res <- rename(res, fdr = adj.P.Val, coefficient = logFC, pvalue = P.Value)
v9 <- plot_volcano(res)

res <- res |> filter(fdr < 0.05)
tmpmeta$group <- tmpmeta$infiltration_type
if (dim(res)[1] > 0) {
  cs <- str_split(pattern = " / |\\+", simplify = T, string = res$contrast)[1,] |> gsub(pattern = "tumor ", replacement = "")
  infil_cols <- c("idb" = "red", "iib" = "blue", "ieb" = "green", "iec"="orange", "iic"="magenta")[cs]
  idx <- tmpmeta[(tmpmeta$group %in% cs) & (tmpmeta$AOI_code == "tumor"),] |> arrange(group) |> pull(aoi_id)
  toplot <- res |> arrange(desc(coefficient)) %$% target
  
  ha <- HeatmapAnnotation(group = tmpmeta[idx,]$group, 
                          pt = tmpmeta[idx,]$Patient_number,
                          run = tmpmeta[idx,]$Run_number,
                          col = list("group"= infil_cols, 
                                     "pt" = pt_cols, 
                                     "run" = run_cols))
  mat <- vtmp[toplot, idx, drop = F] |> apply(1, scale) |> t()
  hm9 <- Heatmap(matrix = mat, 
                top_annotation = ha, 
                name = "Scaled\nlog2(CPM + 1)", 
                cluster_columns = F, cluster_rows = F,
                width = ncol(mat)*unit(3, "mm"), 
                height = nrow(mat)*unit(3, "mm"),
                rect_gp = gpar(col = "white", lwd = 0.1), 
                row_names_gp = gpar(fontsize = 5)
  )
}

# tumor intra-patient 8: ieb vs iec --> limma
tmpcts <- cts_list[[8]]
tmpmeta <- meta_list[[8]]
tmpmeta$group <- paste(tmpmeta$AOI_code, tmpmeta$infiltration_type, sep = "_")
mm <- model.matrix(~1+group, data = tmpmeta)
tumor_ieb <- mm[tmpmeta$group == "tumor_ieb",] |> colMeans()
tumor_iec <- mm[tmpmeta$group == "tumor_iec",] |> colMeans()

tmpcts <- edgeR::DGEList(counts = tmpcts)
tmpcts <- edgeR::calcNormFactors(object = tmpcts)

vtmp <- voom(counts = tmpcts, design = mm)
fit <- lmFit(object = vtmp, design = mm)
fit <- contrasts.fit(fit = fit, contrasts = (tumor_ieb - tumor_iec))
res <- eBayes(fit = fit) |> topTable(number = Inf)

res$target <- rownames(res)
res$contrast <- "tumor ieb / iec"

write.csv(res, "patient8_tumor_ieb_vs_iec_limma.csv")

res <- rename(res, fdr = adj.P.Val, coefficient = logFC, pvalue = P.Value)
v10 <- plot_volcano(res)

res <- res |> filter(fdr < 0.05)
tmpmeta$group <- tmpmeta$infiltration_type
if (dim(res)[1] > 0) {
  cs <- str_split(pattern = " / |\\+", simplify = T, string = res$contrast)[1,] |> gsub(pattern = "tumor ", replacement = "")
  infil_cols <- c("idb" = "red", "iib" = "blue", "ieb" = "green", "iec"="orange", "iic"="magenta")[cs]
  idx <- tmpmeta[(tmpmeta$group %in% cs) & (tmpmeta$AOI_code == "tumor"),] |> arrange(group) |> pull(aoi_id)
  toplot <- res |> arrange(desc(coefficient)) %$% target
  
  ha <- HeatmapAnnotation(group = tmpmeta[idx,]$group, 
                          pt = tmpmeta[idx,]$Patient_number,
                          run = tmpmeta[idx,]$Run_number,
                          col = list("group"= infil_cols, 
                                     "pt" = pt_cols, 
                                     "run" = run_cols))
  mat <- vtmp[toplot, idx, drop = F] |> apply(1, scale) |> t()
  hm10 <- Heatmap(matrix = mat, 
                top_annotation = ha, 
                name = "Scaled\nlog2(CPM + 1)", 
                cluster_columns = F, cluster_rows = F,
                width = ncol(mat)*unit(3, "mm"), 
                height = nrow(mat)*unit(3, "mm"),
                rect_gp = gpar(col = "white", lwd = 0.1), 
                row_names_gp = gpar(fontsize = 5)
  )
}

# Plotting ---------------------------------------------------------------------
pdf(file = "intrapatient_comparisons_volcanos_limma.pdf", width = 8, height = 8)
print(v1)
print(v2)
print(v3)
print(v4)
print(v5)
print(v6)
print(v7)
print(v8)
print(v9)
print(v10)
dev.off()

pdf(file = "intrapatient_comparisons_heatmaps_limma.pdf", width = 12, height = 24)
print(hm1)
draw(hm2)
#draw(hm3)
#draw(hm4)
draw(hm5)
draw(hm6)
#draw(hm7)
#draw(hm8)
draw(hm9)
#draw(hm10)
dev.off()

# This took about 2 hours

