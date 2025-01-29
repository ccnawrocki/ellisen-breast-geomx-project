## DE Analysis ##
### Cole Nawrocki ###

# Clearing environment and making sure we are in the correct micromamba environment.
rm(list = ls())
.libPaths()

# [1] "/Users/ccn22/micromamba/envs/geomx-env/lib/R/library"

# Packages
library(lme4)
library(lmerTest)
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

# Testing function
lmem_de_test <- function(outmod, .contr, .test_type = "univariate") {
  if (.test_type == "univariate") {
    testout <- outmod |> 
      glht(linfct = matrix(contrast_list[[.contr]], nrow = 1), alternative = "two.sided", rhs = 0) |> 
      summary(test = univariate()) %$%test
    out <- data.frame(
      contrast = .contr,
      coefficient = testout$coefficients[[1]],
      sigma = testout$sigma[[1]],
      tstat = testout$tstat[[1]],
      pvalue = testout$pvalues[[1]],
      test_type = testout$type
    )
  }
  
  if (.test_type == "Ftest") {
    testout <- outmod |> 
      glht(linfct = matrix(contrast_list[[.contr]], nrow = 1), alternative = "two.sided", rhs = 0) |> 
      summary(test = Ftest()) %$%test
    out <- data.frame(
      contrast = .contr,
      coefficient = testout$coefficients[[1]],
      sigma = testout$sigma[[1]],
      df1 = testout$df[1],
      df2 = testout$df[2],
      fstat = testout$fstat,
      tstat = testout$tstat[[1]],
      pvalue = testout$pvalue,
      test_type = testout$type
    )
  }
  
  if (.test_type == "Chisqtest") {
    testout <- outmod |> 
      glht(linfct = matrix(contrast_list[[.contr]], nrow = 1), alternative = "two.sided", rhs = 0) |> 
      summary(test = Chisqtest()) %$%test
    out <- data.frame(
      contrast = .contr,
      coefficient = testout$coefficients[[1]],
      sigma = testout$sigma[[1]],
      df = 1,
      fstat = testout$fstat,
      tstat = testout$tstat[[1]],
      pvalue = testout$pvalue,
      test_type = testout$type
    )
  }
  
  if (.test_type == "Satterthwaite") {
    out <- lmerTest::contest(model = outmod, L = contrast_list[[.contr]], rhs = 0, joint = F, confint = F, ddf = "Satterthwaite")
    out <- cbind("contrast" = .contr, "type" = "Satterthwaite", out)
  }
  
  if (.test_type == "Kenward-Roger") {
    out <- lmerTest::contest(model = outmod, L = contrast_list[[.contr]], rhs = 0, joint = F, confint = F, ddf = "Kenward-Roger")
    out <- cbind("contrast" = .contr, "type" = "Kenward-Roger", out)
  }
  
  return(out)
}

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
norm <- read.csv(file = "q3norm.csv", row.names = 1)
meta <- openxlsx::read.xlsx("meta_cleaned_v5.xlsx")
rownames(meta) <- meta$aoi_id
counts <- read.csv(file = "counts.csv", row.names = 1)

# Mixed models
simple_res <- readRDS("4_DE_analysis/lmem_de_simple_model_results_ml.RDS")
minus_idc_res <- readRDS("4_DE_analysis/lmem_de_model_results_ml_minus_idc.RDS")

## Comparison Cohort 1  --------------------------------------------------------
# Please ignore border and core, so for following comparison, please always combine border and core.
# tumor: ii vs. ie+id --> simple lme4 model, limma (patient as FE)
# tumor: ie vs. ii+id --> simple lme4 model, limma (patient as FE)
# tumor: id vs. ii+ie --> simple lme4 model, limma (patient as FE)
# tumor: ie vs. ii --> simple lme4 model, limma (patient as FE)
# tumor: id vs. ii --> simple lme4 model, limma (patient as FE)
# tumor: ie vs. id --> simple lme4 model, limma (patient as FE)

### MIXED MODEL ----------------------------------------------------------------
meta$group <- paste(meta$AOI_code, meta$infiltration_category, sep = "_")
mm <- model.matrix(~1+group, data = meta)
tumor_ii <- mm[meta$group == "tumor_inflamed",] |> colMeans()
tumor_ie <- mm[meta$group == "tumor_excluded",] |> colMeans()
tumor_id <- mm[meta$group == "tumor_desert",] |> colMeans()

contrast_list <- list(
  "tumor ii / ie+id" = tumor_ii - (tumor_ie + tumor_id)/2, 
  "tumor ie / ii+id" = tumor_ie - (tumor_ii + tumor_id)/2,
  "tumor id / ii+ie" = tumor_id - (tumor_ii + tumor_ie)/2, 
  "tumor ie / ii" = tumor_ie - tumor_ii,
  "tumor id / ii" = tumor_id - tumor_ii, 
  "tumor ie / id" = tumor_ie - tumor_id
)

pdf(file = "comparisons_first_cohort_volcanos_mixed_model.pdf", width = 8, height = 8)
for (c in names(contrast_list)) {
  restest <- lapply(X = simple_res, FUN = lmem_de_test, .contr = c)
  restest <- bind_rows(restest, .id = "target")
  restest$fdr <- p.adjust(restest$pvalue, method = "fdr")
  write.csv(restest, paste(gsub(pattern = "/", replacement = "vs", x = c), " mixed model.csv", sep = ""))
  print(plot_volcano(restest))
}
dev.off()

run_cols <- c("1"="grey", "2"="black")
pt_cols <- ggprism::ggprism_data$colour_palettes$colors[1:8]
names(pt_cols) <- 1:8

pdf(file = "comparisons_first_cohort_heatmaps_mixed_model.pdf", width = 24, height = 8)
meta$group <- case_when(meta$infiltration_category == "excluded" ~ "ie", 
                        meta$infiltration_category == "inflamed" ~ "ii", 
                        meta$infiltration_category == "desert" ~ "id")
for (c in names(contrast_list)) {
  restest <- read.csv(paste(gsub(pattern = "/", replacement = "vs", x = c), " mixed model.csv", sep = ""), row.names = 1)
  restest <- restest |> filter(fdr < 0.05)
  if (dim(restest)[1] > 0) {
    cs <- str_split(pattern = " / |\\+", simplify = T, string = c)[1,] |> gsub(pattern = "tumor ", replacement = "")
    infil_cols <- c("id" = "red", "ii" = "blue", "ie" = "green")[cs]
    idx <- meta[(meta$group %in% cs) & (meta$AOI_code == "tumor"),] |> arrange(group) |> pull(aoi_id)
    toplot <- restest |> arrange(desc(coefficient)) %$% target
    
    ha <- HeatmapAnnotation(group = meta[idx,]$group, 
                            pt = meta[idx,]$Patient_number,
                            run = meta[idx,]$Run_number,
                            col = list("group"= infil_cols, 
                                       "pt" = pt_cols, 
                                       "run" = run_cols))
    mat <- norm[toplot, idx, drop = F] |> apply(1, scale) |> t()
    hm <- Heatmap(matrix = mat, 
                  top_annotation = ha, 
                  name = "Scaled\nlog2(q3norm + 1)", 
                  cluster_columns = F, cluster_rows = F,
                  width = ncol(mat)*unit(3, "mm"), 
                  height = nrow(mat)*unit(3, "mm"),
                  rect_gp = gpar(col = "white", lwd = 0.1), 
                  row_names_gp = gpar(fontsize = 6),
                  column_title = c
                  )
    draw(hm)
  }
}

dev.off()

### LIMMA ----------------------------------------------------------------------
cts <- edgeR::DGEList(counts = counts)
cts <- edgeR::calcNormFactors(cts)

meta$group <- paste(meta$AOI_code, meta$infiltration_category, sep = "_")
mm2 <- model.matrix(~1+group+Patient_number, data = meta)
tumor_ii <- mm2[meta$group == "tumor_inflamed",] |> colMeans()
tumor_ie <- mm2[meta$group == "tumor_excluded",] |> colMeans()
tumor_id <- mm2[meta$group == "tumor_desert",] |> colMeans()

contrast_list <- list(
  "tumor ii / ie+id" = tumor_ii - (tumor_ie + tumor_id)/2, 
  "tumor ie / ii+id" = tumor_ie - (tumor_ii + tumor_id)/2,
  "tumor id / ii+ie" = tumor_id - (tumor_ii + tumor_ie)/2, 
  "tumor ie / ii" = tumor_ie - tumor_ii,
  "tumor id / ii" = tumor_id - tumor_ii, 
  "tumor ie / id" = tumor_ie - tumor_id
)

vn <- voom(counts = cts, design = mm2)
limmafit <- lmFit(object = vn, design = mm2)

pdf(file = "comparisons_first_cohort_volcanos_limma.pdf", width = 8, height = 8)
for (c in names(contrast_list)) {
  deout <- contrasts.fit(limmafit, contrasts = contrast_list[[c]])
  deout <- eBayes(deout) |> topTable(number = Inf)
  deout$target <- rownames(deout)
  deout$contrast <- c
  
  write.csv(deout, paste(gsub(pattern = "/", replacement = "vs", x = c), " limma.csv", sep = ""))
  
  deout <- rename(deout, fdr = adj.P.Val, coefficient = logFC, pvalue = P.Value)
  print(plot_volcano(deout))
}
dev.off()

run_cols <- c("1"="grey", "2"="black")
pt_cols <- ggprism::ggprism_data$colour_palettes$colors[1:8]
names(pt_cols) <- 1:8

pdf(file = "comparisons_first_cohort_heatmaps_limma.pdf", width = 25, height = 50)
meta$group <- case_when(meta$infiltration_category == "excluded" ~ "ie", 
                        meta$infiltration_category == "inflamed" ~ "ii", 
                        meta$infiltration_category == "desert" ~ "id")
for (c in names(contrast_list)) {
  restest <- read.csv(paste(gsub(pattern = "/", replacement = "vs", x = c), " limma.csv", sep = ""), row.names = 1)
  restest <- rename(restest, fdr = adj.P.Val, coefficient = logFC, pvalue = P.Value)
  restest <- restest |> filter(fdr < 0.05)
  if (dim(restest)[1] > 0) {
    cs <- str_split(pattern = " / |\\+", simplify = T, string = c)[1,] |> gsub(pattern = "tumor ", replacement = "")
    infil_cols <- c("id" = "red", "ii" = "blue", "ie" = "green")[cs]
    idx <- meta[(meta$group %in% cs) & (meta$AOI_code == "tumor"),] |> arrange(group) |> pull(aoi_id)
    toplot <- restest |> arrange(desc(coefficient)) %$% target
    
    ha <- HeatmapAnnotation(group = meta[idx,]$group, 
                            pt = meta[idx,]$Patient_number,
                            run = meta[idx,]$Run_number,
                            col = list("group"= infil_cols, 
                                       "pt" = pt_cols, 
                                       "run" = run_cols))
    mat <- vn[toplot, idx, drop = F] |> apply(1, scale) |> t()
    hm <- Heatmap(matrix = mat, 
                  top_annotation = ha, 
                  name = "Scaled\nlog2(CPM + 1)", 
                  cluster_columns = F, cluster_rows = F,
                  width = ncol(mat)*unit(3, "mm"), 
                  height = nrow(mat)*unit(1.5, "mm"),
                  rect_gp = gpar(col = "white", lwd = 0.1), 
                  row_names_gp = gpar(fontsize = 3),
                  column_title = c
    )
    draw(hm)
  }
}

dev.off()

## Comparison Cohort 2 ---------------------------------------------------------
# Take core vs border into consideration, please do not combine border and core in any cases.
# tumor ieb vs. iib --> complex lme4 model, limma (patient as FE)
# tumor idb vs. iib --> complex lme4 model, limma (patient as FE)
# tumor ieb vs. idb --> complex lme4 model, limma (patient as FE)
# tumor iec vs. iic --> complex lme4 model, limma (patient as FE)
# tumor ieb vs. iec --> complex lme4 model, limma (patient as FE)
# tumor iib vs. iic --> complex lme4 model, limma (patient as FE)

### MIXED MODEL ----------------------------------------------------------------
meta$infiltration_group <- meta$infiltration_type
meta$infiltration_group <- ifelse(test = meta$infiltration_category == "desert", 
                                  yes = str_sub(meta$infiltration_group, start = 1, end = 2), 
                                  no = meta$infiltration_group)
meta$group <- paste(meta$AOI_code, meta$infiltration_group, sep = "_")
meta <- meta |> filter(infiltration_type != "idc")
mm <- model.matrix(~1+group, data = meta)
tumor_ieb <- mm[meta$group == "tumor_ieb",] |> colMeans()
tumor_iib <- mm[meta$group == "tumor_iib",] |> colMeans()
tumor_idb <- mm[meta$group == "tumor_id",] |> colMeans()
tumor_iec <- mm[meta$group == "tumor_iec",] |> colMeans()
tumor_iic <- mm[meta$group == "tumor_iic",] |> colMeans()

contrast_list <- list(
  "tumor ieb / iib" = tumor_ieb - tumor_iib, 
  "tumor idb / iib" = tumor_idb - tumor_iib, 
  "tumor ieb / idb" = tumor_ieb - tumor_idb, 
  "tumor iec / iic" = tumor_iec - tumor_iic, 
  "tumor ieb / iec" = tumor_ieb - tumor_iec, 
  "tumor iib / iic" = tumor_iib - tumor_iic
)

pdf(file = "comparisons_second_cohort_volcanos_mixed_model.pdf", width = 8, height = 8)
for (c in names(contrast_list)) {
  restest <- lapply(X = minus_idc_res, FUN = lmem_de_test, .contr = c)
  restest <- bind_rows(restest, .id = "target")
  restest$fdr <- p.adjust(restest$pvalue, method = "fdr")
  write.csv(restest, paste(gsub(pattern = "/", replacement = "vs", x = c), " mixed model.csv", sep = ""))
  print(plot_volcano(restest))
}
dev.off()

run_cols <- c("1"="grey", "2"="black")
pt_cols <- ggprism::ggprism_data$colour_palettes$colors[1:8]
names(pt_cols) <- 1:8

pdf(file = "comparisons_second_cohort_heatmaps_mixed_model.pdf", width = 24, height = 8)
meta$group <- meta$infiltration_type
for (c in names(contrast_list)) {
  restest <- read.csv(paste(gsub(pattern = "/", replacement = "vs", x = c), " mixed model.csv", sep = ""), row.names = 1)
  restest <- restest |> filter(fdr < 0.05)
  if (dim(restest)[1] > 0) {
    cs <- str_split(pattern = " / |\\+", simplify = T, string = c)[1,] |> gsub(pattern = "tumor ", replacement = "")
    infil_cols <- c("idb" = "red", "iib" = "blue", "ieb" = "green", "iec"="orange", "iic"="magenta")[cs]
    idx <- meta[(meta$group %in% cs) & (meta$AOI_code == "tumor"),] |> arrange(group) |> pull(aoi_id)
    toplot <- restest |> arrange(desc(coefficient)) %$% target
    
    ha <- HeatmapAnnotation(group = meta[idx,]$group, 
                            pt = meta[idx,]$Patient_number,
                            run = meta[idx,]$Run_number,
                            col = list("group"= infil_cols, 
                                       "pt" = pt_cols, 
                                       "run" = run_cols))
    mat <- norm[toplot, idx, drop = F] |> apply(1, scale) |> t()
    hm <- Heatmap(matrix = mat, 
                  top_annotation = ha, 
                  name = "Scaled\nlog2(q3norm + 1)", 
                  cluster_columns = F, cluster_rows = F,
                  width = ncol(mat)*unit(3, "mm"), 
                  height = nrow(mat)*unit(3, "mm"),
                  rect_gp = gpar(col = "white", lwd = 0.1), 
                  row_names_gp = gpar(fontsize = 6),
                  column_title = c
    )
    draw(hm)
  }
}

dev.off()

### LIMMA ----------------------------------------------------------------------
cts <- edgeR::DGEList(counts = counts[,rownames(meta)])
cts <- edgeR::calcNormFactors(cts)

meta$group <- paste(meta$AOI_code, meta$infiltration_group, sep = "_")
mm2 <- model.matrix(~1+group+Patient_number, data = meta)
tumor_ieb <- mm2[meta$group == "tumor_ieb",] |> colMeans()
tumor_iib <- mm2[meta$group == "tumor_iib",] |> colMeans()
tumor_idb <- mm2[meta$group == "tumor_id",] |> colMeans()
tumor_iec <- mm2[meta$group == "tumor_iec",] |> colMeans()
tumor_iic <- mm2[meta$group == "tumor_iic",] |> colMeans()

contrast_list <- list(
  "tumor ieb / iib" = tumor_ieb - tumor_iib, 
  "tumor idb / iib" = tumor_idb - tumor_iib, 
  "tumor ieb / idb" = tumor_ieb - tumor_idb, 
  "tumor iec / iic" = tumor_iec - tumor_iic, 
  "tumor ieb / iec" = tumor_ieb - tumor_iec, 
  "tumor iib / iic" = tumor_iib - tumor_iic
)

vn <- voom(counts = cts, design = mm2)
limmafit <- lmFit(object = vn, design = mm2)

pdf(file = "comparisons_second_cohort_volcanos_limma.pdf", width = 8, height = 8)
for (c in names(contrast_list)) {
  deout <- contrasts.fit(limmafit, contrasts = contrast_list[[c]])
  deout <- eBayes(deout) |> topTable(number = Inf)
  deout$target <- rownames(deout)
  deout$contrast <- c
  
  write.csv(deout, paste(gsub(pattern = "/", replacement = "vs", x = c), " limma.csv", sep = ""))
  
  deout <- rename(deout, fdr = adj.P.Val, coefficient = logFC, pvalue = P.Value)
  print(plot_volcano(deout))
}
dev.off()

run_cols <- c("1"="grey", "2"="black")
pt_cols <- ggprism::ggprism_data$colour_palettes$colors[1:8]
names(pt_cols) <- 1:8

pdf(file = "comparisons_second_cohort_heatmaps_limma.pdf", width = 25, height = 50)
meta$group <- meta$infiltration_type
for (c in names(contrast_list)) {
  restest <- read.csv(paste(gsub(pattern = "/", replacement = "vs", x = c), " limma.csv", sep = ""), row.names = 1)
  restest <- rename(restest, fdr = adj.P.Val, coefficient = logFC, pvalue = P.Value)
  restest <- restest |> filter(fdr < 0.05)
  if (dim(restest)[1] > 0) {
    cs <- str_split(pattern = " / |\\+", simplify = T, string = c)[1,] |> gsub(pattern = "tumor ", replacement = "")
    infil_cols <- c("idb" = "red", "iib" = "blue", "ieb" = "green", "iec"="orange", "iic"="magenta")[cs]
    idx <- meta[(meta$group %in% cs) & (meta$AOI_code == "tumor"),] |> arrange(group) |> pull(aoi_id)
    toplot <- restest |> arrange(desc(coefficient)) %$% target
    
    ha <- HeatmapAnnotation(group = meta[idx,]$group, 
                            pt = meta[idx,]$Patient_number,
                            run = meta[idx,]$Run_number,
                            col = list("group"= infil_cols, 
                                       "pt" = pt_cols, 
                                       "run" = run_cols))
    mat <- vn[toplot, idx, drop = F] |> apply(1, scale) |> t()
    hm <- Heatmap(matrix = mat, 
                  top_annotation = ha, 
                  name = "Scaled\nlog2(CPM + 1)", 
                  cluster_columns = F, cluster_rows = F,
                  width = ncol(mat)*unit(3, "mm"), 
                  height = nrow(mat)*unit(1.5, "mm"),
                  rect_gp = gpar(col = "white", lwd = 0.1), 
                  row_names_gp = gpar(fontsize = 3),
                  column_title = c
    )
    draw(hm)
  }
}

dev.off()


# This all took about 6 hours

