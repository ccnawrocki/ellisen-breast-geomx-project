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

# limma
meta$group <- paste(meta$AOI_code, meta$infiltration_category, sep = "_")
mm <- model.matrix(~1+group+Patient_number, data = meta)
fibroblast_id <- mm[meta$group == "fibroblast_desert",] |> colMeans()
fibroblast_ii <- mm[meta$group == "fibroblast_inflamed",] |> colMeans()
fibroblast_ie <- mm[meta$group == "fibroblast_excluded",] |> colMeans()

contrast_list <- list(
  "fibroblast id / ii+ie" = fibroblast_id - (fibroblast_ii+fibroblast_ie)/2
)

cts <- DGEList(counts = counts)
cts <- calcNormFactors(cts)
vn <- voom(counts = cts, design = mm, plot = T)
limmafit <- lmFit(object = vn, design = mm)

pdf(file = "volcano_limma_for_trop2_paper_bullet3.pdf", width = 8, height = 8)
for (c in names(contrast_list)) {
  deout <- contrasts.fit(limmafit, contrasts = contrast_list[[c]])
  deout <- eBayes(deout) |> topTable(number = Inf)
  deout$target <- rownames(deout)
  deout$contrast <- c
  
  #write.csv(deout, paste(gsub(pattern = "/", replacement = "vs", x = c), " limma_for_trop2_paper_bullet3.csv", sep = ""))
  
  deout <- rename(deout, fdr = adj.P.Val, coefficient = logFC, pvalue = P.Value)
  print(plot_volcano(deout))
}

dev.off()


