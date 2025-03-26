# These are just a couple specific plots that Bogang wanted for the TROP2 paper.

library(ggprism)
library(ggplot2)

pdf("ellisen_plots.pdf", width = 6, height = 8)

### TROP2 EXPRESSION PLOT
meta <- openxlsx::read.xlsx("meta_cleaned_v5.xlsx")
rownames(meta) <- meta$aoi_id
meta$group <- paste(meta$AOI_code, meta$infiltration_category, sep = "_")
norm <- read.csv("q3norm.csv", row.names = 1)
norm <- t(norm) |> as.data.frame()
norm$aoi_id <- rownames(norm)

trop2 <- inner_join(norm[,c("aoi_id", "TACSTD2")], meta[,c("aoi_id", "AOI_code", "infiltration_category")], by = "aoi_id")
tumor_trop2 <- trop2 |> filter(AOI_code == "tumor") |> rename(TROP2 = TACSTD2)
tumor_trop2$hot_cold_group <- ifelse(test = tumor_trop2$infiltration_category == "desert", yes = "cold", no = "hot")
fibroblast_trop2 <- trop2 |> filter(AOI_code == "fibroblast") |> rename(TROP2 = TACSTD2)
fibroblast_trop2$hot_cold_group <- ifelse(test = fibroblast_trop2$infiltration_category == "desert", yes = "cold", no = "hot")

df_p_val <- data.frame(
  group1 = "cold",
  group2 = "hot",
  label = "*",
  y.position = 7.5
)

ggplot(data = tumor_trop2, mapping = aes(x = hot_cold_group, y = TROP2)) + 
  geom_boxplot(outliers = F, aes(color = hot_cold_group)) +
  geom_jitter(width = 0.1, height = 0) + 
  scale_color_manual(values = c("dodgerblue", "red2")) + 
  ggprism::add_pvalue(data = df_p_val, label.size = 5) +
  guides(y = "prism_offset_minor") + 
  ggprism::theme_prism() + 
  labs(title = "Tumor TROP2") +
  theme(axis.title = element_blank(), title = element_text(hjust = 0.5))

ggplot(data = tumor_trop2, mapping = aes(x = hot_cold_group, y = TROP2)) + 
  geom_jitter(width = 0.1, height = 0) + 
  stat_summary(geom = "crossbar", fun = "mean", aes(color = hot_cold_group)) +
  scale_color_manual(values = c("dodgerblue", "red2")) +
  ggprism::theme_prism() + 
  guides(y = "prism_offset_minor") + 
  ggprism::add_pvalue(data = df_p_val, label.size = 5) +
  labs(title = "Tumor TROP2") +
  theme(axis.title = element_blank(), title = element_text(hjust = 0.5))

df_p_val <- data.frame(
  group1 = "cold",
  group2 = "hot",
  label = "ns",
  y.position = 6
)

ggplot(data = fibroblast_trop2, mapping = aes(x = hot_cold_group, y = TROP2)) + 
  geom_boxplot(outliers = F, aes(color = hot_cold_group)) +
  geom_jitter(width = 0.1, height = 0) + 
  scale_color_manual(values = c("dodgerblue", "red2")) +
  ggprism::theme_prism() + 
  guides(y = "prism_offset_minor") + 
  ggprism::add_pvalue(data = df_p_val, label.size = 5) +
  labs(title = "Fibroblast TROP2") +
  theme(axis.title = element_blank(), title = element_text(hjust = 0.5))

ggplot(data = fibroblast_trop2, mapping = aes(x = hot_cold_group, y = TROP2)) + 
  geom_jitter(width = 0.1, height = 0) + 
  stat_summary(geom = "crossbar", fun = "mean", aes(color = hot_cold_group)) +
  scale_color_manual(values = c("dodgerblue", "red2")) +
  ggprism::theme_prism() + 
  guides(y = "prism_offset_minor") + 
  ggprism::add_pvalue(data = df_p_val, label.size = 5) +
  labs(title = "Fibroblast TROP2") +
  theme(axis.title = element_blank(), title = element_text(hjust = 0.5))

#write.csv(tumor_trop2, "tumor_trop2.csv")
#write.csv(fibroblast_trop2, "fibroblast_trop2.csv")

dev.off()

### FIBROBLAST DECONVOLUTION PLOT
fibdcvnres <- read.csv("6_fibroblast_deconvolution/fibroblast_deconvolution_differential_abundance_testing_results_for_trop2_paper.csv", row.names = 1)
fibdcvnres$group <- ifelse(fibdcvnres$mean_diff < 0, yes = "hot", no = "cold")
fibdcvnres$`padj<0.05` <- ifelse(test = fibdcvnres$padj < 0.05, yes = "*", no = "")
fibdcvnres |> filter(celltype != "hsp_tpCAF") |> 
  ggplot() +
  geom_bar(mapping = aes(y = celltype, x = mean_diff, fill = group), stat = "identity") + 
  geom_text(mapping = aes(y = celltype, x = (mean_diff)+0.025*sign(mean_diff), label = `padj<0.05`), size = 6, vjust = 0.75) +
  scale_fill_manual(values = c("dodgerblue", "red2")) + 
  ggthemes::theme_few() + 
  geom_vline(xintercept = 0, linewidth = 0.2) + 
  labs(title = "Fibroblast deconvolution", subtitle = "cold vs. hot", x = "mean difference (cold - hot)", y = "") + 
  theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#ggsave(filename = "updated_S1B.pdf", device = "pdf", width = 6, height = 6)

### TUMOR VOLCANO PLOT
tumderes <- read.csv("4_DE_analysis/comparisons_first_cohort/tumor id vs ii+ie limma.csv", row.names = 1)
tumderes$significance <- "not significant"
tumderes[tumderes$logFC > 0 & tumderes$adj.P.Val < 0.05,]$significance <- "Padj<0.05, up in cold"
tumderes[tumderes$logFC < 0 & tumderes$adj.P.Val < 0.05,]$significance <- "Padj<0.05, up in hot"
ggplot() + 
  geom_point(data = tumderes, mapping = aes(x = logFC, y = -log10(adj.P.Val), color = significance)) +
  ggrepel::geom_text_repel(data = dplyr::filter(tumderes, target %in% c("MGP", "VIM", "HLA-DRB1", "B2M", "TAP1", "TACSTD2", "CLDN1", "CLDN7") & logFC < 0), 
                           mapping = aes(x = logFC, y = -log10(adj.P.Val), label = target), size = 3,
                           max.overlaps = 5, ylim = c(2, 10), xlim = c(-2, -0.5), force_pull = 0, min.segment.length = 0, position = ggrepel::position_nudge_repel(y = 1), box.padding = 2) +
  ggrepel::geom_text_repel(data = dplyr::filter(tumderes, target %in% c("MGP", "VIM", "HLA-DRB1", "B2M", "TAP1", "TACSTD2", "CLDN1", "CLDN7") & logFC > 0), 
                           mapping = aes(x = logFC, y = -log10(adj.P.Val), label = target), size = 3,
                           max.overlaps = 8, ylim = c(2, 20), xlim = c(0, 3), force_pull = 0, min.segment.length = 0, position = ggrepel::position_nudge_repel(y = 1), box.padding = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 1) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey", linewidth = 1) + 
  theme_prism() + theme(legend.position = c(0.15, 0.85), legend.text = element_text(size = 10)) + guides(y = "prism_offset_minor") + 
  scale_color_manual(values = c("black", "dodgerblue", "red3")) +
  labs(x = "log2FC (cold/hot)", y = "-log10(Padj)", title = "Differentially Expressed\nGenes for PanCK+ Cells")
