library(ggprism)

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

### FIBROBLAST DECONVOLUTION PLOT
fibdcvnres <- read.csv("6_fibroblast_deconvolution/fibroblast_deconvolution_differential_abundance_testing_results_for_trop2_paper.csv", row.names = 1)
fibdcvnres |> filter(celltype != "hsp_tpCAF") |> 
  ggplot() +
  geom_bar(mapping = aes(y = celltype, x = mean_diff, fill = padj < 0.05), stat = "identity") + 
  facet_grid(.~contrast) +
  scale_fill_manual(values = c("blue", "red")) + 
  ggthemes::theme_few() + 
  geom_vline(xintercept = 0, linewidth = 0.2)

dev.off()
