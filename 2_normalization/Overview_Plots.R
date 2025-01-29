# Data Overview Plots 

.libPaths()
# [1] "/Users/ccn22/micromamba/envs/geomx-env/lib/R/library"

# Packages
library(dplyr)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(ggbiplot)
library(magrittr)
library(patchwork)

# Data 
meta <- readxl::read_excel("meta_cleaned.xlsx", sheet = 1)
rownames(meta) <- meta$aoi_id
norm <- read.csv("q3norm.csv", row.names = 1)

pdf(file = "Overview_Plots.pdf", width = 10, height = 8)

# Colors
segment_cols <- c("tumor" = "green4", "immune" = "red", "fibroblast" = "dodgerblue")
run_cols <- c("1"="yellow", "2"="orange")
slide_cols <- grDevices::colors()[grep('(gr(a|e)y)|(white)', grDevices::colors(), invert = T)] |> sample(size = 8, replace = F)
names(slide_cols) <- unique(meta$Patient_number)

# PCA
variances <- apply(X = lq3, MARGIN = 1, FUN = var)
top_genes <- names(variances)[variances >= quantile(variances, 0.95)]

# PCA
mat <- norm[top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
meta <- cbind(meta, pcaout$x[,1:2])
ggplot(meta) + 
  geom_text(mapping = aes(x = PC1, y = PC2, color = AOI_code, label = Patient_number)) + 
  labs(title = "PCA: Labeled by Patient Number") + 
  scale_color_manual(values = segment_cols) + 
  theme_classic()
ggplot(meta) + 
  geom_text(mapping = aes(x = PC1, y = PC2, color = AOI_code, label = Run_number)) + 
  labs(title = "PCA: Labeled by Run Number") +
  scale_color_manual(values = segment_cols) + 
  theme_classic()

# Heatmap
ha <- HeatmapAnnotation(segment = meta$AOI_code, 
                        pt =  meta$Patient_number, 
                        run = meta$Run_number,
                        col = list("segment"=segment_cols, 
                                   "pt"=slide_cols, 
                                   "run"=run_cols))
Heatmap(matrix = mat |> t(), 
        top_annotation = ha, 
        row_names_gp = gpar(fontsize = 3), 
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled\nQ3-Normalized\nCounts", 
        show_column_names = F, 
        show_row_names = F
)

boxplot(norm,
        main = "Q3 Normalized Counts", 
        names = 1:155, xlab = "AOI #")

dev.off()

