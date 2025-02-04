###### GeoMx Data Segment Markers ######
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

# Data: we will use the q3 normalized expression for modeling.
meta <- read.csv("meta.csv", row.names = 1)
norm <- read.csv("q3norm.csv", row.names = 1)
norm <- norm[rownames(norm) != "NegProbe-WTX",]
bcnorm <- read.csv("q3norm_bc.csv", row.names = 1)

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

# Testing each segment's contrast with a an Exact Wald Test (t test with Satterthwaite df approximation).
# This is a somewhat strict test.
# lmem_de <- function(gene) {
#   d <- dplyr::mutate(meta, expr=(norm[gene,] |> unlist()))
#   modelout <- lmerTest::lmer(formula = expr ~ segment + (1+segment|slide.name), data = d)
#   out <- list()
#   for (seg in names(contrast_list)) {
#     seg_oi <- contrast_list[[seg]]
#     other <- Reduce(f = "+",
#                     x = contrast_list[names(contrast_list) != seg]) /
#       (length(contrast_list[names(contrast_list) != seg]))
#     seg_oi_vs_other <- seg_oi - other
#     out[[seg]] <- contest(model = modelout, L = seg_oi_vs_other, joint = F, rhs = 0)
#   }
#   out <- bind_rows(out, .id = "group")
#   out$gene <- gene
#   return(out)
# }
# 
# de <- parallel::mclapply(rownames(norm), lmem_de, mc.cores = 10)
# dedf <- bind_rows(de)
# dedf <- group_by(dedf, group) |> mutate(fdr = p.adjust(method = "fdr", `Pr(>|t|)`))
# write.csv(dedf, file = "3_segment_markers/segment_markers_results.csv")
dedf <- read.csv(file = "3_segment_markers/segment_markers_results.csv", row.names = 1)

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
        cluster_columns = T, 
        name = "Scaled\nlog2(Q3+1)", 
        column_split = meta[idx, ]$segment
)

# Heatmap, with batch-corrected data
idx <- meta |> arrange(segment) |> rownames()
mat <- bcnorm[toplot$gene, idx] |> apply(1, scale) |> t()
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
        cluster_columns = T, 
        name = "Scaled\nBatch\nCorrected\nlog2(Q3+1)", 
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

sessionInfo()

# R version 4.4.1 (2024-06-14)
# Platform: aarch64-apple-darwin20.0.0
# Running under: macOS 15.1.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Users/ccn22/micromamba/envs/geomx-env/lib/libopenblas.0.dylib;  LAPACK version 3.12.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/New_York
# tzcode source: system (macOS)
# 
# attached base packages:
#   [1] stats4    parallel  grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] cowplot_1.1.3            reshape2_1.4.4           knitr_1.48               openxlsx_4.2.7.1         data.table_1.16.0       
# [6] GeoMxWorkflows_1.10.0    GeomxTools_3.8.0         NanoStringNCTools_1.12.0 S4Vectors_0.42.1         Biobase_2.64.0          
# [11] BiocGenerics_0.50.0      ggh4x_0.3.0              ggdendro_0.2.0           ggbiplot_0.6.2           edgeR_4.2.1             
# [16] limma_3.60.4             Seurat_5.2.0             SeuratObject_5.0.2       sp_2.1-4                 SpatialDecon_1.14.0     
# [21] multcomp_1.4-26          TH.data_1.1-2            MASS_7.3-61              survival_3.7-0           mvtnorm_1.3-1           
# [26] lmerTest_3.1-3           lme4_1.1-35.5            Matrix_1.7-0             patchwork_1.3.0          magrittr_2.0.3          
# [31] ComplexHeatmap_2.20.0    ggplot2_3.5.1            stringr_1.5.1            dplyr_1.1.4             
# 
# loaded via a namespace (and not attached):
#   [1] matrixStats_1.4.1       spatstat.sparse_3.1-0   httr_1.4.7              RColorBrewer_1.1-3      doParallel_1.0.17      
# [6] numDeriv_2016.8-1.1     backports_1.5.0         tools_4.4.1             sctransform_0.4.1       utf8_1.2.4             
# [11] R6_2.5.1                lazyeval_0.2.2          uwot_0.2.2              GetoptLong_1.0.5        withr_3.0.1            
# [16] GGally_2.2.1            gridExtra_2.3           progressr_0.14.0        cli_3.6.3               spatstat.explore_3.3-4 
# [21] fastDummies_1.7.5       sandwich_3.1-1          labeling_0.4.3          spatstat.data_3.1-4     ggridges_0.5.6         
# [26] pbapply_1.7-2           askpass_1.2.0           systemfonts_1.1.0       R.utils_2.12.3          parallelly_1.38.0      
# [31] readxl_1.4.3            rstudioapi_0.16.0       generics_0.1.3          shape_1.4.6.1           ica_1.0-3              
# [36] spatstat.random_3.3-2   car_3.1-3               zip_2.3.1               ggbeeswarm_0.7.2        fansi_1.0.6            
# [41] abind_1.4-8             R.methodsS3_1.8.2       lifecycle_1.0.4         yaml_2.3.10             carData_3.0-5          
# [46] Rtsne_0.17              promises_1.3.0          crayon_1.5.3            miniUI_0.1.1.1          lattice_0.22-6         
# [51] magick_2.8.5            pillar_1.9.0            rjson_0.2.23            boot_1.3-31             future.apply_1.11.2    
# [56] codetools_0.2-20        glue_1.8.0              ggiraph_0.8.10          outliers_0.15           spatstat.univar_3.1-1  
# [61] vctrs_0.6.5             png_0.1-8               spam_2.10-0             cellranger_1.1.0        gtable_0.3.5           
# [66] xfun_0.47               mime_0.12               pheatmap_1.0.12         iterators_1.0.14        statmod_1.5.0          
# [71] fitdistrplus_1.2-2      ROCR_1.0-11             nlme_3.1-166            pbkrtest_0.5.3          EnvStats_3.0.0         
# [76] RcppAnnoy_0.0.22        GenomeInfoDb_1.40.1     R.cache_0.16.0          irlba_2.3.5.1           vipor_0.4.7            
# [81] KernSmooth_2.23-24      colorspace_2.1-1        tidyselect_1.2.1        logNormReg_0.5-0        repmis_0.5             
# [86] compiler_4.4.1          plotly_4.10.4           scales_1.3.0            lmtest_0.9-40           digest_0.6.37          
# [91] goftest_1.2-3           spatstat.utils_3.1-2    minqa_1.2.8             rmarkdown_2.28          XVector_0.44.0         
# [96] htmltools_0.5.8.1       pkgconfig_2.0.3         umap_0.2.10.0           fastmap_1.2.0           rlang_1.1.4            
# [101] GlobalOptions_0.1.2     htmlwidgets_1.6.4       ggthemes_5.1.0          UCSC.utils_1.0.0        shiny_1.9.1            
# [106] farver_2.1.2            zoo_1.8-12              jsonlite_1.8.9          R.oo_1.26.0             Formula_1.2-5          
# [111] GenomeInfoDbData_1.2.12 dotCall64_1.1-1         munsell_0.5.1           Rcpp_1.0.13             reticulate_1.39.0      
# [116] stringi_1.8.4           zlibbioc_1.50.0         plyr_1.8.9              ggstats_0.7.0           listenv_0.9.1          
# [121] ggrepel_0.9.6           deldir_2.0-4            Biostrings_2.72.1       splines_4.4.1           tensor_1.5             
# [126] circlize_0.4.16         locfit_1.5-9.6          igraph_2.0.3            uuid_1.2-1              spatstat.geom_3.3-5    
# [131] RcppHNSW_0.6.0          evaluate_1.0.0          BiocManager_1.30.25     ggprism_1.0.5           nloptr_2.1.1           
# [136] foreach_1.5.2           tweenr_2.0.3            httpuv_1.6.15           networkD3_0.4           RANN_2.6.2             
# [141] tidyr_1.3.1             openssl_2.2.2           purrr_1.0.2             polyclip_1.10-7         future_1.34.0          
# [146] clue_0.3-65             scattermore_1.2         ggforce_0.4.2           broom_1.0.7             xtable_1.8-4           
# [151] RSpectra_0.16-2         later_1.3.2             viridisLite_0.4.2       tibble_3.2.1            beeswarm_0.4.0         
# [156] IRanges_2.38.1          cluster_2.1.6           globals_0.16.3          BiocStyle_2.32.1       

