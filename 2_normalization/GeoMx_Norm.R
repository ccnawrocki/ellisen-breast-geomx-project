###### GeoMx Data Normalization ######
## Cole Nawrocki ##

# We want to determine if normalizing via Q3, with limma's vst (log2CPM), was best. 
# We also want to explore limma's batch correction function for data like this. 

# All writing and saving functions are commented out to avoid overwriting.

pdf("Normalization_plots.pdf", width = 12, height = 8)

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

# Colors
segment_cols <- c("tumor" = "green4", "immune" = "red", "fibroblast" = "dodgerblue")
run_cols <- c("R1"="yellow", "R2"="orange")
slide_cols <- grDevices::colors()[grep('(gr(a|e)y)|(white)', grDevices::colors(), invert = T)] |> sample(size = 8, replace = F)
names(slide_cols) <- unique(meta$slide.name)

# Counts and metadata
cts <- read.csv("counts.csv", row.names = 1)
meta <- read.csv("meta.csv", row.names = 1)
all(rownames(meta) == colnames(cts))

### Q3 Normalization ### -------------------------------------------------------
# This is equivalent to the Q3 normalization in the GeomxTools package. 
qs <- apply(X = cts, MARGIN = 2, FUN = quantile, 0.75)
nfs <- qs / EnvStats::geoMean(qs)
q3norm <- sweep(x = cts, MARGIN = 2, STATS = nfs, FUN = "/")

# I will be working with the log-transformation of the Q3 normalized data.
# log2(q3-normalized count + 1)
lq3 <- log2(q3norm + 1)

# Identifying the 5% of genes with the highest variance
variances <- apply(X = lq3, MARGIN = 1, FUN = var)
top_genes <- names(variances)[variances >= quantile(variances, 0.95)]

# PCA
# There is a patient-driven batch effect.
mat <- lq3[top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = meta$segment) + 
  theme_classic() 
ggbiplot(pcaout, var.axes = F, groups = meta$slide.name) + 
  theme_classic()
ggbiplot(pcaout, var.axes = F, groups = meta$Run) + 
  theme_classic()

# Heatmap with these genes
# There is a huge patient-driven batch effect in the tumor. This makes sense. 
# The batch effect is not as pronounced in the other segments, but it is still 
# present. Q3 seems to work well, but the batch effect will need to be accounted
# for in downstream analysis. 
ha <- HeatmapAnnotation(segment = meta$segment, 
                        slide =  meta$slide.name, 
                        run = meta$Run,
                        col = list("segment"=segment_cols, 
                                   "slide"=slide_cols, 
                                   "run"=run_cols))
Heatmap(matrix = mat |> t(), 
        top_annotation = ha,
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled\nlog2(Q3+1)", 
        show_column_names = F, 
        show_row_names = F
)

# Fitting a linear mixed effects model to the data, with segment as the fixed effect and slide as the random effect.
# Then, determining the 5% of genes that have the smallest MSE values (based on Pearson residuals). 
meta$segment %<>% as.factor()
meta$slide.name %<>% as.factor()
lmem_pear_mse <- function(gene) {
  d <- dplyr::mutate(meta, expr=(lq3[gene,] |> unlist()))
  testout <- lme4::lmer(formula = expr ~ segment + (1+segment|slide.name), data = d)
  return((resid(testout, type = "pear")**2 |> sum()) / nrow(d))
}

mses <- parallel::mclapply(rownames(lq3), lmem_pear_mse, mc.cores = 8)
names(mses) <- rownames(lq3)
mses %<>% unlist()

top_genes <- names(mses)[mses <= quantile(mses, 0.05)]

# PCA
# Batch effect is still visible, but a bit less pronounced. This makes sense. 
# The pca is not taking the batch effect into account like the model is. But, 
# we have still selected a "less-affected" subset of genes in a sense. So, we 
# may expect to see marginally-better results. 
mat <- lq3[top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = meta$segment) + 
  theme_classic() 
ggbiplot(pcaout, var.axes = F, groups = meta$slide.name) + 
  theme_classic()
ggbiplot(pcaout, var.axes = F, groups = meta$Run) + 
  theme_classic()

# Heatmap with these genes
# The logic above carries over to the heatmap. 
Heatmap(matrix = mat |> t(), 
        top_annotation = ha, 
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled\nlog2(Q3+1)", 
        show_column_names = F, 
        show_row_names = F
)

# Saving
# I think that we have established that Q3 normalization will be used downstream. 
# Results look decent, and it is the recommended method, so we will save it. 
#write.csv(x = lq3, file = "q3norm.csv")

### Normalization with limma ### -----------------------------------------------
# Creating the object
cts <- DGEList(counts = cts, samples = meta)

# Calculating normalization factors
cts <- calcNormFactors(object = cts)

# Normalizing with voom
norm <- voom(counts = cts, plot = T)

# Fitting the limma model
limfit <- lmFit(norm)

# Identifing the 5% most "variable" genes via dispersion
disps <- cbind(limfit$coefficients, "sigma"=limfit$sigma) |> as.data.frame() 
disp_top_genes <- rownames(disps)[disps$sigma >= quantile(disps$sigma, 0.95)]

# PCA
# Does not look any better than the Q3 results
mat <- norm$E[disp_top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = meta$slide.name) + 
  theme_classic()
ggbiplot(pcaout, var.axes = F, groups = meta$segment) + 
  theme_classic()
ggbiplot(pcaout, var.axes = F, groups = meta$Run) + 
  theme_classic()

# Heatmap with these "variable" genes
# Similar interpretation to the Q3 results above
Heatmap(matrix = mat |> t(), 
        top_annotation = ha, 
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled\nlog2CPM", 
        show_row_names = F
)

# Refitting the limma model with a similar "mixed" model setup as above 
# It is not truly mixed; it takes duplicate correlation into account.
mm <- model.matrix(~segment, data = meta)
cts <- DGEList(counts = cts)
cts <- calcNormFactors(object = cts)

# The voom trend looks marginally better
norm <- voom(cts, design = mm, plot = T)
corfit <- duplicateCorrelation(norm, design = mm, block = meta$slide.name)
limfit <- lmFit(object = norm, design = mm, block = meta$slide.name, correlation = corfit$consensus)

# Identifing the 5% most "variable" genes via dispersion 
disps <- cbind(limfit$coefficients, "sigma"=limfit$sigma) |> as.data.frame() 
disp_top_genes <- rownames(disps)[disps$sigma >= quantile(disps$sigma, 0.95)]

# PCA
# Tough to say if there is any improvement
mat <- norm$E[disp_top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = meta$slide.name) + 
  theme_classic()
ggbiplot(pcaout, var.axes = F, groups = meta$segment) + 
  theme_classic()
ggbiplot(pcaout, var.axes = F, groups = meta$Run) + 
  theme_classic()

# Heatmap with these "variable" genes
# Honestly, it looks worse.
Heatmap(matrix = mat |> t(), 
        top_annotation = ha, 
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled log2CPM", 
        show_row_names = F
)

### Batch Correction ### -------------------------------------------------------
# We will use limma for batch correction of the q3 normalized data. 
# We will pass it patient as the batch. I originally passed some over info to the
# function, but I think that it is best to only be using the function for the 
# one intended purpose: mitigating the patient's influence in the data. 
norm_bc <- limma::removeBatchEffect(x = lq3, 
                                    batch = meta$slide.name, 
                                    #group = meta$segment, 
                                    #covariates = meta$area
)

# We will find variable genes with variance.
variances <- apply(X = norm_bc, MARGIN = 1, FUN = var)
top_genes <- names(variances)[variances >= quantile(variances, 0.95)]

# PCA
mat <- norm_bc[top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = meta$slide.name) + 
  theme_classic() 
ggbiplot(pcaout, var.axes = F, groups = meta$segment) + 
  theme_classic()
ggbiplot(pcaout, var.axes = F, groups = meta$Run) + 
  theme_classic() 

# Heatmap with the "variable" genes
# It looks nice. 
Heatmap(matrix = mat |> t(), 
        top_annotation = ha, 
        row_names_gp = gpar(fontsize = 3), 
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled\nBatch\nCorrected\nlog2(Q3+1)", 
        show_row_names = F
        )

# What if we used the genes found from the true mixed model above? 
# Instead of using MSE, we will use the global F statistic to rank genes based 
# on their model significance when segment is a predictor. Then, we will take 
# the top 5%.
# lmem_fstat <- function(gene) {
#   d <- dplyr::mutate(meta, expr=(lq3[gene,] |> unlist()))
#   testout <- lme4::lmer(formula = expr ~ segment + (1+segment|slide.name), data = d) |> 
#     car::Anova(test.statistic = "F")
#   rownames(testout) <- gene
#   return(testout)
# }
# ftests <- parallel::mclapply(rownames(lq3), lmem_fstat, mc.cores = 10)
# ftests <- bind_rows(ftests, .id = "target")

# Took a while to run, so I am going to save this. Then I can just read it back 
# in whenever.
#write.csv(ftests, "2_normalization/ftest_results.csv")
ftests <- read.csv(file = "2_normalization/ftest_results.csv", row.names = 1)
top_genes <- rownames(ftests)[ftests$`F` >= quantile(ftests$`F`, 0.95)]

# It looks great. This is the way to go: pick genes via modeling, then use the
# batch-corrected data for visualization/clustering. This is the go-to in 
# scRNA-seq.
mat <- norm_bc[top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = meta$slide.name) + 
  theme_classic() 
ggbiplot(pcaout, var.axes = F, groups = meta$segment) + 
  theme_classic() 
ggbiplot(pcaout, var.axes = F, groups = meta$Run) + 
  theme_classic() 

# Totally unsupervised
Heatmap(matrix = mat |> t(), 
        top_annotation = ha, 
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled\nBatch\nCorrected\nlog2(Q3+1)", 
        show_row_names = F
)

# Unsupervised, but split into 3 groups
Heatmap(matrix = mat |> t(), 
        top_annotation = ha, 
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled\nBatch\nCorrected\nlog2(Q3+1)", 
        show_row_names = F, 
        column_split = 3
)

# Supervised
Heatmap(matrix = mat |> t(), 
        top_annotation = ha, 
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled\nBatch\nCorrected\nlog2(Q3+1)", 
        show_row_names = F, 
        column_split = meta$segment
)

# Saving
# I beleive that the batch-corrected data has value for visualization purposes
# only. So, I will save it.
#write.csv(x = norm_bc, file = "q3norm_bc.csv")

### Comparison ### -------------------------------------------------------------
boxplot(norm_bc,
        col = "#9EDAE5", 
        main = "Batch-corrected log2(Q3+1) Normalized Counts", 
        names = 1:155, xlab = "AOI #")
boxplot(norm$E,
        col = "green4", 
        main = "Voom Normalized Counts log2(CPM+0.5)", 
        names = 1:155, xlab = "AOI #")
boxplot(lq3,
        col = "orange", 
        main = "log2(Q3+1) Normalized Counts", 
        names = 1:155, xlab = "AOI #")
boxplot(q3norm,
        col = "purple", 
        main = "Q3 Normalized Counts", 
        names = 1:155, xlab = "AOI #")
boxplot(cts$counts,
        col = "red3", 
        main = "Raw Counts", 
        names = 1:155, xlab = "AOI #")

### Thoughts/Conclusions ### ---------------------------------------------------
# The Q3 normalized counts are what we should use. The PCA and clustering based 
# on the batch-corrected normalized counts will be useful for visualization 
# only. 

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
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/New_York
# tzcode source: system (macOS)
# 
# attached base packages:
#  [1] stats4    parallel  grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] cowplot_1.1.3            reshape2_1.4.4           knitr_1.48               openxlsx_4.2.7.1         data.table_1.16.0       
#  [6] GeoMxWorkflows_1.10.0    GeomxTools_3.8.0         NanoStringNCTools_1.12.0 S4Vectors_0.42.1         Biobase_2.64.0          
# [11] BiocGenerics_0.50.0      ggh4x_0.3.0              ggdendro_0.2.0           ggbiplot_0.6.2           edgeR_4.2.1             
# [16] limma_3.60.4             Seurat_5.2.0             SeuratObject_5.0.2       sp_2.1-4                 SpatialDecon_1.14.0     
# [21] multcomp_1.4-26          TH.data_1.1-2            MASS_7.3-61              survival_3.7-0           mvtnorm_1.3-1           
# [26] lmerTest_3.1-3           lme4_1.1-35.5            Matrix_1.7-0             patchwork_1.3.0          magrittr_2.0.3          
# [31] ComplexHeatmap_2.20.0    ggplot2_3.5.1            stringr_1.5.1            dplyr_1.1.4             
# 
# loaded via a namespace (and not attached):
#   [1] matrixStats_1.4.1       spatstat.sparse_3.1-0   httr_1.4.7              RColorBrewer_1.1-3      doParallel_1.0.17      
#   [6] numDeriv_2016.8-1.1     backports_1.5.0         tools_4.4.1             sctransform_0.4.1       utf8_1.2.4             
#  [11] R6_2.5.1                lazyeval_0.2.2          uwot_0.2.2              GetoptLong_1.0.5        withr_3.0.1            
#  [16] GGally_2.2.1            gridExtra_2.3           progressr_0.14.0        cli_3.6.3               spatstat.explore_3.3-4 
#  [21] fastDummies_1.7.5       sandwich_3.1-1          labeling_0.4.3          spatstat.data_3.1-4     ggridges_0.5.6         
#  [26] pbapply_1.7-2           askpass_1.2.0           systemfonts_1.1.0       R.utils_2.12.3          parallelly_1.38.0      
#  [31] readxl_1.4.3            rstudioapi_0.16.0       generics_0.1.3          shape_1.4.6.1           ica_1.0-3              
#  [36] spatstat.random_3.3-2   car_3.1-3               zip_2.3.1               ggbeeswarm_0.7.2        fansi_1.0.6            
#  [41] abind_1.4-8             R.methodsS3_1.8.2       lifecycle_1.0.4         yaml_2.3.10             carData_3.0-5          
#  [46] Rtsne_0.17              promises_1.3.0          crayon_1.5.3            miniUI_0.1.1.1          lattice_0.22-6         
#  [51] magick_2.8.5            pillar_1.9.0            rjson_0.2.23            boot_1.3-31             future.apply_1.11.2    
#  [56] codetools_0.2-20        glue_1.8.0              ggiraph_0.8.10          outliers_0.15           spatstat.univar_3.1-1  
#  [61] vctrs_0.6.5             png_0.1-8               spam_2.10-0             cellranger_1.1.0        gtable_0.3.5           
#  [66] xfun_0.47               mime_0.12               pheatmap_1.0.12         iterators_1.0.14        statmod_1.5.0          
#  [71] fitdistrplus_1.2-2      ROCR_1.0-11             nlme_3.1-166            pbkrtest_0.5.3          EnvStats_3.0.0         
#  [76] RcppAnnoy_0.0.22        GenomeInfoDb_1.40.1     R.cache_0.16.0          irlba_2.3.5.1           vipor_0.4.7            
#  [81] KernSmooth_2.23-24      colorspace_2.1-1        tidyselect_1.2.1        logNormReg_0.5-0        repmis_0.5             
#  [86] compiler_4.4.1          plotly_4.10.4           scales_1.3.0            lmtest_0.9-40           digest_0.6.37          
#  [91] goftest_1.2-3           spatstat.utils_3.1-2    minqa_1.2.8             rmarkdown_2.28          XVector_0.44.0         
#  [96] htmltools_0.5.8.1       pkgconfig_2.0.3         umap_0.2.10.0           fastmap_1.2.0           rlang_1.1.4            
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

