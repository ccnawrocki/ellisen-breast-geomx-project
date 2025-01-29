##### Ellisen GeoMx Immune Deconvolution #####
### Cole Nawrocki ###

set.seed(2001)
pdf("immune_deconvolution.pdf", width = 12, height = 8)

## PLAN ------------------------------------------------------------------------
# Using the T cell atlas from this paper: https://doi.org/10.1038/s41591-023-02371-y
# to do deconvolution on the immune segments only. 
# The segments are not CD45 positive. Rather, they are CD8 positive. So, we cannot
# just fo immune deconvolution; we need to do CD8 T cell deconvolution. Bogang has
# said that this paper has reliable scRNA-seq data with CD8 T cell subtypes types 
# that he believes are implicated in TNBC tumor infiltration.

# First, I will read the paper, download the data, and make a custom reference profile, 
# using the SpatialDecon package. 

# Second, I will perform CD8 T cell deconvolution with that custom reference profile. 

# Third, I will fit a lme4 model for predicting CD8 T cell type abundance, based on 
# the infiltration groups Bogang defined (fixed effect), accounting for AOI patient
# identity (random effect).


## SETUP -----------------------------------------------------------------------
rm(list = ls())

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
library(multcomp)
library(parallel)
library(SpatialDecon)
library(ggdendro)
library(ggh4x)
library(Seurat)

# Data
meta <- readxl::read_excel("meta_cleaned_v5.xlsx") |> as.data.frame()
rownames(meta) <- meta$aoi_id
cts <- read.csv(file = "counts.csv", row.names = 1)
norm <- read.csv(file = "q3norm.csv", row.names = 1)


## CONSTRUCTING REFERENCE PROFILE ----------------------------------------------
# The data I downloaded is from this web page: https://singlecell.mdanderson.org/TCM/
# I downloaded the CD8 T cells dataset
refds <- readRDS("5_immune_deconvolution/CD8.rds")

# Summary of the dataset
table(refds$TissueType, refds$CancerType) |> as.data.frame() |> 
  ggplot() + geom_bar(mapping = aes(x = Var2, y = Freq, fill = Var1), stat = "identity") + 
  ggthemes::theme_few() + 
  theme(axis.text.x = element_text(angle = 45))

# How many cells are BRCA
(refds$CancerType == "BRCA") |> sum()
# [1] 1535

# There are not many BRCA T cells, so we will use BRCA and Breast together
# The cells are all from cancer patients, not healthy donors
# Subseting the data to what we want
refds <- subset(refds, subset = CancerType %in% c("BRCA", "Breast"))
refds |> dim()
# [1] 70777  3469

table(refds$cell.type)
# CD8_c0_t-Teff          CD8_c1_Tex         CD8_c2_Teff           CD8_c3_Tn         CD8_c4_Tstr         CD8_c5_Tisg          CD8_c6_Tcm 
# 1788                 286                 165                 307                  83                 263                  35 
# CD8_c7_p-Tex   CD8_c8_Teff_KLRG1         CD8_c9_Tsen  CD8_c10_Teff_CD244 CD8_c11_Teff_SEMA4A         CD8_c12_Trm     CD8_c13_Tn_TCF7 
# 30                  23                   9                  37                  29                 301                 113 

# Bogang said that he wanted to see the following cell types: 
# CD8_Tex (exhausted)	CD8_Tstr (stressed)	CD8_Tisg (IFN_response)	CD8_Teff (effector)	CD8_Tn (naïve)
# I will group all Teff together: CD8_c0_t-Teff, CD8_c2_Teff, CD8_c8_Teff_KLRG1, CD8_c10_Teff_CD244, and CD8_c11_Teff_SEMA4A
# I will group all Tex together: CD8_c1_Tex and CD8_c7_p-Tex
# I will group all Tn together: CD8_c3_Tn and CD8_c13_Tn_TCF7
# There is only one Tstr group: CD8_c4_Tstr
# There is only one Tisg group: CD8_c5_Tisg
# The memory cells seem distinct, so they will be left ungrouped
# The Tsen group has only 9 cells in it, so I will exclude these cells

# Before
DimPlot(refds, group.by = "cell.type") + ggprism::scale_color_prism()

# After
refds$type.new <- case_when(
  (refds$cell.type %in% c("CD8_c0_t-Teff", "CD8_c2_Teff", "CD8_c8_Teff_KLRG1", "CD8_c10_Teff_CD244", "CD8_c11_Teff_SEMA4A")) ~ "Teff", 
  (refds$cell.type %in% c("CD8_c1_Tex", "CD8_c7_p-Tex")) ~ "Tex", 
  (refds$cell.type %in% c("CD8_c3_Tn", "CD8_c13_Tn_TCF7")) ~ "Tn", 
  (refds$cell.type %in% c("CD8_c4_Tstr")) ~ "Tstr", 
  (refds$cell.type %in% c("CD8_c5_Tisg")) ~ "Tisg",
  (refds$cell.type %in% c("CD8_c6_Tcm")) ~ "Tcm",
  (refds$cell.type %in% c("CD8_c12_Trm")) ~ "Trm",
)
DimPlot(refds, group.by = "type.new") + ggprism::scale_color_prism()

# Making the reference profile
refds$cell_ID <- colnames(refds)
idx <- !(is.na(refds$type.new))
genestokeep <- openxlsx::read.xlsx("5_immune_deconvolution/41591_2023_2371_MOESM3_ESM.xlsx", sheet = "Table S3", startRow = 2) %$% gene
ref <- SpatialDecon::create_profile_matrix(mtx = refds@assays$RNA@counts[intersect(genestokeep, rownames(refds)),idx], 
                                           cellAnnots = refds@meta.data[idx,], cellTypeCol = "type.new",  
                                           cellNameCol = "cell_ID", 
                                           normalize = T)
glimpse(ref)
# num [1:484, 1:5] 0.637 6.894 4.455 18.814 7.383 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:484] "LMNA" "ZFP36" "RGCC" "CXCR4" ...
# ..$ : chr [1:5] "Tstr" "Teff" "Tisg" "Tn" ...


## DECONVOLUTION ---------------------------------------------------------------
# Subset to just immune AOIs
immnorm <- (norm |> as.matrix())[,meta$AOI_code == "immune"]
immcts <- (cts |> as.matrix())[,meta$AOI_code == "immune"]
immmeta <- meta[meta$AOI_code == "immune",]

# Background
bg = derive_GeoMx_background(norm = immnorm,
                             probepool = rep(1, nrow(immnorm)),
                             negnames = "NegProbe-WTX")

# Deconvolution
immdcvn <- spatialdecon(norm = immnorm, 
                        bg = bg, 
                        X = ref,
                        raw = immcts, 
                        cell_counts = immmeta$Number_of_cells_sampled)

# Visualization
# Something is always up with this same AOI: DSP.1001660011138.F.H07.dcc
ha1 <- HeatmapAnnotation(type = immmeta$infiltration_type, 
                         col = list(type = c("ieb"="dodgerblue", "iib"="purple", "iic"="magenta")))
mat <- immdcvn$prop_of_all
mat[is.na(mat)] <- NA
Heatmap(mat, 
        cluster_columns = F,
        show_column_names = T, 
        width = ncol(immdcvn$beta)*unit(5, "mm"),
        height = nrow(immdcvn$beta)*unit(5, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5),
        name = "Proportion", 
        top_annotation = ha1, 
        column_split = immmeta$infiltration_type) |> draw(column_title = "Immune, Proportion")

# Its library size, gene count, and area are not the smallest
meta$library_size <- colSums(cts)
meta$nGenes <- colSums(cts > 0)

# Without it
immnorm <- (norm |> as.matrix())[,meta$AOI_code == "immune" & meta$aoi_id != "DSP.1001660011138.F.H07.dcc"]
immcts <- (cts |> as.matrix())[,meta$AOI_code == "immune" & meta$aoi_id != "DSP.1001660011138.F.H07.dcc"]
immmeta <- meta[meta$AOI_code == "immune" & meta$aoi_id != "DSP.1001660011138.F.H07.dcc",]

bg = derive_GeoMx_background(norm = immnorm,
                             probepool = rep(1, nrow(immnorm)),
                             negnames = "NegProbe-WTX")

immdcvn <- spatialdecon(norm = immnorm, 
                        bg = bg, 
                        X = ref,
                        raw = immcts, 
                        cell_counts = immmeta$Number_of_cells_sampled)

ha2 <- HeatmapAnnotation(type = immmeta$infiltration_type, 
                         col = list(type = c("ieb"="dodgerblue", "iib"="purple", "iic"="magenta")))
mat2 <- immdcvn$prop_of_all
Heatmap(mat2,
        show_column_names = T, 
        width = ncol(immdcvn$beta)*unit(5, "mm"),
        height = nrow(immdcvn$beta)*unit(5, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5),
        name = "Proportion", 
        top_annotation = ha2, 
        column_split = immmeta$infiltration_type) |> draw(column_title = "Immune, Proportion")
Heatmap(mat2,
        show_column_names = T, 
        width = ncol(immdcvn$beta)*unit(5, "mm"),
        height = nrow(immdcvn$beta)*unit(5, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5),
        name = "Proportion", 
        top_annotation = ha2) |> draw(column_title = "Immune, Proportion")

# Barplot
d <- mat2 |> t() |> as.data.frame()
d$aoi_id <- rownames(d)
immdcvn_long <- tidyr::pivot_longer(data = d, 
                                    names_to = "celltype", 
                                    values_to = "prop",
                                    cols = rownames(immdcvn$prop_of_all),
                                    values_drop_na = T)

aoiorder <- immmeta |> arrange(infiltration_type) |> pull(aoi_id)
immdcvn_long$aoi_id <- factor(immdcvn_long$aoi_id, levels = aoiorder)

p1 <- ggplot(immdcvn_long) + 
  geom_bar(mapping = aes(x = aoi_id, y = prop, fill = celltype, group = aoi_id), stat = "identity", color = NA) + 
  ggprism::scale_fill_prism() + 
  theme_dendro() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(title = "Immune, Proportions")

p1 + geom_text(data=immmeta[aoiorder,],
               aes(label=infiltration_type, x=aoi_id, y=1.075, color=infiltration_type)) + 
  scale_color_manual(values = c("ieb"="dodgerblue", "iib"="purple", "iic"="magenta"))


## MODELING --------------------------------------------------------------------
# Adding data we need for the model
immdcvn_long$ncells <- plyr::mapvalues(x = immdcvn_long$aoi_id, from = immmeta$aoi_id, to = immmeta$Number_of_cells_sampled) |> as.character() |> as.numeric()
immdcvn_long$patient <- plyr::mapvalues(x = immdcvn_long$aoi_id, from = immmeta$aoi_id, to = immmeta$Patient_number)
immdcvn_long$group <- plyr::mapvalues(x = immdcvn_long$aoi_id, from = immmeta$aoi_id, to = immmeta$infiltration_type)

# Need to use random slopes model
table(immmeta$Patient_number, immmeta$infiltration_type)

# Proportion modeled, weighted by sample size
# Here are some sources: 
# https://stats.oarc.ucla.edu/stata/faq/how-does-one-do-regression-when-the-dependent-variable-is-a-proportion/#:~:text=Proportion%20data%20has%20values%20that,link%20and%20the%20binomial%20family.
# https://stats.stackexchange.com/questions/87956/how-to-apply-binomial-glmm-glmer-to-percentages-rather-than-yes-no-counts
proptestfunc <- function(ct, .contr) {
  tryCatch({
    
    testdata <- filter(immdcvn_long, celltype == ct)
    modout <- glmer(formula = prop ~ 1+group + (1+group|patient), weights = testdata$ncells, data = testdata, family = binomial(link = "logit"))
    testout = (multcomp::glht(model = modout, linfct = matrix(contrast_list[[.contr]], nrow = 1, byrow = T),
                           alternative = "two.sided",
                           rhs = 0)) |> summary(test = univariate())
    outs <- data.frame(
      contrast = .contr,
      celltype = ct,
      pval = testout[["test"]][["pvalues"]][1],
      Estimate = testout[["test"]][["coefficients"]][[1]], 
      tstat = testout[["test"]][["tstat"]][[1]], 
      sigma = testout[["test"]][["sigma"]][[1]]
    )
    return(outs)
    
  }, error = function(e){;})
}

mm <- model.matrix(~1+infiltration_type, data = immmeta)
contrast_list <- list(
  "ieb / iib" = mm[immmeta$infiltration_type == "ieb",] |> colMeans() - mm[immmeta$infiltration_type == "iib",] |> colMeans(),
  "ieb / iic" = mm[immmeta$infiltration_type == "ieb",] |> colMeans() - mm[immmeta$infiltration_type == "iic",] |> colMeans(),
  "iib / iic" = mm[immmeta$infiltration_type == "iib",] |> colMeans() - mm[immmeta$infiltration_type == "iic",] |> colMeans()
)

res <- lapply(X = names(contrast_list), FUN = purrr::map, .x = unique(immdcvn_long$celltype), .f = proptestfunc)
res <- bind_rows(res)
res <- group_by(res, contrast) |> mutate(fdr = p.adjust(p = pval, method = "BH"))
knitr::kable(res, format = "simple")

# Feeling better about these results... might need a better way to display them.

# contrast    celltype         pval      Estimate           tstat       sigma         fdr
# ----------  ---------  ----------  ------------  --------------  ----------  ----------
# ieb / iib   Tisg        0.2578228    -0.4614994      -1.1315520   0.4078463   0.2578228
# ieb / iib   Tn          0.0000000    -0.5169059     -74.3263975   0.0069545   0.0000000
# ieb / iib   Tex         0.0180229     0.4454697       2.3651469   0.1883476   0.0270344
# ieb / iic   Tisg        0.5089490    -2.3919428      -0.6604751   3.6215487   0.7634235
# ieb / iic   Tn          0.0000000   499.2624068   71786.4086945   0.0069548   0.0000000
# ieb / iic   Tex         0.7757621     0.2651651       0.2848460   0.9309069   0.7757621
# iib / iic   Tisg        0.5919139    -1.9304435      -0.5360646   3.6011399   0.8429069
# iib / iic   Tn          0.0000000   499.7793127   50814.2737767   0.0098354   0.0000000
# iib / iic   Tex         0.8429069    -0.1803046      -0.1981765   0.9098179   0.8429069


# I think using difference in average proportion will work for now
mm <- model.matrix(~0+infiltration_type, data = immmeta)
contrast_list <- list(
  "ieb / iib" = mm[immmeta$infiltration_type == "ieb",] |> colMeans() - mm[immmeta$infiltration_type == "iib",] |> colMeans(),
  "ieb / iic" = mm[immmeta$infiltration_type == "ieb",] |> colMeans() - mm[immmeta$infiltration_type == "iic",] |> colMeans(),
  "iib / iic" = mm[immmeta$infiltration_type == "iib",] |> colMeans() - mm[immmeta$infiltration_type == "iic",] |> colMeans()
)
meanprops <- immdcvn_long |> group_by(celltype, group) |> summarise(mean_prop = mean(prop))
meanprops <- tidyr::pivot_wider(meanprops, id_cols = "group", names_from = "celltype", values_from = "mean_prop")
calcmeandiff <- function(.ct, .contr) {
  return(data.frame(
    celltype = .ct, 
    contrast = .contr,
    mean_diff = (meanprops[[.ct]]*contrast_list[[.contr]]) |> sum()))
}

meansres <- lapply(X = names(contrast_list), FUN = purrr::map, .x = unique(immdcvn_long$celltype), .f = calcmeandiff)
meansres <- bind_rows(meansres)
res <- inner_join(x = res, y = meansres, by = c("celltype", "contrast"))

# Final visualizaton for testing
ggplot(data = res) + 
  geom_bar(mapping = aes(y = celltype, x = mean_diff, fill = fdr < 0.05), stat = "identity") + 
  facet_grid(.~contrast) +
  scale_fill_manual(values = c("blue", "red")) + 
  ggthemes::theme_few() + 
  geom_vline(xintercept = 0, linewidth = 0.2)

dev.off()

write.csv(res, "immune_deconvolution_differential_abundance_testing_results.csv")

## NOTES -----------------------------------------------------------------------
# This took 7 hours
# Still want to double check with martin about the test

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
#   [1] parallel  grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] Seurat_5.2.0          SeuratObject_5.0.2    sp_2.1-4              ggh4x_0.3.0           ggdendro_0.2.0        SpatialDecon_1.14.0  
# [7] multcomp_1.4-26       TH.data_1.1-2         MASS_7.3-61           survival_3.7-0        mvtnorm_1.3-1         lmerTest_3.1-3       
# [13] lme4_1.1-35.5         Matrix_1.7-0          patchwork_1.3.0       magrittr_2.0.3        ggbiplot_0.6.2        edgeR_4.2.1          
# [19] limma_3.60.4          ComplexHeatmap_2.20.0 ggplot2_3.5.1         stringr_1.5.1         dplyr_1.1.4          
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22         splines_4.4.1            later_1.3.2              tibble_3.2.1             R.oo_1.26.0             
# [6] cellranger_1.1.0         polyclip_1.10-7          fastDummies_1.7.5        lifecycle_1.0.4          doParallel_1.0.17       
# [11] globals_0.16.3           lattice_0.22-6           openxlsx_4.2.7.1         GeomxTools_3.8.0         plotly_4.10.4           
# [16] httpuv_1.6.15            sctransform_0.4.1        zip_2.3.1                spam_2.10-0              spatstat.sparse_3.1-0   
# [21] reticulate_1.39.0        cowplot_1.1.3            pbapply_1.7-2            minqa_1.2.8              RColorBrewer_1.1-3      
# [26] abind_1.4-8              zlibbioc_1.50.0          EnvStats_3.0.0           Rtsne_0.17               R.cache_0.16.0          
# [31] purrr_1.0.2              R.utils_2.12.3           BiocGenerics_0.50.0      sandwich_3.1-1           circlize_0.4.16         
# [36] GenomeInfoDbData_1.2.12  IRanges_2.38.1           S4Vectors_0.42.1         ggrepel_0.9.6            irlba_2.3.5.1           
# [41] spatstat.utils_3.1-2     listenv_0.9.1            pheatmap_1.0.12          goftest_1.2-3            RSpectra_0.16-2         
# [46] spatstat.random_3.3-2    fitdistrplus_1.2-2       parallelly_1.38.0        codetools_0.2-20         tidyselect_1.2.1        
# [51] shape_1.4.6.1            UCSC.utils_1.0.0         farver_2.1.2             matrixStats_1.4.1        stats4_4.4.1            
# [56] spatstat.explore_3.3-4   jsonlite_1.8.9           GetoptLong_1.0.5         progressr_0.14.0         ggridges_0.5.6          
# [61] iterators_1.0.14         systemfonts_1.1.0        foreach_1.5.2            tools_4.4.1              ica_1.0-3               
# [66] Rcpp_1.0.13              glue_1.8.0               gridExtra_2.3            xfun_0.47                ggthemes_5.1.0          
# [71] GenomeInfoDb_1.40.1      withr_3.0.1              numDeriv_2016.8-1.1      fastmap_1.2.0            NanoStringNCTools_1.12.0
# [76] GGally_2.2.1             repmis_0.5               boot_1.3-31              fansi_1.0.6              digest_0.6.37           
# [81] R6_2.5.1                 mime_0.12                ggprism_1.0.5            colorspace_2.1-1         scattermore_1.2         
# [86] tensor_1.5               spatstat.data_3.1-4      R.methodsS3_1.8.2        utf8_1.2.4               tidyr_1.3.1             
# [91] generics_0.1.3           data.table_1.16.0        httr_1.4.7               htmlwidgets_1.6.4        ggstats_0.7.0           
# [96] uwot_0.2.2               pkgconfig_2.0.3          gtable_0.3.5             lmtest_0.9-40            XVector_0.44.0          
# [101] htmltools_0.5.8.1        dotCall64_1.1-1          clue_0.3-65              scales_1.3.0             Biobase_2.64.0          
# [106] png_0.1-8                logNormReg_0.5-0         spatstat.univar_3.1-1    knitr_1.48               rstudioapi_0.16.0       
# [111] reshape2_1.4.4           rjson_0.2.23             uuid_1.2-1               nlme_3.1-166             nloptr_2.1.1            
# [116] zoo_1.8-12               GlobalOptions_0.1.2      KernSmooth_2.23-24       miniUI_0.1.1.1           vipor_0.4.7             
# [121] pillar_1.9.0             vctrs_0.6.5              RANN_2.6.2               promises_1.3.0           xtable_1.8-4            
# [126] cluster_2.1.6            beeswarm_0.4.0           magick_2.8.5             cli_3.6.3                locfit_1.5-9.6          
# [131] compiler_4.4.1           rlang_1.1.4              crayon_1.5.3             future.apply_1.11.2      labeling_0.4.3          
# [136] plyr_1.8.9               ggbeeswarm_0.7.2         ggiraph_0.8.10           stringi_1.8.4            deldir_2.0-4            
# [141] viridisLite_0.4.2        munsell_0.5.1            Biostrings_2.72.1        lazyeval_0.2.2           spatstat.geom_3.3-5     
# [146] RcppHNSW_0.6.0           future_1.34.0            statmod_1.5.0            shiny_1.9.1              ROCR_1.0-11             
# [151] igraph_2.0.3             readxl_1.4.3            
# 
# 
