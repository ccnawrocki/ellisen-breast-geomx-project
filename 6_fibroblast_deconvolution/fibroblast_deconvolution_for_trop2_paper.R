set.seed(2001)
pdf("fibroblast_deconvolution_for_trop2_paper.pdf", width = 12, height = 8)

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
norm <- read.csv(file = "q3norm_nonlog.csv", row.names = 1) # The data needs to be linear (not log-transformed)

# Reference profile
ref <- read.csv("6_fibroblast_deconvolution/Custom_profileMatrix.csv", row.names = 1) |> as.matrix()

# Subset to just fibroblast AOIs
fibnorm <- (norm |> as.matrix())[,meta$AOI_code == "fibroblast"]
fibcts <- (cts |> as.matrix())[,meta$AOI_code == "fibroblast"]
fibmeta <- meta[meta$AOI_code == "fibroblast",]

## DECON -----------------------------------------------------------------------
# Background
bg = derive_GeoMx_background(norm = fibnorm,
                             probepool = rep(1, nrow(fibnorm)),
                             negnames = "NegProbe-WTX")

# Deconvolution
fibdcvn <- spatialdecon(norm = fibnorm, 
                        bg = bg, 
                        X = ref,
                        raw = fibcts
)

# Visualization
ha1 <- HeatmapAnnotation(type = fibmeta$infiltration_type, 
                         col = list(type = c("ieb"="dodgerblue", "iib"="purple", "idb"="orange")))
mat <- fibdcvn$prop_of_all

Heatmap(mat,
        show_column_names = T, 
        width = ncol(fibdcvn$beta)*unit(5, "mm"),
        height = nrow(fibdcvn$beta)*unit(5, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5),
        name = "Proportion", 
        top_annotation = ha1, 
        column_split = fibmeta$infiltration_type) |> draw(column_title = "Fibroblast, Proportion")

## MODELING --------------------------------------------------------------------
d <- mat |> t() |> as.data.frame()
d$aoi_id <- rownames(d)
fibdcvn_long <- tidyr::pivot_longer(data = d, 
                                    names_to = "celltype", 
                                    values_to = "prop",
                                    cols = rownames(fibdcvn$prop_of_all),
                                    values_drop_na = T)
aoiorder <- fibmeta |> arrange(infiltration_type) |> pull(aoi_id)
fibdcvn_long$aoi_id <- factor(fibdcvn_long$aoi_id, levels = aoiorder)

# We need discrete values to do the modeling. Thus, we will have to round the cell
# counts to the nearest whole number.

# Adding some other data we need for the model
fibdcvn_long$ncells <- plyr::mapvalues(x = fibdcvn_long$aoi_id, from = fibmeta$aoi_id, to = fibmeta$Number_of_cells_sampled) |> as.character() |> as.numeric()
fibdcvn_long$patient <- plyr::mapvalues(x = fibdcvn_long$aoi_id, from = fibmeta$aoi_id, to = fibmeta$Patient_number)
fibdcvn_long$group <- plyr::mapvalues(x = fibdcvn_long$aoi_id, from = fibmeta$aoi_id, to = fibmeta$infiltration_type)

# Getting the discrete counts
fibdcvn_long$count <- (fibdcvn_long$ncells*fibdcvn_long$prop) |> round(digits = 0)

# Should use random slopes model, since some patients exist across the groups, 
# according to NanoString. But, looking here, the study design is very unbalanced. 
# Also, there are only 8 patients, and patient 4 only has 1 aoi. Using a random
# slopes model will likely overfit the data or it will fail to converge.
table(fibmeta$Patient_number, fibmeta$infiltration_type)
#   idb ieb iib
# 1   0   7   0
# 2   0   1   3
# 3   3   0   0
# 4   0   0   1
# 5   0   0   4
# 6   1   2   1
# 7   0   3   0
# 8   0   1   1

# Another thing to note is that when looking at the bar plots shown above, it just 
# does not seem like the makeup of the aois within patient groups is very constant.
# So, we will test a bunch of models, but I think that a simpler model that does
# not lead to convergence issues will be best.

# Proportion modeled, weighted by sample size
# Here are some sources: 
# https://stats.oarc.ucla.edu/stata/faq/how-does-one-do-regression-when-the-dependent-variable-is-a-proportion/#:~:text=Proportion%20data%20has%20values%20that,link%20and%20the%20binomial%20family.
# https://stats.stackexchange.com/questions/87956/how-to-apply-binomial-glmm-glmer-to-percentages-rather-than-yes-no-counts

# In the end, I think that a simple poisson regression model is best. Note that 
# I am using a Wald test via the multcomp package.
proptestfunc <- function(ct, .contr) {
  tryCatch({
    
    testdata <- filter(fibdcvn_long, celltype == ct)
    
    # Poisson regression, patient excluded from model
    modout <- glm(formula = count ~ 1+group+offset(log(ncells)), data = testdata, family = poisson(link = "log"))
    
    testout <- (multcomp::glht(model = modout, linfct = matrix(contrast_list[[.contr]], nrow = 1, byrow = T),
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

fibmeta$infiltration_type %<>% as.factor()
all(levels(fibmeta$infiltration_type) == levels(fibdcvn_long$group))
# [1] TRUE

mm <- model.matrix(~1+infiltration_type, data = fibmeta)
contrast_list <- list(
  "idb vs iib and ieb" = mm[fibmeta$infiltration_type == "idb",] |> colMeans() - (mm[fibmeta$infiltration_type == "iib",] |> colMeans() + mm[fibmeta$infiltration_type == "ieb",] |> colMeans())/2
)

# Note that there are no warning messages when using the simplest poisson model.
res <- lapply(X = names(contrast_list), FUN = purrr::map, .x = unique(fibdcvn_long$celltype), .f = proptestfunc)
res <- bind_rows(res)
res <- group_by(res, contrast) |> mutate(fdr = p.adjust(p = pval, method = "BH"))
res$padj <- p.adjust(p = res$pval, method = "bonferroni")
knitr::kable(res, format = "simple")

# contrast             celltype          pval      Estimate        tstat          sigma         fdr       padj
# -------------------  ----------  ----------  ------------  -----------  -------------  ----------  ---------
# ibd vs iib and ieb   iCAF         0.9994269    -9.4601591   -0.0007183   1.316985e+04   0.9999991   1.00e+00
# ibd vs iib and ieb   mCAF         0.0000000     1.6179919    9.2039499   1.757932e-01   0.0000000   0.00e+00 *
# ibd vs iib and ieb   hsp_tpCAF    0.9999991    -0.1921863   -0.0000012   1.651461e+05   0.9999991   1.00e+00
# ibd vs iib and ieb   IDO_CAF      0.9867631   -17.1595087   -0.0165908   1.034280e+03   0.9999991   1.00e+00
# ibd vs iib and ieb   apCAF        0.0000000    -1.2516910   -7.3158900   1.710921e-01   0.0000000   0.00e+00 *
# ibd vs iib and ieb   dCAF         0.0000000     1.2244045    6.8921757   1.776514e-01   0.0000000   0.00e+00 *
# ibd vs iib and ieb   Pericyte     0.9990501    -9.5094652   -0.0011905   7.987917e+03   0.9999991   1.00e+00
# ibd vs iib and ieb   tpCAF        0.0000025     1.7567644    4.7117229   3.728497e-01   0.0000061   2.46e-05 *
# ibd vs iib and ieb   vCAF         0.9994269    -9.4601591   -0.0007183   1.316985e+04   0.9999991   1.00e+00
# ibd vs iib and ieb   rCAF         0.8784055    -0.1727508   -0.1529909   1.129158e+00   0.9999991   1.00e+00

## VIZ -------------------------------------------------------------------------
mm <- model.matrix(~0+infiltration_type, data = fibmeta)
contrast_list <- list(
  "idb vs iib and ieb" = mm[fibmeta$infiltration_type == "idb",] |> colMeans() - (mm[fibmeta$infiltration_type == "iib",] |> colMeans() + mm[fibmeta$infiltration_type == "ieb",] |> colMeans())/2
)
meanprops <- fibdcvn_long |> group_by(celltype, group) |> summarise(mean_prop = mean(count/ncells))
meanprops <- tidyr::pivot_wider(meanprops, id_cols = "group", names_from = "celltype", values_from = "mean_prop")
calcmeandiff <- function(.ct, .contr) {
  return(data.frame(
    celltype = .ct, 
    contrast = .contr,
    mean_diff = (meanprops[[.ct]]*contrast_list[[.contr]]) |> sum()))
}

meansres <- lapply(X = names(contrast_list), FUN = purrr::map, .x = unique(fibdcvn_long$celltype), .f = calcmeandiff)
meansres <- bind_rows(meansres)
res <- inner_join(x = res, y = meansres, by = c("celltype", "contrast"))

# Final visualizaton for testing
# Bonferroni correction is probably the way to go here. 
ggplot(data = res) + 
  geom_bar(mapping = aes(y = celltype, x = mean_diff, fill = padj < 0.05), stat = "identity") + 
  facet_grid(.~contrast) +
  scale_fill_manual(values = c("blue", "red")) + 
  ggthemes::theme_few() + 
  geom_vline(xintercept = 0, linewidth = 0.2)

write.csv(res, "fibroblast_deconvolution_differential_abundance_testing_results_for_trop2_paper.csv")

dev.off()

## SESSION ---------------------------------------------------------------------

# sessionInfo()
# 
# R version 4.4.1 (2024-06-14)
# Platform: aarch64-apple-darwin20.0.0
# Running under: macOS 15.3.1
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
# [11] globals_0.16.3           lattice_0.22-6           GeomxTools_3.8.0         plotly_4.10.4            httpuv_1.6.15           
# [16] sctransform_0.4.1        spam_2.10-0              spatstat.sparse_3.1-0    reticulate_1.39.0        cowplot_1.1.3           
# [21] pbapply_1.7-2            minqa_1.2.8              RColorBrewer_1.1-3       abind_1.4-8              zlibbioc_1.50.0         
# [26] EnvStats_3.0.0           Rtsne_0.17               R.cache_0.16.0           purrr_1.0.2              R.utils_2.12.3          
# [31] BiocGenerics_0.50.0      sandwich_3.1-1           circlize_0.4.16          GenomeInfoDbData_1.2.12  IRanges_2.38.1          
# [36] S4Vectors_0.42.1         ggrepel_0.9.6            irlba_2.3.5.1            spatstat.utils_3.1-2     listenv_0.9.1           
# [41] pheatmap_1.0.12          goftest_1.2-3            RSpectra_0.16-2          spatstat.random_3.3-2    fitdistrplus_1.2-2      
# [46] parallelly_1.38.0        codetools_0.2-20         tidyselect_1.2.1         shape_1.4.6.1            UCSC.utils_1.0.0        
# [51] farver_2.1.2             matrixStats_1.4.1        stats4_4.4.1             spatstat.explore_3.3-4   jsonlite_1.8.9          
# [56] GetoptLong_1.0.5         progressr_0.14.0         emmeans_1.10.6           ggridges_0.5.6           iterators_1.0.14        
# [61] systemfonts_1.1.0        foreach_1.5.2            tools_4.4.1              ica_1.0-3                Rcpp_1.0.13             
# [66] glue_1.8.0               gridExtra_2.3            xfun_0.47                ggthemes_5.1.0           GenomeInfoDb_1.40.1     
# [71] withr_3.0.1              numDeriv_2016.8-1.1      fastmap_1.2.0            NanoStringNCTools_1.12.0 GGally_2.2.1            
# [76] repmis_0.5               boot_1.3-31              fansi_1.0.6              digest_0.6.37            estimability_1.5.1      
# [81] R6_2.5.1                 mime_0.12                colorspace_2.1-1         scattermore_1.2          tensor_1.5              
# [86] spatstat.data_3.1-4      R.methodsS3_1.8.2        utf8_1.2.4               tidyr_1.3.1              generics_0.1.3          
# [91] data.table_1.16.0        httr_1.4.7               htmlwidgets_1.6.4        ggstats_0.7.0            uwot_0.2.2              
# [96] pkgconfig_2.0.3          gtable_0.3.5             lmtest_0.9-40            XVector_0.44.0           htmltools_0.5.8.1       
# [101] dotCall64_1.1-1          clue_0.3-65              scales_1.3.0             Biobase_2.64.0           png_0.1-8               
# [106] logNormReg_0.5-0         spatstat.univar_3.1-1    knitr_1.48               rstudioapi_0.16.0        reshape2_1.4.4          
# [111] rjson_0.2.23             uuid_1.2-1               coda_0.19-4.1            nlme_3.1-166             nloptr_2.1.1            
# [116] zoo_1.8-12               GlobalOptions_0.1.2      KernSmooth_2.23-24       miniUI_0.1.1.1           vipor_0.4.7             
# [121] pillar_1.9.0             vctrs_0.6.5              RANN_2.6.2               promises_1.3.0           xtable_1.8-4            
# [126] cluster_2.1.6            beeswarm_0.4.0           magick_2.8.5             cli_3.6.3                locfit_1.5-9.6          
# [131] compiler_4.4.1           rlang_1.1.4              crayon_1.5.3             future.apply_1.11.2      labeling_0.4.3          
# [136] plyr_1.8.9               ggbeeswarm_0.7.2         ggiraph_0.8.10           stringi_1.8.4            deldir_2.0-4            
# [141] viridisLite_0.4.2        munsell_0.5.1            Biostrings_2.72.1        lazyeval_0.2.2           spatstat.geom_3.3-5     
# [146] RcppHNSW_0.6.0           future_1.34.0            statmod_1.5.0            shiny_1.9.1              ROCR_1.0-11             
# [151] igraph_2.0.3             readxl_1.4.3            
