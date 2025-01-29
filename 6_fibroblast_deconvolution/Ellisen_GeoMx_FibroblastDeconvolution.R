##### Ellisen GeoMx Fibroblast Deconvolution #####
### Cole Nawrocki ###

set.seed(2001)
pdf("fibroblast_deconvolution.pdf", width = 12, height = 8)

## PLAN ------------------------------------------------------------------------
# Using the data from this paper: https://www.nature.com/articles/s41467-023-39762-1#Sec20
# to do deconvolution on the fibroblast segments only. 

# First, I will read the paper, download the data, and make a custom reference profile, 
# using the SpatialDecon package. 

# Second, I will perform fibroblast deconvolution with that custom reference profile. 

# Third, I will fit a lme4 model for predicting fibroblast cell type abundance, based on 
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
refds <- readRDS("6_fibroblast_deconvolution/BREAST_fibro_tumour.rds")

# Summary of the dataset
table(refds$CAFtype, refds$dataset) |> as.data.frame() |> 
  ggplot() + geom_bar(mapping = aes(x = Var2, y = Freq, fill = Var1), stat = "identity") + 
  ggthemes::theme_few() + 
  theme(axis.text.x = element_text(angle = 45))

# The cells are all from cancer patients, not healthy donors
# Subseting the data to what we want
refds |> dim()
# [1] 17101 16704

table(refds$CAFtype)
# mCAF      iCAF      vCAF  Pericyte     tpCAF hsp_tpCAF   IDO_CAF     apCAF      rCAF      dCAF 
# 4525      3439      2886      2389       786       722       665       793       373       126 

# Bogang did not say which types he wanted, so I will leave them all in for now.
DimPlot(refds, group.by = "CAFtype") + ggprism::scale_color_prism()

# Making the reference profile
refds$cell_ID <- colnames(refds)
genestokeep <- read.csv("6_fibroblast_deconvolution/41467_2023_39762_MOESM6_ESM.csv") |> 
  filter(p_val_adj < 0.05) |> group_by(cluster) |> top_n(n = 50, wt = avg_log2FC) |> pull(gene)
ref <- SpatialDecon::create_profile_matrix(mtx = refds@assays$RNA@counts[intersect(genestokeep, rownames(refds)),], 
                                           cellAnnots = refds@meta.data, cellTypeCol = "CAFtype",  
                                           cellNameCol = "cell_ID", 
                                           normalize = T)
glimpse(ref)
# num [1:455, 1:10] 1.99 21.43 239.43 2.92 1.91 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:455] "MMP11" "POSTN" "COL1A1" "COMP" ...
# ..$ : chr [1:10] "iCAF" "mCAF" "hsp_tpCAF" "IDO_CAF" ...

## DECONVOLUTION ---------------------------------------------------------------
# Subset to just fibroblast AOIs
fibnorm <- (norm |> as.matrix())[,meta$AOI_code == "fibroblast"]
fibcts <- (cts |> as.matrix())[,meta$AOI_code == "fibroblast"]
fibmeta <- meta[meta$AOI_code == "fibroblast",]

# Background
bg = derive_GeoMx_background(norm = fibnorm,
                             probepool = rep(1, nrow(fibnorm)),
                             negnames = "NegProbe-WTX")

# Deconvolution
fibdcvn <- spatialdecon(norm = fibnorm, 
                        bg = bg, 
                        X = ref,
                        raw = fibcts, 
                        cell_counts = fibmeta$Number_of_cells_sampled)

# Visualization
ha1 <- HeatmapAnnotation(type = fibmeta$infiltration_type, 
                         col = list(type = c("ieb"="dodgerblue", "iib"="purple", "idb"="orange")))
mat <- fibdcvn$prop_of_all
Heatmap(mat, 
        cluster_columns = F,
        show_column_names = T, 
        width = ncol(fibdcvn$beta)*unit(5, "mm"),
        height = nrow(fibdcvn$beta)*unit(5, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5),
        name = "Proportion", 
        top_annotation = ha1, 
        column_split = fibmeta$infiltration_type) |> draw(column_title = "Fibroblast, Proportion")
Heatmap(mat,
        show_column_names = T, 
        width = ncol(fibdcvn$beta)*unit(5, "mm"),
        height = nrow(fibdcvn$beta)*unit(5, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5),
        name = "Proportion", 
        top_annotation = ha1, 
        column_split = fibmeta$infiltration_type) |> draw(column_title = "Fibroblast, Proportion")
Heatmap(mat,
        show_column_names = T, 
        width = ncol(fibdcvn$beta)*unit(5, "mm"),
        height = nrow(fibdcvn$beta)*unit(5, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5),
        name = "Proportion", 
        top_annotation = ha1) |> draw(column_title = "Fibroblast, Proportion")


# Barplot
d <- mat |> t() |> as.data.frame()
d$aoi_id <- rownames(d)
fibdcvn_long <- tidyr::pivot_longer(data = d, 
                                    names_to = "celltype", 
                                    values_to = "prop",
                                    cols = rownames(fibdcvn$prop_of_all),
                                    values_drop_na = T)

aoiorder <- fibmeta |> arrange(infiltration_type) |> pull(aoi_id)
fibdcvn_long$aoi_id <- factor(fibdcvn_long$aoi_id, levels = aoiorder)

p1 <- ggplot(fibdcvn_long) + 
  geom_bar(mapping = aes(x = aoi_id, y = prop, fill = celltype, group = aoi_id), stat = "identity", color = NA) + 
  ggprism::scale_fill_prism() + 
  theme_dendro() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(title = "Fibroblast, Proportions")

p1 + geom_text(data=fibmeta[aoiorder,],
               aes(label=infiltration_type, x=aoi_id, y=1.075, color=infiltration_type)) + 
  scale_color_manual(values = c("ieb"="dodgerblue", "iib"="purple", "idb"="orange"))

## MODELING --------------------------------------------------------------------
# Adding data we need for the model
fibdcvn_long$ncells <- plyr::mapvalues(x = fibdcvn_long$aoi_id, from = fibmeta$aoi_id, to = fibmeta$Number_of_cells_sampled) |> as.character() |> as.numeric()
fibdcvn_long$patient <- plyr::mapvalues(x = fibdcvn_long$aoi_id, from = fibmeta$aoi_id, to = fibmeta$Patient_number)
fibdcvn_long$group <- plyr::mapvalues(x = fibdcvn_long$aoi_id, from = fibmeta$aoi_id, to = fibmeta$infiltration_type)

# Need to use random slopes model
table(fibmeta$Patient_number, fibmeta$infiltration_type)

# Proportion modeled, weighted by sample size
# Here are some sources: 
# https://stats.oarc.ucla.edu/stata/faq/how-does-one-do-regression-when-the-dependent-variable-is-a-proportion/#:~:text=Proportion%20data%20has%20values%20that,link%20and%20the%20binomial%20family.
# https://stats.stackexchange.com/questions/87956/how-to-apply-binomial-glmm-glmer-to-percentages-rather-than-yes-no-counts
proptestfunc <- function(ct, .contr) {
  tryCatch({
    
    testdata <- filter(fibdcvn_long, celltype == ct)
    modout <- glmer(formula = round(ncells*prop, digits = 0) ~ 1+group + (1+group|patient), weights = testdata$ncells, data = testdata, family = binomial(link = "logit"))
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

mm <- model.matrix(~1+infiltration_type, data = fibmeta)
contrast_list <- list(
  "ieb / iib" = mm[fibmeta$infiltration_type == "ieb",] |> colMeans() - mm[fibmeta$infiltration_type == "iib",] |> colMeans(),
  "ieb / idb" = mm[fibmeta$infiltration_type == "ieb",] |> colMeans() - mm[fibmeta$infiltration_type == "idb",] |> colMeans(),
  "iib / idb" = mm[fibmeta$infiltration_type == "iib",] |> colMeans() - mm[fibmeta$infiltration_type == "idb",] |> colMeans()
)

proptestfunc(ct = "fibroblast", .contr = contrast_list$`ieb / iib`)
testdata <- filter(fibdcvn_long, celltype == "fibroblast")
modout <- glmer(formula = round(ncells*prop, digits = 0) ~ 1+group + (1+group|patient), weights = testdata$ncells, data = testdata, family = binomial(link = "logit"))
res <- lapply(X = names(contrast_list), FUN = purrr::map, .x = unique(fibdcvn_long$celltype), .f = proptestfunc)
res <- bind_rows(res)
res <- group_by(res, contrast) |> mutate(fdr = p.adjust(p = pval, method = "BH"))
knitr::kable(res, format = "simple")

# Feeling better about these results... might need a better way to display them.

# contrast    celltype          pval        Estimate           tstat          sigma         fdr
# ----------  ----------  ----------  --------------  --------------  -------------  ----------
# ieb / iib   mCAF         0.7710380      -2.2196673   -2.910175e-01   7.627265e+00   0.9999913
# ieb / iib   hsp_tpCAF    0.6228166       3.0266817    4.918624e-01   6.153513e+00   0.9965066
# ieb / iib   IDO_CAF      0.2207077      -1.2710705   -1.224649e+00   1.037906e+00   0.4763920
# ieb / iib   apCAF        0.2381960      -0.9061108   -1.179508e+00   7.682108e-01   0.4763920
# ieb / iib   dCAF         0.0074988       2.2530026    2.673843e+00   8.426084e-01   0.0299950
# ieb / iib   Pericyte     0.9999913      35.9069628    1.100000e-05   3.274577e+06   0.9999913
# ieb / iib   tpCAF        0.0000000    2364.1384586    3.308021e+05   7.146700e-03   0.0000000
# ieb / iib   rCAF         0.9999793     -75.2097531   -2.590000e-05   2.904085e+06   0.9999913
# ieb / idb   mCAF         0.0775484      -9.4992293   -1.765092e+00   5.381721e+00   0.1550968
# ieb / idb   hsp_tpCAF    0.1840084      -2.9458924   -1.328514e+00   2.217434e+00   0.2944134
# ieb / idb   IDO_CAF      0.9702384      20.3259050    3.730930e-02   5.447952e+02   0.9999863
# ieb / idb   apCAF        0.0026661       1.7277293    3.003829e+00   5.751756e-01   0.0071095
# ieb / idb   dCAF         0.0002882      -1.4757300   -3.625714e+00   4.070178e-01   0.0011527
# ieb / idb   Pericyte     0.9999774      28.1416107    2.830000e-05   9.946694e+05   0.9999863
# ieb / idb   tpCAF        0.0000000      -1.2269547   -2.428663e+02   5.052000e-03   0.0000000
# ieb / idb   rCAF         0.9999863     -51.5843856   -1.710000e-05   3.009823e+06   0.9999863
# iib / idb   mCAF         0.1683372      -7.2795621   -1.377566e+00   5.284364e+00   0.3366744
# iib / idb   hsp_tpCAF    0.2949366      -5.9725740   -1.047353e+00   5.702540e+00   0.4718986
# iib / idb   IDO_CAF      0.9683782      21.5969755    3.964240e-02   5.447952e+02   0.9999982
# iib / idb   apCAF        0.0000005       2.6338401    5.031159e+00   5.235056e-01   0.0000020
# iib / idb   dCAF         0.0000199      -3.7287326   -4.266320e+00   8.739928e-01   0.0000530
# iib / idb   Pericyte     0.9999982      -7.7653521   -2.300000e-06   3.422313e+06   0.9999982
# iib / idb   tpCAF        0.0000000   -2365.3654133   -4.679293e+05   5.055000e-03   0.0000000
# iib / idb   rCAF         0.9999762      23.6253676    2.990000e-05   7.907757e+05   0.9999982

# I think using difference in average proportion will work for now
mm <- model.matrix(~0+infiltration_type, data = fibmeta)
contrast_list <- list(
  "ieb / iib" = mm[fibmeta$infiltration_type == "ieb",] |> colMeans() - mm[fibmeta$infiltration_type == "iib",] |> colMeans(),
  "ieb / idb" = mm[fibmeta$infiltration_type == "ieb",] |> colMeans() - mm[fibmeta$infiltration_type == "idb",] |> colMeans(),
  "iib / idb" = mm[fibmeta$infiltration_type == "iib",] |> colMeans() - mm[fibmeta$infiltration_type == "idb",] |> colMeans()
)
meanprops <- fibdcvn_long |> group_by(celltype, group) |> summarise(mean_prop = mean(prop))
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
ggplot(data = res) + 
  geom_bar(mapping = aes(y = celltype, x = mean_diff, fill = fdr < 0.05), stat = "identity") + 
  facet_grid(.~contrast) +
  scale_fill_manual(values = c("blue", "red")) + 
  ggthemes::theme_few() + 
  geom_vline(xintercept = 0, linewidth = 0.2)

dev.off()

write.csv(res, "fibroblast_deconvolution_differential_abundance_testing_results.csv")

## NOTES -----------------------------------------------------------------------
# This took 2 hours

