## Differential Abundance Testing Checks ##
### Cole Nawrocki ###

set.seed(2001)
rm(list = ls())
.libPaths()
# [1] "/Users/ccn22/micromamba/envs/geomx-env/lib/R/library"

library(dplyr)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(magrittr)
library(patchwork)
library(lme4)
library(lmerTest)
library(multcomp)
library(parallel)
library(SpatialDecon)
library(Seurat)

meta <- readxl::read_excel("meta_cleaned_v5.xlsx") |> as.data.frame()
rownames(meta) <- meta$aoi_id
cts <- read.csv(file = "counts.csv", row.names = 1)
norm <- read.csv(file = "q3norm.csv", row.names = 1)

refds <- readRDS("6_fibroblast_deconvolution/BREAST_fibro_tumour.rds")
refds$cell_ID <- colnames(refds)
genestokeep <- read.csv("6_fibroblast_deconvolution/41467_2023_39762_MOESM6_ESM.csv") |> 
  filter(p_val_adj < 0.05) |> group_by(cluster) |> top_n(n = 50, wt = avg_log2FC) |> pull(gene)
ref <- SpatialDecon::create_profile_matrix(mtx = refds@assays$RNA@counts[intersect(genestokeep, rownames(refds)),], 
                                           cellAnnots = refds@meta.data, cellTypeCol = "CAFtype",  
                                           cellNameCol = "cell_ID", 
                                           normalize = T)

fibnorm <- (norm |> as.matrix())[,meta$AOI_code == "fibroblast"]
fibcts <- (cts |> as.matrix())[,meta$AOI_code == "fibroblast"]
fibmeta <- meta[meta$AOI_code == "fibroblast",]

bg = derive_GeoMx_background(norm = fibnorm,
                             probepool = rep(1, nrow(fibnorm)),
                             negnames = "NegProbe-WTX")
fibdcvn <- spatialdecon(norm = fibnorm, 
                        bg = bg, 
                        X = ref,
                        raw = fibcts, 
                        cell_counts = fibmeta$Number_of_cells_sampled)

# Binomial logistic regression with prop as the response variable 
################################################################################
mat <- fibdcvn$prop_of_all
d <- mat |> t() |> as.data.frame()
d$aoi_id <- rownames(d)
fibdcvn_long <- tidyr::pivot_longer(data = d, 
                                    names_to = "celltype", 
                                    values_to = "prop",
                                    cols = rownames(fibdcvn$prop_of_all),
                                    values_drop_na = T)

fibdcvn_long$ncells <- plyr::mapvalues(x = fibdcvn_long$aoi_id, from = fibmeta$aoi_id, to = fibmeta$Number_of_cells_sampled) |> as.character() |> as.numeric()
fibdcvn_long$patient <- plyr::mapvalues(x = fibdcvn_long$aoi_id, from = fibmeta$aoi_id, to = fibmeta$Patient_number)
fibdcvn_long$group <- plyr::mapvalues(x = fibdcvn_long$aoi_id, from = fibmeta$aoi_id, to = fibmeta$infiltration_type)

glimpse(fibdcvn_long)
# Rows: 280
# Columns: 6
# $ aoi_id   <chr> "DSP.1001660011137.E.A05.dcc", "DSP.1001660011137.E.A05.dcc", "DSP.1001660011137.E.A05.dcc", "DSP.1001660011137.E.A05.dcc"…
# $ celltype <chr> "iCAF", "mCAF", "hsp_tpCAF", "IDO_CAF", "apCAF", "dCAF", "Pericyte", "tpCAF", "vCAF", "rCAF", "iCAF", "mCAF", "hsp_tpCAF",…
# $ prop     <dbl> 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.97319520, 0.02680480, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.…
# $ ncells   <dbl> 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 26…
# $ patient  <chr> "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1…
# $ group    <chr> "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ie…

mm <- model.matrix(~1+infiltration_type, data = fibmeta)
contrast_list <- list(
  "ieb / iib" = mm[fibmeta$infiltration_type == "ieb",] |> colMeans() - mm[fibmeta$infiltration_type == "iib",] |> colMeans(),
  "ieb / idb" = mm[fibmeta$infiltration_type == "ieb",] |> colMeans() - mm[fibmeta$infiltration_type == "idb",] |> colMeans(),
  "iib / idb" = mm[fibmeta$infiltration_type == "iib",] |> colMeans() - mm[fibmeta$infiltration_type == "idb",] |> colMeans()
)

proptestfunc <- function(ct, .contr) {
  tryCatch({
    
    testdata <- filter(fibdcvn_long, celltype == ct)
    modout <- glmer(formula = prop ~ 1+group + (1+group|patient), weights = testdata$ncells, data = testdata, family = binomial(link = "logit"))
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

res <- lapply(X = names(contrast_list), FUN = purrr::map, .x = unique(fibdcvn_long$celltype), .f = proptestfunc)
res <- bind_rows(res)
res <- group_by(res, contrast) |> mutate(fdr = p.adjust(p = pval, method = "BH"))
res$padj <- p.adjust(p = res$pval, method = "bonferroni")
knitr::kable(res, format = "simple")

# contrast    celltype          pval        Estimate        tstat          sigma         fdr        padj
# ----------  ----------  ----------  --------------  -----------  -------------  ----------  ----------
# ieb / iib   mCAF         0.7717131      -2.2113550   -0.2901348   7.621819e+00   0.9999792   1.0000000
# ieb / iib   hsp_tpCAF    0.6081341       3.0646941    0.5127388   5.977106e+00   0.9730145   1.0000000
# ieb / iib   IDO_CAF      0.2207077      -1.2710704   -1.2246484   1.037906e+00   0.4763473   1.0000000
# ieb / iib   apCAF        0.2381736      -0.9060961   -1.1795640   7.681618e-01   0.4763473   1.0000000
# ieb / iib   dCAF         0.0075665       2.2535751    2.6708260   8.437746e-01   0.0302659   0.1815957    ** No longer significant with bonferroni
# ieb / iib   Pericyte     0.9998708      25.2971311    0.0001619   1.562501e+05   0.9999792   1.0000000
# ieb / iib   tpCAF        0.0016235    2359.4177930    3.1516487   7.486297e+02   0.0129881   0.0389644    *
# ieb / iib   rCAF         0.9999792     -75.7953537   -0.0000261   2.904085e+06   0.9999792   1.0000000
# ieb / idb   mCAF         0.0761787      -9.4905110   -1.7733018   5.351887e+00   0.2031432   1.0000000
# ieb / idb   hsp_tpCAF    0.1809003      -2.9379132   -1.3379881   2.195769e+00   0.2894405   1.0000000
# ieb / idb   IDO_CAF      0.9660961      21.1380638    0.0425050   4.973077e+02   0.9999870   1.0000000
# ieb / idb   apCAF        0.0026653       1.7276969    3.0039119   5.751490e-01   0.0106613   0.0639678    ** No longer significant with bonferroni 
# ieb / idb   dCAF         0.0002895      -1.4759983   -3.6244814   4.072302e-01   0.0023163   0.0069490    *
# ieb / idb   Pericyte     0.9993012      20.6570656    0.0008758   2.358618e+04   0.9999870   1.0000000
# ieb / idb   tpCAF        0.1304966      -1.2363696   -1.5121464   8.176256e-01   0.2609933   1.0000000
# ieb / idb   rCAF         0.9999870     -51.2707680   -0.0000162   3.157991e+06   0.9999870   1.0000000
# iib / idb   mCAF         0.1668948      -7.2791560   -1.3822503   5.266163e+00   0.3337897   1.0000000
# iib / idb   hsp_tpCAF    0.2821058      -6.0026073   -1.0756009   5.580701e+00   0.4513692   1.0000000
# iib / idb   IDO_CAF      0.9640588      22.4091342    0.0450609   4.973076e+02   0.9999842   1.0000000
# iib / idb   apCAF        0.0000005       2.6337929    5.0312864   5.234830e-01   0.0000039   0.0000117    *
# iib / idb   dCAF         0.0000203      -3.7295734   -4.2616259   8.751527e-01   0.0000812   0.0004871    *
# iib / idb   Pericyte     0.9999766      -4.6400655   -0.0000294   1.580203e+05   0.9999842   1.0000000
# iib / idb   tpCAF        0.0016143   -2360.6541626   -3.1533019   7.486293e+02   0.0043049   0.0387443    *
# iib / idb   rCAF         0.9999842      24.5245857    0.0000198   1.240645e+06   0.9999842   1.0000000

# Many warnings, including thing about non-integer successes

# Same as above, with pseudocount of 0.001
################################################################################
proptestfunc <- function(ct, .contr) {
  tryCatch({
    
    testdata <- filter(fibdcvn_long, celltype == ct)
    modout <- glmer(formula = (prop+0.001) ~ 1+group + (1+group|patient), weights = testdata$ncells, data = testdata, family = binomial(link = "logit"))
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

res <- lapply(X = names(contrast_list), FUN = purrr::map, .x = unique(fibdcvn_long$celltype), .f = proptestfunc)
res <- bind_rows(res)
res <- group_by(res, contrast) |> mutate(fdr = p.adjust(p = pval, method = "BH"))
res$padj <- p.adjust(p = res$pval, method = "bonferroni")
knitr::kable(res, format = "simple")

# contrast    celltype          pval      Estimate        tstat          sigma         fdr        padj
# ----------  ----------  ----------  ------------  -----------  -------------  ----------  ----------
# ieb / iib   mCAF         0.8992031     0.5398740    0.1266681   4.262113e+00   0.9999923   1.0000000
# ieb / iib   hsp_tpCAF    0.8737636     0.5393552    0.1588798   3.394738e+00   0.9999923   1.0000000
# ieb / iib   IDO_CAF      0.2110834    -1.2717276   -1.2505920   1.016900e+00   0.8443335   1.0000000
# ieb / iib   dCAF         0.0081889     2.2516995    2.6441818   8.515676e-01   0.0655109   0.1965328
# ieb / iib   Pericyte     0.9999923    31.4396773    0.0000096   3.274577e+06   0.9999923   1.0000000
# ieb / iib   tpCAF        0.9326162    14.8654850    0.0845538   1.758111e+02   0.9999923   1.0000000
# ieb / iib   vCAF         0.9793480   -11.7093209   -0.0258864   4.523356e+02   0.9999923   1.0000000
# ieb / iib   rCAF         0.9973889   -17.0148068   -0.0032725   5.199301e+03   0.9999923   1.0000000
# ieb / idb   mCAF         0.1036849    -4.5230879   -1.6272459   2.779597e+00   0.2764930   1.0000000
# ieb / idb   hsp_tpCAF    0.0795411    -2.5955751   -1.7533549   1.480348e+00   0.2764930   1.0000000
# ieb / idb   IDO_CAF      0.9869307    19.1377531    0.0163807   1.168311e+03   0.9999978   1.0000000
# ieb / idb   dCAF         0.0002769    -1.4782104   -3.6359805   4.065507e-01   0.0022154   0.0066462    *
# ieb / idb   Pericyte     0.9999931    45.0035402    0.0000087   5.177561e+06   0.9999978   1.0000000
# ieb / idb   tpCAF        0.1485045    -1.2098374   -1.4448342   8.373538e-01   0.2970089   1.0000000
# ieb / idb   vCAF         0.9753064   -14.0013569   -0.0309538   4.523308e+02   0.9999978   1.0000000
# ieb / idb   rCAF         0.9999978    14.3704372    0.0000028   5.143237e+06   0.9999978   1.0000000
# iib / idb   mCAF         0.0668561    -5.0629619   -1.8326403   2.762660e+00   0.2674244   1.0000000
# iib / idb   hsp_tpCAF    0.3280719    -3.1349303   -0.9780050   3.205434e+00   0.6561437   1.0000000
# iib / idb   IDO_CAF      0.9860623    20.4094806    0.0174692   1.168311e+03   0.9999982   1.0000000
# iib / idb   dCAF         0.0000229    -3.7299098   -4.2349687   8.807408e-01   0.0001829   0.0005486    *
# iib / idb   Pericyte     0.9999982    13.5638629    0.0000022   6.126173e+06   0.9999982   1.0000000
# iib / idb   tpCAF        0.9271468   -16.0753224   -0.0914352   1.758111e+02   0.9999982   1.0000000
# iib / idb   vCAF         0.3224564    -2.2920361   -0.9894227   2.316539e+00   0.6561437   1.0000000
# iib / idb   rCAF         0.9999951    31.3852440    0.0000061   5.143235e+06   0.9999982   1.0000000

# Same warnings as above

# Binomial logistic regression with (rounded cell counts / n cells) as the response variable 
################################################################################

# I noticed here that the cell library sizes (total number of cells in an AOI) are different than the numbers I gave
# If the library sizes are not all scaled by the same factor, then the relative weights in the model are changed... right?
# Should I be using the 

mat <- fibdcvn$cell.counts$cell.counts |> round(digits = 0)
d <- mat |> t() |> as.data.frame()
d$aoi_id <- rownames(d)
fibdcvn_long <- tidyr::pivot_longer(data = d, 
                                    names_to = "celltype", 
                                    values_to = "count",
                                    cols = rownames(fibdcvn$cell.counts$cell.counts),
                                    values_drop_na = T)

fibdcvn_long <- group_by(fibdcvn_long, aoi_id) |> mutate(ncells = sum(count))
fibdcvn_long$patient <- plyr::mapvalues(x = fibdcvn_long$aoi_id, from = fibmeta$aoi_id, to = fibmeta$Patient_number)
fibdcvn_long$group <- plyr::mapvalues(x = fibdcvn_long$aoi_id, from = fibmeta$aoi_id, to = fibmeta$infiltration_type)

glimpse(fibdcvn_long)
# Rows: 280
# Columns: 6
# Groups: aoi_id [28]
# $ aoi_id   <chr> "DSP.1001660011137.E.A05.dcc", "DSP.1001660011137.E.A05.dcc", "DSP.1001660011137.E.A05.dcc", "DSP.1001660011137.E.A05.dcc"…
# $ celltype <chr> "iCAF", "mCAF", "hsp_tpCAF", "IDO_CAF", "apCAF", "dCAF", "Pericyte", "tpCAF", "vCAF", "rCAF", "iCAF", "mCAF", "hsp_tpCAF",…
# $ count    <dbl> 0, 0, 0, 0, 49, 1, 0, 0, 0, 0, 0, 0, 0, 0, 37, 3, 0, 0, 0, 0, 0, 0, 0, 0, 14, 3, 0, 0, 0, 0, 0, 0, 0, 0, 7, 2, 0, 0, 0, 0,…
# $ ncells   <dbl> 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 9,…
# $ patient  <chr> "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1…
# $ group    <chr> "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ie…

mm <- model.matrix(~1+infiltration_type, data = fibmeta)
contrast_list <- list(
  "ieb / iib" = mm[fibmeta$infiltration_type == "ieb",] |> colMeans() - mm[fibmeta$infiltration_type == "iib",] |> colMeans(),
  "ieb / idb" = mm[fibmeta$infiltration_type == "ieb",] |> colMeans() - mm[fibmeta$infiltration_type == "idb",] |> colMeans(),
  "iib / idb" = mm[fibmeta$infiltration_type == "iib",] |> colMeans() - mm[fibmeta$infiltration_type == "idb",] |> colMeans()
)

proptestfunc <- function(ct, .contr) {
  tryCatch({
    
    testdata <- filter(fibdcvn_long, celltype == ct)
    modout <- glmer(formula = (count/ncells) ~ 1+group + (1+group|patient), weights = testdata$ncells, data = testdata, family = binomial(link = "logit"))
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

res <- lapply(X = names(contrast_list), FUN = purrr::map, .x = unique(fibdcvn_long$celltype), .f = proptestfunc)
res <- bind_rows(res)
res <- group_by(res, contrast) |> mutate(fdr = p.adjust(p = pval, method = "BH"))
res$padj <- p.adjust(p = res$pval, method = "bonferroni")
knitr::kable(res, format = "simple")

# contrast    celltype          pval        Estimate        tstat          sigma         fdr        padj
# ----------  ----------  ----------  --------------  -----------  -------------  ----------  ----------
# ieb / iib   mCAF         0.9999638     -25.0488293   -0.0000453   5.525711e+05   0.9999999   1.0000000
# ieb / iib   hsp_tpCAF    0.8408544       1.1931612    0.2008007   5.942016e+00   0.9999999   1.0000000
# ieb / iib   IDO_CAF      0.5951685      -0.5111470   -0.5313611   9.619578e-01   0.9999999   1.0000000
# ieb / iib   apCAF        0.6910605      -0.2332319   -0.3974163   5.868705e-01   0.9999999   1.0000000
# ieb / iib   dCAF         0.3926804       4.3846385    0.8547665   5.129633e+00   0.9999999   1.0000000
# ieb / iib   Pericyte     0.9996816    1821.9796494    0.0003990   4.566180e+06   0.9999999   1.0000000
# ieb / iib   tpCAF        0.9999934      32.7006571    0.0000083   3.926668e+06   0.9999999   1.0000000
# ieb / iib   vCAF         0.9999999       1.1346794    0.0000002   6.239931e+06   0.9999999   1.0000000
# ieb / iib   rCAF         0.9999941     -31.5008082   -0.0000074   4.252851e+06   0.9999999   1.0000000
# ieb / idb   mCAF         0.9999570     -29.7630223   -0.0000539   5.525711e+05   0.9999993   1.0000000
# ieb / idb   hsp_tpCAF    0.1568168      -1.8811870   -1.4158590   1.328654e+00   0.3528378   1.0000000
# ieb / idb   IDO_CAF      0.9983519      19.1472509    0.0020656   9.269556e+03   0.9999993   1.0000000
# ieb / idb   apCAF        0.0010228       1.7166162    3.2841660   5.226947e-01   0.0046028   0.0276169    *
# ieb / idb   dCAF         0.0000499      -1.3562660   -4.0561447   3.343732e-01   0.0004490   0.0013470    *
# ieb / idb   Pericyte     0.9995059      20.3539029    0.0006192   3.287107e+04   0.9999993   1.0000000
# ieb / idb   tpCAF        0.0960229      -1.2949339   -1.6644484   7.779958e-01   0.2880686   1.0000000
# ieb / idb   vCAF         0.9999933     -35.4965390   -0.0000083   4.252851e+06   0.9999993   1.0000000
# ieb / idb   rCAF         0.9999993      -4.0223627   -0.0000009   4.696147e+06   0.9999993   1.0000000
# iib / idb   mCAF         0.0439321      -4.7141930   -2.0147384   2.339854e+00   0.1976943   1.0000000
# iib / idb   hsp_tpCAF    0.6012603      -3.0743482   -0.5225889   5.882919e+00   0.9999936   1.0000000
# iib / idb   IDO_CAF      0.9983079      19.6583979    0.0021207   9.269556e+03   0.9999936   1.0000000
# iib / idb   apCAF        0.0000001       1.9498481    5.3342022   3.655370e-01   0.0000009   0.0000026    *
# iib / idb   dCAF         0.2634477      -5.7409045   -1.1182794   5.133694e+00   0.7903430   1.0000000
# iib / idb   Pericyte     0.9996852   -1801.6257464   -0.0003945   4.566298e+06   0.9999936   1.0000000
# iib / idb   tpCAF        0.9999931     -33.9955910   -0.0000087   3.926668e+06   0.9999936   1.0000000
# iib / idb   vCAF         0.9999936     -36.6312184   -0.0000080   4.566180e+06   0.9999936   1.0000000
# iib / idb   rCAF         0.9999890      27.4784455    0.0000138   1.991747e+06   0.9999936   1.0000000

# Still warnings, but the things about the non-integer successes is gone

# Same as above, with pseudocount of 0.001 
################################################################################
proptestfunc <- function(ct, .contr) {
  tryCatch({
    
    testdata <- filter(fibdcvn_long, celltype == ct)
    modout <- glmer(formula = ((count/ncells)+0.001) ~ 1+group + (1+group|patient), weights = testdata$ncells, data = testdata, family = binomial(link = "logit"))
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

res <- lapply(X = names(contrast_list), FUN = purrr::map, .x = unique(fibdcvn_long$celltype), .f = proptestfunc)
res <- bind_rows(res)
res <- group_by(res, contrast) |> mutate(fdr = p.adjust(p = pval, method = "BH"))
res$padj <- p.adjust(p = res$pval, method = "bonferroni")
knitr::kable(res, format = "simple")

# contrast    celltype          pval      Estimate           tstat          sigma         fdr        padj
# ----------  ----------  ----------  ------------  --------------  -------------  ----------  ----------
# ieb / iib   mCAF         0.9044706    -7.8172496   -1.200159e-01   6.513514e+01   0.9999626   1.0000000
# ieb / iib   hsp_tpCAF    0.9501714     0.2605053    6.249150e-02   4.168648e+00   0.9999626   1.0000000 
# ieb / iib   dCAF         0.1766632     3.3540698    1.351101e+00   2.482472e+00   0.5299897   1.0000000
# ieb / iib   Pericyte     0.9999626   214.1481335    4.690000e-05   4.566180e+06   0.9999626   1.0000000
# ieb / iib   vCAF         0.9957738    -0.7822303   -5.296800e-03   1.476791e+02   0.9999626   1.0000000
# ieb / iib   rCAF         0.0000000   -94.3253299   -1.502222e+04   6.279100e-03   0.0000000   0.0000000   *
# ieb / idb   mCAF         0.8569737   -11.7329351   -1.802278e-01   6.510058e+01   0.9999732   1.0000000
# ieb / idb   hsp_tpCAF    0.1411097    -1.8719927   -1.471671e+00   1.272018e+00   0.2822194   1.0000000
# ieb / idb   dCAF         0.0000534    -1.3555220   -4.040087e+00   3.355181e-01   0.0001603   0.0009618   *
# ieb / idb   Pericyte     0.9999732   251.7774534    3.360000e-05   7.502999e+06   0.9999732   1.0000000
# ieb / idb   vCAF         0.9285691   -10.7376048   -8.964530e-02   1.197788e+02   0.9999732   1.0000000
# ieb / idb   rCAF         0.0000000   -65.3008232   -1.470841e+04   4.439700e-03   0.0000000   0.0000000   *
# iib / idb   mCAF         0.0760611    -3.9156856   -1.774012e+00   2.207248e+00   0.1521222   1.0000000
# iib / idb   hsp_tpCAF    0.6028598    -2.1324980   -5.202924e-01   4.098653e+00   0.9042896   1.0000000
# iib / idb   dCAF         0.0584683    -4.7095918   -1.892170e+00   2.488989e+00   0.1521222   1.0000000
# iib / idb   Pericyte     0.9999966    37.6293199    4.300000e-06   8.783222e+06   0.9999966   1.0000000
# iib / idb   vCAF         0.9082624    -9.9553745   -1.152305e-01   8.639528e+01   0.9999966   1.0000000
# iib / idb   rCAF         0.0000000    29.0245067    6.536713e+03   4.440200e-03   0.0000000   0.0000000   *

# Warnings about the non-integer successes comes back

# Linear regression on log2(CPH + 0.01)
################################################################################
mat <- fibdcvn$cell.counts$cells.per.100
d <- mat |> t() |> as.data.frame()
d$aoi_id <- rownames(d)
fibdcvn_long <- tidyr::pivot_longer(data = d, 
                                    names_to = "celltype", 
                                    values_to = "cph",
                                    cols = rownames(fibdcvn$cell.counts$cells.per.100),
                                    values_drop_na = T)

fibdcvn_long <- mutate(fibdcvn_long, y = log2(cph + 0.01))
fibdcvn_long$patient <- plyr::mapvalues(x = fibdcvn_long$aoi_id, from = fibmeta$aoi_id, to = fibmeta$Patient_number)
fibdcvn_long$group <- plyr::mapvalues(x = fibdcvn_long$aoi_id, from = fibmeta$aoi_id, to = fibmeta$infiltration_type)

glimpse(fibdcvn_long)
# Rows: 280
# Columns: 6
# $ aoi_id   <chr> "DSP.1001660011137.E.A05.dcc", "DSP.1001660011137.E.A05.dcc", "DSP.1001660011137.E.A05.dcc", "DSP.1001660011137.E.A05.dcc"…
# $ celltype <chr> "iCAF", "mCAF", "hsp_tpCAF", "IDO_CAF", "apCAF", "dCAF", "Pericyte", "tpCAF", "vCAF", "rCAF", "iCAF", "mCAF", "hsp_tpCAF",…
# $ cph      <dbl> 0.000000, 0.000000, 0.000000, 0.000000, 58.507282, 1.611471, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0…
# $ y        <dbl> -6.6438562, -6.6438562, -6.6438562, -6.6438562, 5.8707909, 0.6973036, -6.6438562, -6.6438562, -6.6438562, -6.6438562, -6.6…
# $ patient  <chr> "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1…
# $ group    <chr> "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ie…

mm <- model.matrix(~1+infiltration_type, data = fibmeta)
contrast_list <- list(
  "ieb / iib" = mm[fibmeta$infiltration_type == "ieb",] |> colMeans() - mm[fibmeta$infiltration_type == "iib",] |> colMeans(),
  "ieb / idb" = mm[fibmeta$infiltration_type == "ieb",] |> colMeans() - mm[fibmeta$infiltration_type == "idb",] |> colMeans(),
  "iib / idb" = mm[fibmeta$infiltration_type == "iib",] |> colMeans() - mm[fibmeta$infiltration_type == "idb",] |> colMeans()
)

proptestfunc <- function(ct, .contr) {
  tryCatch({
    
    testdata <- filter(fibdcvn_long, celltype == ct)
    modout <- lme4::lmer(formula = y ~ 1+group + (1+group|patient), data = testdata, REML = F)
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

res <- lapply(X = names(contrast_list), FUN = purrr::map, .x = unique(fibdcvn_long$celltype), .f = proptestfunc)
res <- bind_rows(res)
res <- group_by(res, contrast) |> mutate(fdr = p.adjust(p = pval, method = "BH"))
res$padj <- p.adjust(p = res$pval, method = "bonferroni")
knitr::kable(res, format = "simple")

# contrast    celltype          pval     Estimate        tstat       sigma         fdr        padj
# ----------  ----------  ----------  -----------  -----------  ----------  ----------  ----------
# ieb / iib   iCAF         0.8541364    0.0000000    0.1838433   0.0000000   0.8541364   1.0000000
# ieb / iib   mCAF         0.7595930   -0.5060553   -0.3060153   1.6536927   0.8439922   1.0000000
# ieb / iib   hsp_tpCAF    0.6758058   -0.4872647   -0.4181933   1.1651661   0.8439922   1.0000000
# ieb / iib   IDO_CAF      0.0913556   -2.9815782   -1.6882901   1.7660344   0.2950531   1.0000000
# ieb / iib   apCAF        0.1180212   -2.1955559   -1.5631334   1.4045864   0.2950531   1.0000000
# ieb / iib   dCAF         0.0345193    2.9315182    2.1139525   1.3867474   0.2950531   1.0000000
# ieb / iib   Pericyte     0.3434702    0.6485523    0.9473309   0.6846100   0.4906717   1.0000000
# ieb / iib   tpCAF        0.0635150    2.4898859    1.8555690   1.3418450   0.2950531   1.0000000
# ieb / iib   vCAF         0.2341579   -1.2106713   -1.1897162   1.0176135   0.4504760   1.0000000
# ieb / iib   rCAF         0.2702856   -1.3350709   -1.1024051   1.2110529   0.4504760   1.0000000
# ieb / idb   iCAF         0.0078103    0.0000000    2.6601642   0.0000000   0.0264924   0.2343077    ** No longer significant with bonferroni
# ieb / idb   mCAF         0.0001465   -4.9963992   -3.7968691   1.3159261   0.0014654   0.0043961    *
# ieb / idb   hsp_tpCAF    0.3334713   -1.5430307   -0.9671454   1.5954485   0.4763876   1.0000000
# ieb / idb   IDO_CAF      0.1491250    2.9307807    1.4426290   2.0315554   0.2485417   1.0000000
# ieb / idb   apCAF        0.4843972   -1.1004226   -0.6992477   1.5737236   0.5433804   1.0000000
# ieb / idb   dCAF         0.0583101   -3.2566093   -1.8933596   1.7200162   0.1392090   1.0000000
# ieb / idb   Pericyte     0.4890424    0.6485523    0.6918327   0.9374409   0.5433804   1.0000000
# ieb / idb   tpCAF        0.0079477   -4.8769705   -2.6542825   1.8373969   0.0264924   0.2384317    ** No longer significant with bonferroni
# ieb / idb   vCAF         0.0696045   -1.7198697   -1.8144758   0.9478603   0.1392090   1.0000000
# ieb / idb   rCAF         0.9999941   -0.0000058   -0.0000073   0.7955004   0.9999941   1.0000000
# iib / idb   iCAF         0.0155944    0.0000000    2.4182729   0.0000000   0.0311888   0.4678313    ** No longer significant with bonferroni
# iib / idb   mCAF         0.0052197   -4.4903439   -2.7931548   1.6076244   0.0165565   0.1565901    ** No longer significant with bonferroni
# iib / idb   hsp_tpCAF    0.5259845   -1.0557660   -0.6341477   1.6648584   0.6574806   1.0000000
# iib / idb   IDO_CAF      0.0066226    5.9123589    2.7152488   2.1774649   0.0165565   0.1986780    ** No longer significant with bonferroni
# iib / idb   apCAF        0.3254219    1.0951333    0.9833770   1.1136454   0.4648885   1.0000000
# iib / idb   dCAF         0.0002194   -6.1881275   -3.6955094   1.6744992   0.0010972   0.0065834    *
# iib / idb   Pericyte     1.0000000    0.0000000    0.0000000   0.9782181   1.0000000   1.0000000
# iib / idb   tpCAF        0.0001219   -7.3668564   -3.8422661   1.9173207   0.0010972   0.0036571    *
# iib / idb   vCAF         0.6765932   -0.5091983   -0.4171165   1.2207579   0.7517702   1.0000000
# iib / idb   rCAF         0.3221230    1.3350651    0.9901046   1.3484082   0.4648885   1.0000000

# Still many warnings

# Poisson regression on rounded counts
################################################################################
mat <- fibdcvn$cell.counts$cell.counts |> round(digits = 0)
d <- mat |> t() |> as.data.frame()
d$aoi_id <- rownames(d)
fibdcvn_long <- tidyr::pivot_longer(data = d, 
                                    names_to = "celltype", 
                                    values_to = "count",
                                    cols = rownames(fibdcvn$cell.counts$cell.counts),
                                    values_drop_na = T)

fibdcvn_long <- group_by(fibdcvn_long, aoi_id) |> mutate(ncells = sum(count))
fibdcvn_long$patient <- plyr::mapvalues(x = fibdcvn_long$aoi_id, from = fibmeta$aoi_id, to = fibmeta$Patient_number)
fibdcvn_long$group <- plyr::mapvalues(x = fibdcvn_long$aoi_id, from = fibmeta$aoi_id, to = fibmeta$infiltration_type)

glimpse(fibdcvn_long)
# Rows: 280
# Columns: 6
# Groups: aoi_id [28]
# $ aoi_id   <chr> "DSP.1001660011137.E.A05.dcc", "DSP.1001660011137.E.A05.dcc", "DSP.1001660011137.E.A05.dcc", "DSP.1001660011137.E.A05.dcc"…
# $ celltype <chr> "iCAF", "mCAF", "hsp_tpCAF", "IDO_CAF", "apCAF", "dCAF", "Pericyte", "tpCAF", "vCAF", "rCAF", "iCAF", "mCAF", "hsp_tpCAF",…
# $ count    <dbl> 0, 0, 0, 0, 49, 1, 0, 0, 0, 0, 0, 0, 0, 0, 37, 3, 0, 0, 0, 0, 0, 0, 0, 0, 14, 3, 0, 0, 0, 0, 0, 0, 0, 0, 7, 2, 0, 0, 0, 0,…
# $ ncells   <dbl> 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 9,…
# $ patient  <chr> "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1…
# $ group    <chr> "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ieb", "ie…

mm <- model.matrix(~1+infiltration_type, data = fibmeta)
contrast_list <- list(
  "ieb / iib" = mm[fibmeta$infiltration_type == "ieb",] |> colMeans() - mm[fibmeta$infiltration_type == "iib",] |> colMeans(),
  "ieb / idb" = mm[fibmeta$infiltration_type == "ieb",] |> colMeans() - mm[fibmeta$infiltration_type == "idb",] |> colMeans(),
  "iib / idb" = mm[fibmeta$infiltration_type == "iib",] |> colMeans() - mm[fibmeta$infiltration_type == "idb",] |> colMeans()
)

proptestfunc <- function(ct, .contr) {
  tryCatch({
    
    testdata <- filter(fibdcvn_long, celltype == ct)
    modout <- glmer(formula = count ~ 1+group + (1+group|patient) + offset(log(ncells)), data = testdata, family = poisson())
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

res <- lapply(X = names(contrast_list), FUN = purrr::map, .x = unique(fibdcvn_long$celltype), .f = proptestfunc)
res <- bind_rows(res)
res <- group_by(res, contrast) |> mutate(fdr = p.adjust(p = pval, method = "BH"))
res$padj <- p.adjust(p = res$pval, method = "bonferroni")
knitr::kable(res, format = "simple")

# contrast    celltype          pval        Estimate         tstat          sigma         fdr        padj
# ----------  ----------  ----------  --------------  ------------  -------------  ----------  ----------
# ieb / iib   mCAF         0.0000000   -1.451900e+05   -50.1291747   2.896318e+03   0.0000000   0.0000000   *
# ieb / iib   hsp_tpCAF    0.8948527    7.902830e-01     0.1321663   5.979458e+00   0.9999999   1.0000000
# ieb / iib   IDO_CAF      0.6204773   -5.044570e-01    -0.4951740   1.018747e+00   0.9999999   1.0000000
# ieb / iib   apCAF        0.8265060    2.265390e-02     0.2191849   1.033552e-01   0.9999999   1.0000000
# ieb / iib   dCAF         0.3814988    4.096178e+00     0.8751380   4.680608e+00   0.9999999   1.0000000
# ieb / iib   Pericyte     0.9995298    1.250480e+04     0.0005892   2.122169e+07   0.9999999   1.0000000
# ieb / iib   tpCAF        0.9999980    3.522439e+01     0.0000025   1.393921e+07   0.9999999   1.0000000
# ieb / iib   vCAF         0.9999999    2.445388e+00     0.0000001   2.432150e+07   0.9999999   1.0000000
# ieb / iib   rCAF         0.9999978   -4.975808e+01    -0.0000028   1.793560e+07   0.9999999   1.0000000
# ieb / idb   mCAF         0.0000000   -1.451943e+05   -50.1307681   2.896312e+03   0.0000000   0.0000000   *
# ieb / idb   hsp_tpCAF    0.1580183   -1.858502e+00    -1.4117678   1.316436e+00   0.2844330   1.0000000
# ieb / idb   IDO_CAF      0.9655318    1.987618e+01     0.0432130   4.599587e+02   0.9999991   1.0000000
# ieb / idb   apCAF        0.0002117    7.043790e-01     3.7045930   1.901367e-01   0.0006352   0.0057167   *
# ieb / idb   dCAF         0.0001600   -1.103301e+00    -3.7750781   2.922590e-01   0.0006352   0.0043189   *
# ieb / idb   Pericyte     0.9999765    2.699893e+01     0.0000295   9.151871e+05   0.9999991   1.0000000
# ieb / idb   tpCAF        0.0778886   -1.247546e+00    -1.7630704   7.075988e-01   0.1752494   1.0000000
# ieb / idb   vCAF         0.9999977   -3.381577e+01    -0.0000028   1.188172e+07   0.9999991   1.0000000
# ieb / idb   rCAF         0.9999991   -2.074162e+01    -0.0000011   1.840067e+07   0.9999991   1.0000000
# iib / idb   mCAF         0.2448721   -4.316959e+00    -1.1628951   3.712251e+00   0.8005540   1.0000000
# iib / idb   hsp_tpCAF    0.6541057   -2.648785e+00    -0.4480657   5.911598e+00   0.9999986   1.0000000
# iib / idb   IDO_CAF      0.9646575    2.038064e+01     0.0443097   4.599587e+02   0.9999986   1.0000000
# iib / idb   apCAF        0.0003961    6.817251e-01     3.5426677   1.924327e-01   0.0035649   0.0106947   *
# iib / idb   dCAF         0.2668513   -5.199478e+00    -1.1103428   4.682769e+00   0.8005540   1.0000000
# iib / idb   Pericyte     0.9995313   -1.247780e+04    -0.0005874   2.124141e+07   0.9999986   1.0000000
# iib / idb   tpCAF        0.9999979   -3.647194e+01    -0.0000026   1.393921e+07   0.9999986   1.0000000
# iib / idb   vCAF         0.9999986   -3.626116e+01    -0.0000017   2.122169e+07   0.9999986   1.0000000
# iib / idb   rCAF         0.9999944    2.901646e+01     0.0000071   4.110847e+06   0.9999986   1.0000000

# Still some warnings

# Negative binomial regression on rounded counts
################################################################################

# Only use this if there is overdispersion? 

proptestfunc <- function(ct, .contr) {
  tryCatch({
    
    testdata <- filter(fibdcvn_long, celltype == ct)
    modout <- glmer.nb(formula = count ~ 1+group + (1+group|patient) + offset(log(ncells)), data = testdata)
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

res <- lapply(X = names(contrast_list), FUN = purrr::map, .x = unique(fibdcvn_long$celltype), .f = proptestfunc)
res <- bind_rows(res)
res <- group_by(res, contrast) |> mutate(fdr = p.adjust(p = pval, method = "BH"))
res$padj <- p.adjust(p = res$pval, method = "bonferroni")
knitr::kable(res, format = "simple")

# contrast    celltype          pval        Estimate           tstat          sigma         fdr        padj
# ----------  ----------  ----------  --------------  --------------  -------------  ----------  ----------
# ieb / iib   mCAF         0.0000000     -36.4693704   -4806.3501785   7.587700e-03   0.0000000   0.0000000   *
# ieb / iib   hsp_tpCAF    0.2438767      -1.7865883      -1.1653515   1.533090e+00   0.6503379   1.0000000
# ieb / iib   IDO_CAF      0.6622370      -0.4539175      -0.4368267   1.039125e+00   0.9999999   1.0000000
# ieb / iib   apCAF        0.8263590       0.0226480       0.2193737   1.032394e-01   0.9999999   1.0000000
# ieb / iib   dCAF         0.3724902       4.0671557       0.8918185   4.560520e+00   0.7449805   1.0000000
# ieb / iib   tpCAF        0.0001329    4205.8598059       3.8210094   1.100720e+03   0.0005316   0.0031898   *
# ieb / iib   vCAF         0.9999999       3.1306167       0.0000001   2.152223e+07   0.9999999   1.0000000
# ieb / iib   rCAF         0.9750600     -87.8433288      -0.0312628   2.809838e+03   0.9999999   1.0000000
# ieb / idb   mCAF         0.0000000     -40.7985480   -7604.0919064   5.365300e-03   0.0000000   0.0000000   *
# ieb / idb   hsp_tpCAF    0.3791521      -1.6826146      -0.8794597   1.913237e+00   0.6066434   1.0000000
# ieb / idb   IDO_CAF      0.9992064      20.8377357       0.0009947   2.094945e+04   0.9999960   1.0000000
# ieb / idb   apCAF        0.0001478       0.7043715       3.7947922   1.856153e-01   0.0005911   0.0035464   *
# ieb / idb   dCAF         0.0010546      -1.0809996      -3.2755306   3.300227e-01   0.0028124   0.0253113   *
# ieb / idb   tpCAF        0.1889349      -1.3304896      -1.3137363   1.012752e+00   0.3778699   1.0000000
# ieb / idb   vCAF         0.9999960     -32.5480686      -0.0000050   6.546608e+06   0.9999960   1.0000000
# ieb / idb   rCAF         0.9805635     -68.4545997      -0.0243625   2.809836e+03   0.9999960   1.0000000
# iib / idb   mCAF         0.0000000      -4.3291775    -806.9375923   5.364900e-03   0.0000000   0.0000000   *
# iib / idb   hsp_tpCAF    0.9501988       0.1039737       0.0624571   1.664722e+00   0.9999986   1.0000000
# iib / idb   IDO_CAF      0.9991891      21.2916532       0.0010163   2.094945e+04   0.9999986   1.0000000
# iib / idb   apCAF        0.0002866       0.6817235       3.6271245   1.879515e-01   0.0005732   0.0068783   *
# iib / idb   dCAF         0.2591966      -5.1481553      -1.1282921   4.562786e+00   0.4147145   1.0000000
# iib / idb   tpCAF        0.0001323   -4207.1902955      -3.8222202   1.100719e+03   0.0003527   0.0031741   *
# iib / idb   vCAF         0.9999986     -35.6786853      -0.0000017   2.050240e+07   0.9999986   1.0000000
# iib / idb   rCAF         0.0000004      19.3887292       5.0545372   3.835906e+00   0.0000017   0.0000104   *

# Same warnings as above

# Notes 
################################################################################

# Looking at the bar plot, the original model and the Poisson model seem to make
# the most sense, by eye.

# Am I choosing pseudo-count correctly? 

# The pseudocount does seem to limit the estimates. Is this a good sign? 


