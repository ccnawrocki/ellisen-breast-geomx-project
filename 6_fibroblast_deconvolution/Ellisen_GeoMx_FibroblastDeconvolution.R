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

# When I give the spatialdecon function counts of nuclei, it estimates 1 or less counts at times.
# Also, counts per 100 that is output is not accurate; columns do not add to 100.
# Finally, giving the cell numbers does not change anything about the proportions or the the betas 
# that are predicted...
# Therefore, I will use proportion that it output by the function and actual nuclei  counts to make 
# discrete counts of cells. Then, I will model on that.


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
                        raw = fibcts 
                        #cell_counts = fibmeta$Number_of_cells_sampled
                        )

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

fibdcvn_long$patient <- plyr::mapvalues(x = fibdcvn_long$aoi_id, from = fibmeta$aoi_id, to = fibmeta$Patient_number) |> as.factor()
aoiorder <- fibmeta |> arrange(Patient_number) |> pull(aoi_id)
fibdcvn_long$aoi_id <- factor(fibdcvn_long$aoi_id, levels = aoiorder)

p2 <- ggplot(fibdcvn_long) + 
  geom_bar(mapping = aes(x = aoi_id, y = prop, fill = celltype, group = aoi_id), stat = "identity", color = NA) + 
  ggprism::scale_fill_prism() + 
  theme_dendro() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(title = "Fibroblast, Proportions")

p2 + geom_text(data=fibmeta[aoiorder,],
               aes(label=infiltration_type, x=aoi_id, y=1.075, color=infiltration_type)) +
  scale_color_manual(values = c("ieb"="dodgerblue", "iib"="purple", "idb"="orange")) + 
  geom_text(data=fibmeta[aoiorder,],
            aes(label=Patient_number, x=aoi_id, y=1.175)) 

aoitree <- hclust(dist(fibdcvn$prop_of_all |> t()))
dend <- aoitree |> as.dendrogram()

ddata_x <- dendro_data(dend)
ddata_x$labels$patient <- plyr::mapvalues(x = ddata_x$labels$label, from = fibmeta$aoi_id, to = fibmeta$Patient_number)
aoiorder <- aoitree$labels[aoitree$order]

p1 <- ggplot(fibdcvn_long) + 
  geom_bar(mapping = aes(x = aoi_id, y = prop, fill = celltype, group = aoi_id), stat = "identity", color = NA) + 
  ggprism::scale_fill_prism() + 
  theme_dendro() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_x_dendrogram(hclust = aoitree, position = "top", labels = NULL)
p1 + geom_point(data=ddata_x$labels, inherit.aes = F,
               aes(color = patient, x = x, y = 1.05), size = 8, shape = 15) + 
  ggprism::scale_color_prism()


## MODELING --------------------------------------------------------------------
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

    # Binomial logistic regression, random slopes
    #modout <- glmer(formula = (count/ncells) ~ 1+group + (1+group|patient), weights = testdata$ncells, data = testdata, family = binomial(link = "logit"))
    
    # Binomial logistic regression, random intercepts
    #modout <- glmer(formula = (count/ncells) ~ 1+group + (1+1|patient), weights = testdata$ncells, data = testdata, family = binomial(link = "logit"))
    
    # Binomial logistic regression, patient as fixed effect (use other model matrix)
    #modout <- glm(formula = (count/ncells) ~ 1+group+patient, weights = testdata$ncells, data = testdata, family = binomial(link = "logit"))
    
    # Binomial logistic regression, patient excluded from model
    #modout <- glm(formula = (count/ncells) ~ 1+group, weights = testdata$ncells, data = testdata, family = binomial(link = "logit"))
    
    # Poisson regression, random slopes
    #modout <- glmer(formula = count ~ 1+group+(1+group|patient)+offset(log(ncells)), data = testdata, family = poisson(link = "log"))
    
    # Poisson regression, random intercepts
    #modout <- glmer(formula = count ~ 1+group+(1+1|patient)+offset(log(ncells)), data = testdata, family = poisson(link = "log"))
    
    # Poisson regression, patient as fixed effect (use other model matrix)
    #modout <- glm(formula = count ~ 1+group+patient+offset(log(ncells)), data = testdata, family = poisson(link = "log"))
    
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
levels(fibmeta$infiltration_type) <- levels(fibdcvn_long$group)
all(levels(fibmeta$infiltration_type) == levels(fibdcvn_long$group))
# [1] TRUE

fibmeta$Patient_number %<>% as.factor()
levels(fibmeta$Patient_number) <- levels(fibdcvn_long$patient)
all(levels(fibmeta$Patient_number) == levels(fibdcvn_long$patient))
# [1] TRUE

mm <- model.matrix(~1+infiltration_type, data = fibmeta)
#mm <- model.matrix(~1+infiltration_type+Patient_number, data = fibmeta)
contrast_list <- list(
  "ieb / iib" = mm[fibmeta$infiltration_type == "ieb",] |> colMeans() - mm[fibmeta$infiltration_type == "iib",] |> colMeans(),
  "ieb / idb" = mm[fibmeta$infiltration_type == "ieb",] |> colMeans() - mm[fibmeta$infiltration_type == "idb",] |> colMeans(),
  "iib / idb" = mm[fibmeta$infiltration_type == "iib",] |> colMeans() - mm[fibmeta$infiltration_type == "idb",] |> colMeans()
)

# Note that there are no warning messages when using the simplest poisson model.
res <- lapply(X = names(contrast_list), FUN = purrr::map, .x = unique(fibdcvn_long$celltype), .f = proptestfunc)
res <- bind_rows(res)
res <- group_by(res, contrast) |> mutate(fdr = p.adjust(p = pval, method = "BH"))
res$padj <- p.adjust(p = res$pval, method = "bonferroni")
knitr::kable(res, format = "simple")

# Feeling better about these results. The estimates are relatively controlled.

# contrast    celltype          pval     Estimate        tstat          sigma         fdr        padj
# ----------  ----------  ----------  -----------  -----------  -------------  ----------  ----------
# ieb / iib   iCAF         0.9980584   18.6482999    0.0024334   7.663473e+03   0.9999993   1.0000000
# ieb / iib   mCAF         0.0573043    0.4736253    1.9009825   2.491477e-01   0.1146086   1.0000000
# ieb / iib   hsp_tpCAF    0.9999993    0.1123544    0.0000009   1.218317e+05   0.9999993   1.0000000
# ieb / iib   IDO_CAF      0.0000409   -1.2052220   -4.1023971   2.937848e-01   0.0004089   0.0012267
# ieb / iib   apCAF        0.1278427   -0.1177830   -1.5226642   7.735330e-02   0.2130712   1.0000000
# ieb / iib   dCAF         0.0002592    0.8450481    3.6530260   2.313283e-01   0.0012958   0.0077750
# ieb / iib   Pericyte     0.9967820   18.7469122    0.0040332   4.648131e+03   0.9999993   1.0000000
# ieb / iib   tpCAF        0.0005297    2.0952338    3.4652496   6.046415e-01   0.0017658   0.0158922
# ieb / iib   vCAF         0.9980584   18.6482999    0.0024334   7.663473e+03   0.9999993   1.0000000
# ieb / iib   rCAF         0.0153338   -2.5427262   -2.4243968   1.048808e+00   0.0383346   0.4600153
# ieb / idb   iCAF         0.9988105   18.7843090    0.0014908   1.260011e+04   0.9999988   1.0000000
# ieb / idb   mCAF         0.0000000   -1.3811793   -7.3005013   1.891896e-01   0.0000000   0.0000000
# ieb / idb   hsp_tpCAF    0.9999988    0.2483635    0.0000014   1.723059e+05   0.9999988   1.0000000
# ieb / idb   IDO_CAF      0.9872279   16.5568977    0.0160081   1.034280e+03   0.9999988   1.0000000
# ieb / idb   apCAF        0.0000000    1.1927995    6.8203915   1.748872e-01   0.0000000   0.0000000
# ieb / idb   dCAF         0.0000067   -0.8018804   -4.5041413   1.780318e-01   0.0000222   0.0001999
# ieb / idb   Pericyte     0.9980286   18.8829213    0.0024708   7.642352e+03   0.9999988   1.0000000
# ieb / idb   tpCAF        0.0121024   -0.7091475   -2.5091435   2.826253e-01   0.0302561   0.3630729
# ieb / idb   vCAF         0.9988105   18.7843090    0.0014908   1.260011e+04   0.9999988   1.0000000
# ieb / idb   rCAF         0.4372549   -1.0986123   -0.7768370   1.414212e+00   0.8745099   1.0000000
# iib / idb   iCAF         0.9999926    0.1360091    0.0000092   1.474759e+04   0.9999994   1.0000000
# iib / idb   mCAF         0.0000000   -1.8548046   -7.7654442   2.388536e-01   0.0000000   0.0000000
# iib / idb   hsp_tpCAF    0.9999994    0.1360091    0.0000008   1.796625e+05   0.9999994   1.0000000
# iib / idb   IDO_CAF      0.9862983   17.7621197    0.0171734   1.034280e+03   0.9999994   1.0000000
# iib / idb   apCAF        0.0000000    1.3105825    7.4494707   1.759296e-01   0.0000000   0.0000000
# iib / idb   dCAF         0.0000000   -1.6469285   -6.8278485   2.412075e-01   0.0000000   0.0000000
# iib / idb   Pericyte     0.9999879    0.1360091    0.0000152   8.944868e+03   0.9999994   1.0000000
# iib / idb   tpCAF        0.0000055   -2.8043813   -4.5436170   6.172134e-01   0.0000138   0.0001659
# iib / idb   vCAF         0.9999926    0.1360091    0.0000092   1.474759e+04   0.9999994   1.0000000
# iib / idb   rCAF         0.1685401    1.4441139    1.3769099   1.048808e+00   0.3370802   1.0000000

# I think using difference in average proportion will work for visualizations.
mm <- model.matrix(~0+infiltration_type, data = fibmeta)
contrast_list <- list(
  "ieb / iib" = mm[fibmeta$infiltration_type == "ieb",] |> colMeans() - mm[fibmeta$infiltration_type == "iib",] |> colMeans(),
  "ieb / idb" = mm[fibmeta$infiltration_type == "ieb",] |> colMeans() - mm[fibmeta$infiltration_type == "idb",] |> colMeans(),
  "iib / idb" = mm[fibmeta$infiltration_type == "iib",] |> colMeans() - mm[fibmeta$infiltration_type == "idb",] |> colMeans()
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

# Using limma... still getting warning messages
# Have to fit using the abundance score, not the proportion
mat <- (fibdcvn$beta+1) |> log2()
mm <- model.matrix(~0+infiltration_type, data = fibmeta)
fit <- limma::lmFit(mat, design = mm)
deres <- data.frame()
for (c in names(contrast_list)) {
  deout <- contrasts.fit(fit, contrasts = contrast_list[[c]])
  deout <- eBayes(deout) |> topTable(number = Inf)
  deout$celltype <- rownames(deout)
  deout$contrast <- c
  deres <- rbind(deres, deout)
}

ggplot(data = deres) + 
  geom_bar(mapping = aes(y = celltype, x = logFC, fill = adj.P.Val < 0.05), stat = "identity") + 
  facet_grid(.~contrast) +
  scale_fill_manual(values = c("blue", "red")) + 
  ggthemes::theme_few() + 
  geom_vline(xintercept = 0, linewidth = 0.2)

dev.off()

write.csv(res, "fibroblast_deconvolution_differential_abundance_testing_results.csv")
write.csv(deres, "fibroblast_deconvolution_differential_abundance_testing_results_limma.csv")

# Session 

# sessionInfo()
# 
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
# [96] htmltools_0.5.8.1       pkgconfig_2.0.3         umap_0.2.10.0           highr_0.11              fastmap_1.2.0          
# [101] rlang_1.1.4             GlobalOptions_0.1.2     htmlwidgets_1.6.4       ggthemes_5.1.0          UCSC.utils_1.0.0       
# [106] shiny_1.9.1             farver_2.1.2            zoo_1.8-12              jsonlite_1.8.9          R.oo_1.26.0            
# [111] Formula_1.2-5           GenomeInfoDbData_1.2.12 dotCall64_1.1-1         munsell_0.5.1           Rcpp_1.0.13            
# [116] reticulate_1.39.0       stringi_1.8.4           zlibbioc_1.50.0         plyr_1.8.9              ggstats_0.7.0          
# [121] listenv_0.9.1           ggrepel_0.9.6           deldir_2.0-4            Biostrings_2.72.1       splines_4.4.1          
# [126] tensor_1.5              circlize_0.4.16         locfit_1.5-9.6          igraph_2.0.3            uuid_1.2-1             
# [131] spatstat.geom_3.3-5     RcppHNSW_0.6.0          evaluate_1.0.0          BiocManager_1.30.25     ggprism_1.0.5          
# [136] nloptr_2.1.1            foreach_1.5.2           tweenr_2.0.3            httpuv_1.6.15           networkD3_0.4          
# [141] RANN_2.6.2              tidyr_1.3.1             openssl_2.2.2           purrr_1.0.2             polyclip_1.10-7        
# [146] future_1.34.0           clue_0.3-65             scattermore_1.2         ggforce_0.4.2           broom_1.0.7            
# [151] xtable_1.8-4            RSpectra_0.16-2         later_1.3.2             viridisLite_0.4.2       tibble_3.2.1           
# [156] beeswarm_0.4.0          IRanges_2.38.1          cluster_2.1.6           globals_0.16.3          BiocStyle_2.32.1

