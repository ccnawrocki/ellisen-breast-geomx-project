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
norm <- read.csv(file = "q3norm_nonlog.csv", row.names = 1) # The data needs to be linear (not log-transformed)


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
                        raw = immcts 
                        #cell_counts = immmeta$Number_of_cells_sampled
)

# Visualization
ha1 <- HeatmapAnnotation(type = immmeta$infiltration_type, 
                         col = list(type = c("ieb"="dodgerblue", "iib"="purple", "iic"="magenta")))
mat <- immdcvn$prop_of_all
Heatmap(mat, 
        cluster_columns = F,
        show_column_names = T, 
        width = ncol(immdcvn$beta)*unit(5, "mm"),
        height = nrow(immdcvn$beta)*unit(5, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5),
        name = "Proportion", 
        top_annotation = ha1, 
        column_split = immmeta$infiltration_type) |> draw(column_title = "Immune, Proportion")
Heatmap(mat,
        show_column_names = T, 
        width = ncol(immdcvn$beta)*unit(5, "mm"),
        height = nrow(immdcvn$beta)*unit(5, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5),
        name = "Proportion", 
        top_annotation = ha1, 
        column_split = immmeta$infiltration_type) |> draw(column_title = "Immune, Proportion")
Heatmap(mat,
        show_column_names = T, 
        width = ncol(immdcvn$beta)*unit(5, "mm"),
        height = nrow(immdcvn$beta)*unit(5, "mm"), 
        rect_gp = gpar(col = "white", lwd = 0.5),
        name = "Proportion", 
        top_annotation = ha1) |> draw(column_title = "Immune, Proportion")


# Barplot
d <- mat |> t() |> as.data.frame()
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

immdcvn_long$patient <- plyr::mapvalues(x = immdcvn_long$aoi_id, from = immmeta$aoi_id, to = immmeta$Patient_number) |> as.factor()
aoiorder <- immmeta |> arrange(Patient_number) |> pull(aoi_id)
immdcvn_long$aoi_id <- factor(immdcvn_long$aoi_id, levels = aoiorder)

p2 <- ggplot(immdcvn_long) + 
  geom_bar(mapping = aes(x = aoi_id, y = prop, fill = celltype, group = aoi_id), stat = "identity", color = NA) + 
  ggprism::scale_fill_prism() + 
  theme_dendro() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(title = "immroblast, Proportions")

p2 + geom_text(data=immmeta[aoiorder,],
               aes(label=infiltration_type, x=aoi_id, y=1.075, color=infiltration_type)) +
  scale_color_manual(values = c("ieb"="dodgerblue", "iib"="purple", "iic"="magenta")) + 
  geom_text(data=immmeta[aoiorder,],
            aes(label=Patient_number, x=aoi_id, y=1.175)) 

aoitree <- hclust(dist(immdcvn$prop_of_all |> t()))
dend <- aoitree |> as.dendrogram()

ddata_x <- dendro_data(dend)
ddata_x$labels$patient <- plyr::mapvalues(x = ddata_x$labels$label, from = immmeta$aoi_id, to = immmeta$Patient_number)
aoiorder <- aoitree$labels[aoitree$order]

p1 <- ggplot(immdcvn_long) + 
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
immdcvn_long$ncells <- plyr::mapvalues(x = immdcvn_long$aoi_id, from = immmeta$aoi_id, to = immmeta$Number_of_cells_sampled) |> as.character() |> as.numeric()
immdcvn_long$patient <- plyr::mapvalues(x = immdcvn_long$aoi_id, from = immmeta$aoi_id, to = immmeta$Patient_number)
immdcvn_long$group <- plyr::mapvalues(x = immdcvn_long$aoi_id, from = immmeta$aoi_id, to = immmeta$infiltration_type)

# Getting the discrete counts
immdcvn_long$count <- (immdcvn_long$ncells*immdcvn_long$prop) |> round(digits = 0)

# Should use random slopes model, since some patients exist across the groups, 
# according to NanoString. But, looking here, the study design is very unbalanced. 
# Also, there are only 8 patients, and patient 2 only has 1 aoi. Using a random
# slopes model will likely overfit the data or it will fail to converge.
table(immmeta$Patient_number, immmeta$infiltration_type)
#   ieb iib iic
# 1   6   0   0
# 2   0   1   0
# 3   1   2   4
# 4   0   5   0
# 5   0   3   0
# 6   5   0   0
# 7   8   0   0
# 8   2   3   0

# Another thing to note is that when looking at the bar plots shown above, it just 
# does not seem like the makeup of the aois within patient groups is super constant.
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
    
    testdata <- filter(immdcvn_long, celltype == ct)
    
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

immmeta$infiltration_type %<>% as.factor()
levels(immmeta$infiltration_type) <- levels(immdcvn_long$group)
all(levels(immmeta$infiltration_type) == levels(immdcvn_long$group))
# [1] TRUE

immmeta$Patient_number %<>% as.factor()
levels(immmeta$Patient_number) <- levels(immdcvn_long$patient)
all(levels(immmeta$Patient_number) == levels(immdcvn_long$patient))
# [1] TRUE

mm <- model.matrix(~1+infiltration_type, data = immmeta)
#mm <- model.matrix(~1+infiltration_type+Patient_number, data = immmeta)
contrast_list <- list(
  "ieb / iib" = mm[immmeta$infiltration_type == "ieb",] |> colMeans() - mm[immmeta$infiltration_type == "iib",] |> colMeans(),
  "ieb / iic" = mm[immmeta$infiltration_type == "ieb",] |> colMeans() - mm[immmeta$infiltration_type == "iic",] |> colMeans(),
  "iib / iic" = mm[immmeta$infiltration_type == "iib",] |> colMeans() - mm[immmeta$infiltration_type == "iic",] |> colMeans()
)

# Note that there are no warning messages when using the simplest poisson model.
res <- lapply(X = names(contrast_list), FUN = purrr::map, .x = unique(immdcvn_long$celltype), .f = proptestfunc)
res <- bind_rows(res)
res <- group_by(res, contrast) |> mutate(fdr = p.adjust(p = pval, method = "BH"))
res$padj <- p.adjust(p = res$pval, method = "bonferroni")
knitr::kable(res, format = "simple")

# Feeling better about these results. The estimates are relatively controlled.

# contrast    celltype         pval     Estimate        tstat          sigma         fdr        padj
# ----------  ---------  ----------  -----------  -----------  -------------  ----------  ----------
# ieb / iib   Tstr        0.9999971   -0.3637888   -0.0000036   9.988671e+04   0.9999971   1.0000000
# ieb / iib   Teff        0.9894071   19.0600221    0.0132767   1.435604e+03   0.9999971   1.0000000
# ieb / iib   Tisg        0.0055470   -0.2094899   -2.7734253   7.553470e-02   0.0092449   0.0832043
# ieb / iib   Tn          0.0000000   -0.7024109   -7.5728738   9.275350e-02   0.0000000   0.0000000
# ieb / iib   Tex         0.0000170    0.2160409    4.3016537   5.022280e-02   0.0000424   0.0002543
# ieb / iic   Tstr        0.9999961   -0.8083575   -0.0000049   1.639532e+05   0.9999961   1.0000000
# ieb / iic   Teff        0.9946583   18.6154534    0.0066949   2.780538e+03   0.9999961   1.0000000
# ieb / iic   Tisg        0.0000000   -1.0542932   -9.9031484   1.064604e-01   0.0000000   0.0000000
# ieb / iic   Tn          0.9861349   17.7761974    0.0173782   1.022903e+03   0.9999961   1.0000000
# ieb / iic   Tex         0.0000359    0.5121141    4.1322390   1.239314e-01   0.0000898   0.0005389
# iib / iic   Tstr        0.9999979   -0.4445687   -0.0000026   1.708525e+05   0.9999979   1.0000000
# iib / iic   Teff        0.9998866   -0.4445687   -0.0001421   3.129273e+03   0.9999979   1.0000000
# iib / iic   Tisg        0.0000000   -0.8448033   -7.3914512   1.142946e-01   0.0000000   0.0000000
# iib / iic   Tn          0.9855871   18.4786083    0.0180649   1.022903e+03   0.9999979   1.0000000
# iib / iic   Tex         0.0214648    0.2960732    2.2997074   1.287439e-01   0.0536620   0.3219721

# I think using difference in average proportion will work for visualizations.
mm <- model.matrix(~0+infiltration_type, data = immmeta)
contrast_list <- list(
  "ieb / iib" = mm[immmeta$infiltration_type == "ieb",] |> colMeans() - mm[immmeta$infiltration_type == "iib",] |> colMeans(),
  "ieb / iic" = mm[immmeta$infiltration_type == "ieb",] |> colMeans() - mm[immmeta$infiltration_type == "iic",] |> colMeans(),
  "iib / iic" = mm[immmeta$infiltration_type == "iib",] |> colMeans() - mm[immmeta$infiltration_type == "iic",] |> colMeans()
)
meanprops <- immdcvn_long |> group_by(celltype, group) |> summarise(mean_prop = mean(count/ncells))
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
# Bonferroni correction is probably the way to go here. 
ggplot(data = res) + 
  geom_bar(mapping = aes(y = celltype, x = mean_diff, fill = padj < 0.05), stat = "identity") + 
  facet_grid(.~contrast) +
  scale_fill_manual(values = c("blue", "red")) + 
  ggthemes::theme_few() + 
  geom_vline(xintercept = 0, linewidth = 0.2)

# Using limma... still getting warning messages
# Have to fit using the abundance score, not the proportion
mat <- (immdcvn$beta+1) |> log2()
mm <- model.matrix(~0+infiltration_type, data = immmeta)
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

write.csv(res, "immune_deconvolution_differential_abundance_testing_results.csv")
write.csv(deres, "immune_deconvolution_differential_abundance_testing_results_limma.csv")

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

