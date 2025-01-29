## GeoDiff Pipeline and Downstream DE ##
### Cole Nawrocki ###

# I am trying an alternative pipeline so that we can use the lme4 mixed-effects model for differential expression analysis. 
# Link to the vignette I followed: https://www.bioconductor.org/packages/release/bioc/vignettes/GeoDiff/inst/doc/Workflow_WTA_kidney.html

# Main points: 
## I am performing a background-based normalization, called "Poisson threshold normalization" 
## This method includes new QC as well
## I tried modeling on the normalized data, using the following model formula: y ~ 1 + group + (1 + 1|coreid), but it did not
## improve power.
## "coreid" is a factor variable that defines the piece of tissue (core biopsy) from which the aoi comes from. This seemed more 
## appropriate than patient for the random intercept, since the questions we are attempting to answer are about microenvironment, 
## not about patients in general. aois that are from the same microenvironment should be "more correlated" than aois that are from
## the same patient. However, using a random slopes model actually seems best and is supported by the NanoString GeoMx-Tools
## vignette. Formula: y ~ 1 + group + (1 + group|Patient_number)
## I normalized for each sequencing run separately as to mitigate batch effects. 

# Clearing environment and making sure we are in the correct micromamba environment.
rm(list = ls())
.libPaths()

# Packages
library(lme4)
library(lmerTest)
library(parallel)
library(multcomp)
library(dplyr)
library(stringr)
library(magrittr)
library(GeoDiff)
library(GeomxTools)
library(ggplot2)
options(mc.cores = 10)
set.seed(2001)

## Alternative Pipeline --------------------------------------------------------
d <- readRDS("geomx_obj_unprocessed.RDS") # This should be the raw probe-level data from before QC or transformation
featureType(d) 
# [1] "Probe"

paste("## of Negative Probes:", sum(fData(d)$Negative))
# [1] "## of Negative Probes: 139"

d@phenoData@data$run <- substr(x = d@phenoData@data$`slide name`, start = 10, stop = 11) # Adding the run variable

d <- fitPoisBG(d)
d <- fitPoisBG(d, groupvar = "run")
d_diag <- diagPoisBG(d, split = T)
notes(d_diag)$disper_sp # The vignette says that a dispersion >2 is worrisome, so we are okay. If it was above 2, we might remove outliers and re-run the background model.
# [1] 1.296283
which(assayDataElement(d_diag, "low_outlier") == 1, arr.ind = T)
which(assayDataElement(d_diag, "up_outlier") == 1, arr.ind = T)

all0probeidx <- which(rowSums(exprs(d))==0) # Remove any probes with 0 counts
if (length(all0probeidx) > 0) {
  d <- d[-all0probeidx, ]
}
d <- aggreprobe(d, use = "cor") # Aggregate the probes
d <- BGScoreTest(d, removeoutlier = T, split = F) # Identify outlier targets

sum(fData(d)[["pvalues"]] < 1e-3, na.rm = TRUE) # ~14000 targets retained
# [1] 13676

d <- fitNBth(d, split = TRUE) # Estimating the size factors
features_high <- rownames(fData(d))[fData(d)$feature_high_fitNBth == 1] # The 1500 features that are well above background

bgMean <- mean(fData(d)$featfact, na.rm = TRUE) # The threshold is a bit above the mean background level
notes(d)[["threshold"]]
# [1] 497.597
bgMean
# [1] 466.7266

cor(d$sizefact, d$sizefact_fitNBth) # These are supposed to be somewhat correlated. The value of r is low here, but that might just be because of some strange aois.
plot(d$sizefact, d$sizefact_fitNBth, xlab = "Background Size Factor",
     ylab = "Signal Size Factor")
abline(a = 0, b = 1)

# These visualizations just show that the size factors are correlated with quantiles for the data
posdat <- d[-which(fData(d)$CodeClass == "Negative"), ]
posdat <- exprs(posdat)

quan <- sapply(c(0.75, 0.8, 0.9, 0.95), function(y)
  apply(posdat, 2, function(x) quantile(x, probs = y)))

corrs <- apply(quan, 2, function(x) cor(x, d$sizefact_fitNBth))
names(corrs) <- c(0.75, 0.8, 0.9, 0.95)

corrs

quan75 <- apply(posdat, 2, function(x) quantile(x, probs = 0.75))

d <- QuanRange(d, split = FALSE, probs = c(0.75, 0.8, 0.9, 0.95))

corrs <- apply(pData(d)[, as.character(c(0.75, 0.8, 0.9, 0.95))], 2, function(x)
  cor(x, d$sizefact_fitNBth))

names(corrs) <- c(0.75, 0.8, 0.9, 0.95)

corrs

# Filtering out ROIs with low signal compared to background. I used a less stringent cutoff than the vignette (0.5 instead of 2).
ROIs_high <- sampleNames(d)[which((quantile(fData(d)[["para"]][, 1],
                                            probs = 0.9, na.rm = TRUE) -
                                     notes(d)[["threshold"]])*d$sizefact_fitNBth>0.5)]

# Actually normalizing the data
d <- fitPoisthNorm(object = d,
                   split = TRUE,
                   ROIs_high = ROIs_high,
                   threshold_mean = bgMean,
                   sizescalebythreshold = TRUE)
saveRDS(d, "geomx_obj_geodiff_processed.RDS")
dataobj <- d

norm <- assayDataElement(d[fData(d)[["pvalues"]] < 1e-3, ROIs_high], "normmat_sp") |> na.omit()
write.csv(norm, "ptnorm.csv")

annot <- pData(d)[ROIs_high,]

dat_plot <- cbind(annot[, "segment", drop = F], t(norm))
dat_plot <- cbind(dat_plot, ROI_ID = ROIs_high)

dat_plot <- reshape2::melt(dat_plot, id.vars = c("ROI_ID", "segment"))

ggplot(dat_plot, aes(x = value)) +
  geom_density(aes(group = ROI_ID), fill = NA) +
  facet_wrap(~segment) +
  ggtitle("Poisson threshold normalization") +
  labs(x = "Poisson Threshold Normalized Value (log2)") + 
  ggthemes::theme_few()

ggplot() +
  ggridges::geom_density_ridges_gradient(data = dat_plot, mapping = aes(x = value, y = ROI_ID, fill = segment)) +  
  facet_wrap(.~segment, scales = "free") + 
  ggtitle("Poisson threshold normalization") +
  labs(x = "Poisson Threshold Normalized Value (log2)") + 
  ggthemes::theme_few() + 
  theme(axis.text.y = element_blank())
  
d <- normalize(d,
               norm_method = "quant", 
               desiredQuantile = .75,
               toElt = "q_norm")
qnorm <- assayDataElement(d[fData(d)[["pvalues"]] < 1e-3, ROIs_high], "q_norm") |> na.omit()
qnorm <- log2(1+qnorm)

dat_plot <- cbind(annot[, "segment", drop = F], t(qnorm))
dat_plot <- cbind(dat_plot, ROI_ID = ROIs_high)

dat_plot <- reshape2::melt(dat_plot, id.vars = c("ROI_ID", "segment"))

ggplot(dat_plot, aes(x = value)) +
  geom_density(aes(group = ROI_ID), fill = NA) +
  facet_wrap(~segment) +
  ggtitle("Q3 normalization") +
  labs(x = "Q3 Normalized Value (log2 + 1)") + 
  ggthemes::theme_few()

ggplot() +
  ggridges::geom_density_ridges_gradient(data = dat_plot, mapping = aes(x = value, y = ROI_ID, fill = segment)) +  
  facet_wrap(.~segment, scales = "free") + 
  ggtitle("Q3 normalization") +
  labs(x = "Q3 Normalized Value (log2 + 1)") + 
  ggthemes::theme_few() + 
  theme(axis.text.y = element_blank())

vnorm <- (assayDataElement(d[fData(d)[["pvalues"]] < 1e-3, ROIs_high], "exprs") |> na.omit() |> edgeR::DGEList() |> edgeR::calcNormFactors() |> limma::voom())$E

dat_plot <- cbind(annot[, "segment", drop = F], t(vnorm))
dat_plot <- cbind(dat_plot, ROI_ID = ROIs_high)

dat_plot <- reshape2::melt(dat_plot, id.vars = c("ROI_ID", "segment"))

ggplot(dat_plot, aes(x = value)) +
  geom_density(aes(group = ROI_ID), fill = NA) +
  facet_wrap(~segment) +
  ggtitle("Voom normalization") +
  labs(x = "Voom Normalized Value (log2)") + 
  ggthemes::theme_few()

ggplot() +
  ggridges::geom_density_ridges_gradient(data = dat_plot, mapping = aes(x = value, y = ROI_ID, fill = segment)) +  
  facet_wrap(.~segment, scales = "free") + 
  ggtitle("Voom normalization") +
  labs(x = "Voom Normalized Value (log2)") + 
  ggthemes::theme_few() + 
  theme(axis.text.y = element_blank())

## Expanding the metadata and cleaning it --------------------------------------
annot$aoi_id <- rownames(annot) |> gsub(pattern = "-", replacement = ".")
extra_meta <- readxl::read_excel("BreastCaGeoMx_v7CN.xlsx", sheet = 2)
extra_meta$aoi_id <- gsub(pattern = "-", replacement = ".", x = extra_meta$`Sample Sheet Code`) |> paste(".dcc", sep = "")

meta_combined <- inner_join(x = annot, y = extra_meta, by = "aoi_id")
rownames(meta_combined) <- meta_combined$aoi_id
meta_combined <- mutate(meta_combined, aoi_codes_advanced = as.factor(`AOI coded_Advanced`))

label_map <- c(
  "1"="ieb",
  "2"="iec",
  "3"="iib",
  "4"="iic",
  "5"="idb",
  "6"="idc"
)

meta_combined$infiltration_type <- str_sub(string = meta_combined$aoi_codes_advanced, start = -1)
meta_combined$infiltration_type <- plyr::mapvalues(x = meta_combined$infiltration_type, from = names(label_map), to = label_map)

meta_cleaned <- dplyr::select(.data = meta_combined, c(15, 22, 24, 25, 27, 28, 30:38, 40, 41))
colnames(meta_cleaned) %<>% gsub("\\(in sq. micron\\)", "um2", x = .)
colnames(meta_cleaned) %<>% gsub(" ", "_", x = .)
colnames(meta_cleaned) %<>% gsub("\\+", "_positive_", x = .)

meta_cleaned$infiltration_category <- substr(x = meta_cleaned$infiltration_type, start = 2, stop = 2)
meta_cleaned$infiltration_category[meta_cleaned$infiltration_category == "e"] <- "excluded"
meta_cleaned$infiltration_category[meta_cleaned$infiltration_category == "i"] <- "inflamed"
meta_cleaned$infiltration_category[meta_cleaned$infiltration_category == "d"] <- "desert"
meta_cleaned$group <- ifelse(meta_cleaned$infiltration_category == "desert", 
                             yes = meta_cleaned$infiltration_type |> str_sub(start = 1, end = 2), 
                             no = meta_cleaned$infiltration_type)
meta_cleaned$group <- paste(meta_cleaned$AOI_code, meta_cleaned$group, sep = "_")

openxlsx::write.xlsx(x = meta_cleaned, file = "meta_cleaned_v7.xlsx")
write.csv(x = meta_cleaned, file = "meta_cleaned_v7.csv")

## lme4, LRT -------------------------------------------------------------------
defun <- function(data_df) {
  h0 <- lmer(y ~ 1 + (1+group|Patient_number), data_df, REML = FALSE)
  h1 <- lmer(y ~ 1 + group + (1+group|Patient_number), data_df, REML = FALSE)
  beta <- fixef(h1)[["groupfibroblast_iib"]]
  pval <- anova(h0, h1)["h1", 'Pr(>Chisq)']
  out <- data.frame(beta = beta, pval = pval)
  return(out)
}

# Set up the data... just modeling using the aois in the comparison. It would be complicated and super memory intense otherwise.
colnames(norm) <- gsub(pattern = "-", replacement = ".", x = colnames(norm))
norm <- t(norm) |> as.data.frame()
norm$aoi_id <- rownames(norm)
alldata <- inner_join(x = meta_cleaned, y = norm, by = "aoi_id")
alldata$Patient_number %<>% as.factor()
alldata <- tidyr::pivot_longer(data = alldata, cols = 24:ncol(alldata), names_to = "FeatureName", values_to = "y")

d <- alldata |> filter(group %in% c("fibroblast_ieb", "fibroblast_iib"))
d <- d |> split(f = d$FeatureName)
res <- parallel::mclapply(X = d, FUN = defun, mc.cores = 10)
res <- bind_rows(res, .id = "target")
res$fdr <- p.adjust(res$pval, 'BH')

ggplot() + geom_point(data = res, mapping = aes(x = beta, y = -log10(pval)))
ggplot() + geom_point(data = res, mapping = aes(x = beta, y = -log10(fdr)))

## lme4, Wald test -------------------------------------------------------------
# Modeling using all of the data.
d <- alldata |> split(f = alldata$FeatureName)
mm <- model.matrix(~1+group, data = d$A2M)
contr <- (mm[d$A2M$group == "fibroblast_ieb", ] |> colMeans()) - (mm[d$A2M$group == "fibroblast_iib", ] |> colMeans())
lmem_de <- function(data_df) {
  modelout <- lmer(y ~ 1 + group + (1+group|Patient_number), data_df, REML = FALSE)
  out <- contest(model = modelout, L = contr, joint = F, rhs = 0)
  return(out)
}

res2 <- parallel::mclapply(X = d, FUN = lmem_de, mc.cores = 10)
res2 <- bind_rows(res2, .id = "target")
res2$fdr <- p.adjust(res2$`Pr(>|t|)`, 'BH')
ggplot() + geom_point(data = res2, mapping = aes(x = Estimate, y = -log10(`Pr(>|t|)`)))
ggplot() + geom_point(data = res2, mapping = aes(x = Estimate, y = -log10(fdr)))

# These results seem promising. 
write.csv(x = res2, "geodiff_lme4_fibieb_vs_fibiib.csv")

## Comparing the two tests -----------------------------------------------------
intersect((res |> filter(pval < 0.05) |> pull(target)), (res2 |> filter(`Pr(>|t|)` < 0.05) |> pull(target))) |> length()
intersect((res |> filter(fdr < 0.05) |> pull(target)), (res2 |> filter(fdr < 0.05) |> pull(target))) |> length()

## NB Wald test ----------------------------------------------------------------
# This took ~2 hours to run just for the aois in one comparison.
# I tried with a couple different offsets. The results were all similar.
cts <- assayDataElement(dataobj[fData(dataobj)[["pvalues"]] < 1e-2, ROIs_high], "exprs") |> na.omit()
colnames(cts) <- gsub(pattern = "-", replacement = ".", x = colnames(cts))
cts_obj <- edgeR::DGEList(counts = cts) |> edgeR::calcNormFactors()
effectivelibsize <- cts_obj$samples$lib.size*cts_obj$samples$norm.factors
meta_cleaned$effectivelibsize <- effectivelibsize
meta_cleaned$ln_effectivelibsize <- log(effectivelibsize, base = exp(1))
meta_cleaned$effsize <- meta_cleaned$sizefact_fitNBth*colSums(cts)
meta_cleaned$ln_effsize <- log(meta_cleaned$effsize, base = exp(1))
cts <- t(cts) |> as.data.frame()
cts$aoi_id <- rownames(cts)
alldata <- inner_join(x = meta_cleaned, y = cts, by = "aoi_id")
alldata$Patient_number %<>% as.factor()
alldata <- tidyr::pivot_longer(data = alldata, cols = 24:ncol(alldata), names_to = "FeatureName", values_to = "y")

mm <- model.matrix(~group, data = d$A2M)
contr <- (mm[d$A2M$group == "tumor_iic", ] |> colMeans()) - (mm[d$A2M$group == "tumor_id", ] |> colMeans())
nb_de <- function(data_df) {
  tryCatch({
    modelout <- lme4::glmer.nb(y ~ offset(ln_effectivelibsize) + group + (group|Patient_number), data_df, verbose = F)
    out <- summary(modelout)$coefficients["grouptumor_iic", , drop = F]
    #testout <- multcomp::glht(modelout, linfct = matrix(contr, nrow = 1)) |> summary(test = Chisqtest())
    #coef <- testout$test$coefficients
    #p <- testout$test$pvalue
    #out <- data.frame("coef"=coef, "p"=p)
    return(out)
  }, error = function(e){;})
}
d <- alldata |> filter(group %in% c("tumor_iic", "tumor_id"))
d <- d |> split(f = d$FeatureName)
res3 <- parallel::mclapply(X = d, FUN = nb_de, mc.cores = 10)
res3 <- purrr::map(res3, as.data.frame) |> bind_rows(.id = "target")
res3$fdr <- p.adjust(res3$`Pr(>|z|)`, 'BH')
ggplot() + geom_point(data = res3, mapping = aes(x = Estimate, y = -log10(`Pr(>|z|)`)))
ggplot() + geom_point(data = res3, mapping = aes(x = Estimate, y = -log10(fdr)))

nb_de <- function(data_df) {
  tryCatch({
    modelout <- lme4::glmer.nb(y ~ offset(ln_effsize) + group + (group|Patient_number), data_df, verbose = F)
    out <- summary(modelout)$coefficients["grouptumor_iic", , drop = F]
    #testout <- multcomp::glht(modelout, linfct = matrix(contr, nrow = 1)) |> summary(test = Chisqtest())
    #coef <- testout$test$coefficients
    #p <- testout$test$pvalue
    #out <- data.frame("coef"=coef, "p"=p)
    return(out)
  }, error = function(e){;})
}
res4 <- parallel::mclapply(X = d, FUN = nb_de, mc.cores = 10)
res4 <- purrr::map(res4, as.data.frame) |> bind_rows(.id = "target")
res4$fdr <- p.adjust(res4$`Pr(>|z|)`, 'BH')
ggplot() + geom_point(data = res4, mapping = aes(x = Estimate, y = -log10(`Pr(>|z|)`)))
ggplot() + geom_point(data = res4, mapping = aes(x = Estimate, y = -log10(fdr)))

## limma with the new normalized data ------------------------------------------
norm <- assayDataElement(dataobj[fData(dataobj)[["pvalues"]] < 1e-2, ROIs_high], "normmat_sp") |> na.omit()
mm <- model.matrix(~group, data = meta_cleaned)
colnames(mm) <- gsub(pattern = "group", replacement = "", colnames(mm))
correl = limma::duplicateCorrelation(object = norm, design = mm, block = meta_cleaned$Patient_number)
lmfit <- limma::lmFit(object = norm, design = mm, block = meta_cleaned$Patient_number, correlation = correl$consensus.correlation)
contr <- (mm[meta_cleaned$group == "fibroblast_ieb", ] |> colMeans()) - (mm[meta_cleaned$group == "fibroblast_iib", ] |> colMeans())
res5 <- limma::contrasts.fit(fit = lmfit, contrasts = contr)
res5 <- limma::eBayes(res5)
res5 <- limma::topTable(res5, number = Inf)
ggplot() + geom_point(data = res5, mapping = aes(x = logFC, y = -log10(adj.P.Val)))

## Notes and lessons learned ---------------------------------------------------
# Poisson threshold normalization works well. Might want to lower threshold to 0.5, then start using this method.
# Pretty sure that you can model on the normalized data. It is just a stronger transformation than Q3. 
# Using with limma, it looks better (more significant values). 
# Using with lme4, I was able to model with random slopes and intercepts (1+group|Patient_number), which is more "correct" than what
# limma does. limma models with random intercepts only (1+1|Patient_number).
# Using coreid for the random intercept did not ameliorate anything. 
# I think that using limma with this new normalization method is a good option. 
# I think that using lme4 with this new method is also a good option. For GO, using values that are significant (un-adjusted p) is 
# fine. I think that care needs to be taken when analyzing singular genes that are not significant for adjusted p.
# The normalization method no longer outputs "NegProbe-WTX" which is needed for spatial deconvolution (I think), so I have to 
# figure that out still.
