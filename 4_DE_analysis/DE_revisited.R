## DE Modeling Revisited ##
### Cole Nawrocki ###

# Notes: 
# It seems that using q3 eff lib sizes will yield way less results.
# Do I need to be normalizing each compartment separately?
# What if I use the voom precision weights as weights in the lme4 model? 
# The design makes the results better in limma? How is this happening?

# Clearing environment and making sure we are in the correct micromamba environment.
rm(list = ls())
.libPaths()

# [1] "/Users/ccn22/micromamba/envs/geomx-env/lib/R/library"

# Packages
library(lme4)
library(lmerTest)
library(parallel)
library(multcomp)
library(dplyr)
library(stringr)
library(magrittr)
library(ggplot2)
library(ComplexHeatmap)
options(mc.cores = 10)
set.seed(2001)

# Good trick
theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

# Data
meta <- openxlsx::read.xlsx("meta_cleaned_v5.xlsx")
rownames(meta) <- meta$aoi_id
cts <- read.csv(file = "counts.csv", row.names = 1)
counts <- edgeR::DGEList(counts = cts)
counts <- edgeR::calcNormFactors(object = counts, method = "TMM")
all(colnames(counts$counts) == rownames(meta))
# [1] TRUE

meta$efflibsize <- counts$samples$lib.size*counts$samples$norm.factors
qs <- apply(X = cts, MARGIN = 2, FUN = quantile, 0.75)
nfs <- qs / EnvStats::geoMean(qs)
q3norm <- sweep(x = cts, MARGIN = 2, STATS = nfs, FUN = "/")
meta$q3efflibsize <- colSums(q3norm)
meta$group <- ifelse(test = meta$infiltration_category == "desert", yes = "cold", no = "hot")

cts <- t(cts) |> as.data.frame()
cts$aoi_id <- rownames(cts)
alldata <- inner_join(x = meta, y = cts, by = "aoi_id")
alldata$Patient_number %<>% as.factor()
alldata$group %<>% as.factor()
alldata <- tidyr::pivot_longer(data = alldata, cols = 19:ncol(alldata), names_to = "FeatureName", values_to = "y")

norm <- read.csv(file = "q3norm.csv", row.names = 1)
norm <- t(norm) |> as.data.frame()
norm$aoi_id <- rownames(norm)
alldata <- inner_join(x = meta, y = norm, by = "aoi_id")
alldata$Patient_number %<>% as.factor()
alldata$group %<>% as.factor()
alldata <- tidyr::pivot_longer(data = alldata, cols = 19:ncol(alldata), names_to = "FeatureName", values_to = "y")

## FIGURING OUT THE MODEL ------------------------------------------------------
# Tumor
dt <- alldata |> filter(AOI_code == "tumor")
dt <- dt |> split(f = dt$FeatureName)

# Using TROP2 as an example
dd <- dt$TACSTD2

# Example contrast
mm <- model.matrix(~ 1+group, data = dd)
contr <- (mm[dd$group == "cold", ] |> colMeans()) - (mm[dd$group == "hot", ] |> colMeans())

# Modeling
h0 <- glmer.nb(y ~ 1+(1+1|Patient_number)+offset(log(efflibsize)), dd)
h1 <- glmer.nb(y ~ 1+group+(1+1|Patient_number)+offset(log(efflibsize)), dd)
h0 <- lme4::lmer(y ~ 1+(1+1|Patient_number), dd, REML = F)
h1 <- lme4::lmer(y ~ 1+group+(1+1|Patient_number), dd, REML = F)

# LRT
beta <- fixef(h1)%*%contr
pval <- anova(h0, h1)['h1', 'Pr(>Chisq)']
data.frame(beta=beta, pval=pval)

# Wald
h1 |> glht(linfct = matrix(contr, nrow = 1), alternative = "two.sided", rhs = 0) |> summary(test = Ftest())
h1 |> glht(linfct = matrix(contr, nrow = 1), alternative = "two.sided", rhs = 0) |> summary(test = Chisqtest())

# Wald, emmeans
library(emmeans)
ref_grid(h1)
EMM <- emmeans(object = h1, specs = ~0+group, lmer.df = "satterthwaite")
EMM
testout <- emmeans::contrast(object = EMM, list("test"=c(1,-1))) |> summary()
testout
emmeans::emmip(EMM, formula = ~1+group)
plot(EMM)

# Fibroblast
df <- alldata |> filter(AOI_code == "fibroblast")
df <- df |> split(f = df$FeatureName)

# Using TROP2 as an example
dd <- df$TACSTD2

# Example contrast
mm <- model.matrix(~ 1+group, data = dd)
contr <- (mm[dd$group == "cold", ] |> colMeans()) - (mm[dd$group == "hot", ] |> colMeans())

# Modeling
h1 <- glmer.nb(y ~ 1+group+(1+1|Patient_number)+offset(log(efflibsize)), dd)
h0 <- glmer.nb(y ~ 1+(1+1|Patient_number)+offset(log(efflibsize)), dd)

# LRT
beta <- fixef(h1)%*%contr
pval <- anova(h0, h1)['h1', 'Pr(>Chisq)']
data.frame(beta=beta, pval=pval)

# Wald
h1 |> glht(linfct = matrix(contr, nrow = 1), alternative = "two.sided", rhs = 0) |> summary(test = Ftest())
h1 |> glht(linfct = matrix(contr, nrow = 1), alternative = "two.sided", rhs = 0) |> summary(test = Chisqtest())

# Wald, emmeans
library(emmeans)
ref_grid(h1)
EMM <- emmeans(h1, "group")
EMM
testout <- emmeans::contrast(object = EMM, list("test"=contr)) |> summary()
testout
emmeans::emmip(EMM, formula = ~group)
plot(EMM)

## MODELING --------------------------------------------------------------------
# Function for doing the modeling
de_model <- function(data_df) {
  h1 <- glmer.nb(y ~ 1+group+(1+1|Patient_number)+offset(log(q3efflibsize)), data_df)
  h0 <- glmer.nb(y ~ 1+(1+1|Patient_number)+offset(log(q3efflibsize)), data_df)
  return(list("h0" = h0, "h1" = h1))
}

de_model <- function(data_df) {
  h1 <- lme4::lmer(y ~ 1+group+(1+1|Patient_number), data_df, REML = F)
  h0 <- lme4::lmer(y ~ 1+(1+1|Patient_number), data_df, REML = F)
  return(list("h0" = h0, "h1" = h1))
}

library(furrr)
future::plan(multisession, workers = 10)
tumormodels <- furrr::future_map(.x = dt, .f = de_model, .progress = T)
saveRDS(tumormodels, file = "tumordemodels_q3_lmer.RDS")
tumormodels <- readRDS("tumordemodels_tmm.RDS")
de_test <- function(mods) {
  beta <- fixef(mods$h1)%*%contr
  pval <- anova(mods$h0, mods$h1)["mods$h1", "Pr(>Chisq)"]
  return(data.frame(beta=beta, pval=pval))
}

contrast_list <- list("cold / hot" = contr)
res <- mclapply(tumormodels, de_test, mc.cores = 10)
res <- bind_rows(res, .id = "target")
res$fdr <- p.adjust(p = res$pval, method = "BH")

ggplot(res) + geom_point(mapping = aes(x = beta, y = -log10(fdr)))

# Start: 
fibroblastmodels <- mclapply(X = dt, FUN = de_model, mc.cores = 10)

# Function for doing the testing
de_test <- function(outmod, .contr, .test_type = "univariate") {
  if (.test_type == "univariate") {
    testout <- outmod |> 
      glht(linfct = matrix(contrast_list[[.contr]], nrow = 1), alternative = "two.sided", rhs = 0) |> 
      summary(test = univariate()) %$%test
    out <- data.frame(
      contrast = .contr,
      coefficient = testout$coefficients[[1]],
      sigma = testout$sigma[[1]],
      tstat = testout$tstat[[1]],
      pvalue = testout$pvalues[[1]],
      test_type = testout$type
    )
  }
  
  if (.test_type == "Ftest") {
    testout <- outmod |> 
      glht(linfct = matrix(contrast_list[[.contr]], nrow = 1), alternative = "two.sided", rhs = 0) |> 
      summary(test = Ftest()) %$%test
    out <- data.frame(
      contrast = .contr,
      coefficient = testout$coefficients[[1]],
      sigma = testout$sigma[[1]],
      df1 = testout$df[1],
      df2 = testout$df[2],
      fstat = testout$fstat,
      tstat = testout$tstat[[1]],
      pvalue = testout$pvalue,
      test_type = testout$type
    )
  }
  
  if (.test_type == "Chisqtest") {
    testout <- outmod |> 
      glht(linfct = matrix(contrast_list[[.contr]], nrow = 1), alternative = "two.sided", rhs = 0) |> 
      summary(test = Chisqtest()) %$%test
    out <- data.frame(
      contrast = .contr,
      coefficient = testout$coefficients[[1]],
      sigma = testout$sigma[[1]],
      df = 1,
      fstat = testout$fstat,
      tstat = testout$tstat[[1]],
      pvalue = testout$pvalue,
      test_type = testout$type
    )
  }
  
  if (.test_type == "Satterthwaite") {
    out <- lmerTest::contest(model = outmod, L = contrast_list[[.contr]], rhs = 0, joint = F, confint = F, ddf = "Satterthwaite")
    out <- cbind("contrast" = .contr, "type" = "Satterthwaite", out)
  }
  
  if (.test_type == "Kenward-Roger") {
    out <- lmerTest::contest(model = outmod, L = contrast_list[[.contr]], rhs = 0, joint = F, confint = F, ddf = "Kenward-Roger")
    out <- cbind("contrast" = .contr, "type" = "Kenward-Roger", out)
  }
  
  return(out)
}

## MODELING --------------------------------------------------------------------
# I will do the modeling for the entire dataset on eristwo, then I will load it back in
# I have a conda environment on eristwo called lme4-env that I activate in an interactive session

# Saving the data (it is already all set up)
saveRDS(d, file = "tumor_split_data_for_modeling.RDS")

# Now, switching over to my terminal to use the HPC.

# $ scp /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/analysis/ellisen-breast/split_data_for_modeling.RDS ccn22@eristwo.partners.org:/PHShome/ccn22/scratch/split_data_for_modeling.RDS 
# $ ssh -Y ccn22@eristwo.partners.org
# $ srun --pty -p interactive -c 35 --mem=50000MB -t 300 /bin/bash
# $ cd scratch
# $ conda activate lme4-env
# $ R

# library(lme4)
# library(lmerTest)
# library(parallel)
# library(magrittr)
# options(mc.cores = 35)
# set.seed(2001)
# 
# d <- readRDS(file = "tumor_split_data_for_modeling.RDS")
# 
# de_model <- function(data_df) {
#   h1 <- glmer.nb(y ~ 1+group+(1+1|Patient_number)+offset(log(efflibsize)), dd)
#   h0 <- glmer.nb(y ~ 1+(1+1|Patient_number)+offset(log(efflibsize)), dd)
#   print(unique(data_df$FeatureName))
#   return(list("h0" = h0, "h1" = h1))
# }
# 
# modeling_res <- parallel::mclapply(X = d, FUN = de_model, mc.cores = 30)
# saveRDS(modeling_res, file = "tumor_de_model_results.RDS")
# 
# q()

# $ conda deactivate
# $ exit
# $ exit
# $ scp ccn22@eristwo.partners.org:/PHShome/ccn22/scratch/tumor_de_model_results.RDS /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/analysis/ellisen-breast/tumor_de_model_results.RDS

# It worked. It is ~ 900 MB, but this is reasonable for 14k mixed models.
full_res <- readRDS(file = "lmem_de_model_results_ml.RDS")
object.size(full_res)
# 870199896 bytes
870199896/1e9
# [1] 0.8701999

# I know that Bogang wants some other comparisons done that exclude the one tumor_idc aoi. So, I will remove this aoi and fit another model. 
d0 <- purrr::map(.x = d, .f = filter, (infiltration_type != "idc"))

# Saving the data (it is already all set up)
saveRDS(d0, file = "split_data_for_modeling_minus_idc.RDS")

# Now, switching over to my terminal to use the HPC.

# $ scp /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/analysis/ellisen-breast/split_data_for_modeling_minus_idc.RDS ccn22@eristwo.partners.org:/PHShome/ccn22/scratch/split_data_for_modeling_minus_idc.RDS 
# $ ssh -Y ccn22@eristwo.partners.org
# $ srun --pty -p interactive -c 35 --mem=50000MB -t 300 /bin/bash
# $ cd scratch
# $ conda activate lme4-env
# $ R

# library(lme4)
# library(lmerTest)
# library(parallel)
# library(magrittr)
# options(mc.cores = 35)
# set.seed(2001)
# 
# d <- readRDS(file = "split_data_for_modeling_minus_idc.RDS")
# 
# lmem_de_model <- function(data_df) {
#   lme4mod <- lme4::lmer(y ~ 1+group+(1+group|Patient_number), data_df, REML = T)
#   return(lme4mod)
# }
# 
# modeling_res <- parallel::mclapply(X = d, FUN = lmem_de_model, mc.cores = 35)
# saveRDS(modeling_res, file = "lmem_de_model_results_reml_minus_idc.RDS")
# 
# q()

# $ conda deactivate
# $ exit
# $ exit
# $ scp ccn22@eristwo.partners.org:/PHShome/ccn22/scratch/lmem_de_model_results_reml_minus_idc.RDS /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/analysis/ellisen-breast/lmem_de_model_results_reml_minus_idc.RDS 

# NOTE: I ran this twice: once using REML and once using ML.

# It also worked. It is of similar size.
modeling_res_minus_idc <- readRDS(file = "lmem_de_model_results_ml_minus_idc.RDS")
object.size(modeling_res_minus_idc)
# 870085416 bytes
870085416/1e9
# [1] 0.8700854

## TESTING EXAMPLE -------------------------------------------------------------
# Doing an example: fibroblast_ieb / fibroblast_id
contrast_list <- list(
  "tumor_iic / tumor_iec" = (mm[meta$group == "tumor_iic", ] |> colMeans()) - (mm[meta$group == "tumor_iec", ] |> colMeans()),
  "fibroblast_ieb / fibroblast_id" = (mm[meta$group == "fibroblast_ieb", ] |> colMeans()) - (mm[meta$group == "fibroblast_id", ] |> colMeans())
)

testing_res <- mclapply(X = full_res, FUN = lmem_de_test, .contr = "fibroblast_ieb / fibroblast_id", .test_type = "Chisqtest")
testing_res <- bind_rows(testing_res, .id = "target")
testing_res$fdr <- p.adjust(p = testing_res$pvalue, method = "BH")

# We have some DEGs
(testing_res$fdr < 0.05) |> sum()
# [1] 276

# Corresponding volcano plot... looks good
testing_res$significance <- "NS"
testing_res$significance[testing_res$pvalue < 0.05] <- "p < 0.05"
testing_res$significance[testing_res$fdr < 0.05] <- "fdr < 0.05"
testing_res$significance <- factor(testing_res$significance, levels = c("NS", "p < 0.05", "fdr < 0.05"))
ggplot(data = testing_res) + 
  geom_hline(yintercept = -log10(0.05), linewidth = 0.2) +
  geom_vline(xintercept = c(-1, 0, 1), linewidth = 0.2) +
  geom_point(mapping = aes(x = coefficient, y = -log10(pvalue), color = significance), alpha = 0.75, stroke = 0.1) + 
  geom_text(data = data.frame(x = c(min(testing_res$coefficient)+0.1, max(testing_res$coefficient)-0.1), y = c(0, 0), group = c("id", "ieb")), mapping = aes(x = x, y = y, label = group)) +
  ggrepel::geom_text_repel(data = filter(testing_res, fdr < 0.05), mapping = aes(x = coefficient, y = -log10(pvalue), label = target), size = 2, max.overlaps = 25, segment.size = 0.2) +
  scale_color_manual(values = c("black", "orange", "red")) + 
  ggthemes::theme_few() + 
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5)) + 
  labs(title = "Fibroblast: ieb vs. id") + 
  guides(color = guide_legend(override.aes = list(size = 5)))

## NOTES AND FOLLOW UP --------------------------------------------------------- 

# According to the referenced articles, along with the documentation for glht from multcomp, 
# the test I am performing is the Wald test. This is the same test that DESeq2 performs. It 
# is also used in glmmSeq. Specifically, I am letting glht decide whether or not to do the 
# approximate Wald test.

# Setting REML = FALSE in the modeling step (using maximum likelihood, not restricted maximum 
# likelihood) seems to be okay. There is debate on if using REML is okay for Wald tests, but 
# there does not seem to be resource that say you should not use ML. It might be wise to simply 
# stick with ML therefore. Apparently, you should NOT do REML if doing an LRT.

# In the end, since there are arguments to made for all tests, and the Wald test is used in some 
# widely-accepted packages, I believe that using it is okay. 

# The nebula package calculates a Chi-Square statistic when computing p values. It is an easy, fast 
# computation leading to the same result as doing an approximate Wald test. The Hotelling's T-Square 
# statistic is a possible substitute for the F test that can perform better in the multivariate use 
# case.

# Referenced articles: 
# https://pages.stat.wisc.edu/~ane/st572/notes/lec21.pdf
# https://sesug.org/proceedings/sesug_2023_final_papers/Statistics,_Analytics_and_Reporting/SESUG2023_Paper_122_Final_PDF.pdf#:~:text=Both%20approaches%2C%20Wald%20and%20F%2C%20test%20the,fidelity%20to%20the%20claimed%20false%20rejection%20rate.&text=REML%20estimation%20reparameterizes%20the%20data%20into%20linear,0%20and%20involve%20only%20the%20random%20effects.
# https://phillipalday.com/stats/mixedmodels.html
# https://www.reddit.com/r/AskStatistics/comments/fm06z6/restricted_vs_full_maximum_likelihood_estimation/
# https://cran.r-project.org/web/packages/multcomp/multcomp.pdf
# https://cran.r-project.org/web/packages/glmmSeq/vignettes/glmmSeq.html#Hypothesis_testing
# https://online.stat.psu.edu/stat505/lesson/7/7.1/7.1.3

# Some other R packages for fitting these sorts of models are: 
# - glmmTMB
# - INLA
# - GLMMadaptive

## BILLING ---------------------------------------------------------------------
# This entire process of learning this and doing the modeling took 12 hours
