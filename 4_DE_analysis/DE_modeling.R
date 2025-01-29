## DE Modeling ##
### Cole Nawrocki ###

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

# Data
norm <- read.csv(file = "q3norm.csv", row.names = 1)
norm <- t(norm) |> as.data.frame()
norm$aoi_id <- rownames(norm)
meta <- openxlsx::read.xlsx("meta_cleaned_v5.xlsx")
meta$infiltration_group <- meta$infiltration_type
meta$infiltration_group <- ifelse(test = meta$infiltration_category == "desert", 
                                  yes = str_sub(meta$infiltration_group, start = 1, end = 2), 
                                  no = meta$infiltration_group)
meta$group <- paste(meta$AOI_code, meta$infiltration_group, sep = "_")
alldata <- inner_join(x = meta, y = norm, by = "aoi_id")
alldata$Patient_number %<>% as.factor()
alldata$group %<>% as.factor()
alldata <- tidyr::pivot_longer(data = alldata, cols = 18:ncol(alldata), names_to = "FeatureName", values_to = "y")
d <- alldata |> split(f = alldata$FeatureName)

## FIGURING OUT THE MODEL ------------------------------------------------------
# Just doing for BRCA1 as an example
dd <- d$BRCA1

# Example contrast
mm <- model.matrix(~ 1+group, data = meta)
contr <- (mm[meta$group == "tumor_iic", ] |> colMeans()) - (mm[meta$group == "tumor_iec", ] |> colMeans())

# Modeling with lmerTest
modelout <- lmerTest::lmer(y ~ 1+group+(1+group|Patient_number), dd, REML = T)
object.size(modelout)
# 187128 bytes

modelout |> glht(linfct = matrix(contr, nrow = 1), alternative = "two.sided", rhs = 0) |> summary(test = univariate()) %$%test %$%tstat[[1]]
# [1] 0.4343607

# Modeling with lme4
modelout <- lme4::lmer(y ~ 1+group+(1+group|Patient_number), dd, REML = T)
object.size(modelout)
# 67776 bytes

modelout |> glht(linfct = matrix(contr, nrow = 1), alternative = "two.sided", rhs = 0) |> summary(test = univariate()) %$%test %$%tstat[[1]]
# [1] 0.4343607

# We get the same result and the output is 3x smaller... for our purposes, using lme4 will be sufficient, so we will continue with that.
# Above, we used the Wald test and let glht decide which test to do. It did an approximate Wald test, using the z distribution. We can do other tests as follows:

# Making it do a student's t test. Since we have 1 contrast, the global F test is the same as doing a t test
modelout |> glht(linfct = matrix(contr, nrow = 1), alternative = "two.sided", rhs = 0) |> summary(test = Ftest())
# General Linear Hypotheses
# 
# Linear Hypotheses:
#        Estimate
# 1 == 0   0.1765
# 
# Global Test:
#        F DF1 DF2 Pr(>F)
# 1 0.1887   1  77 0.6652

# The same logic applies for an approximate Wald test (z test). Since  we have 1 contrast, the global Chi-Square test is the same as doing a z test.
modelout |> glht(linfct = matrix(contr, nrow = 1), alternative = "two.sided", rhs = 0) |> summary(test = Chisqtest())
# General Linear Hypotheses
# 
# Linear Hypotheses:
#        Estimate
# 1 == 0   0.1765
# 
# Global Test:
#    Chisq DF Pr(>Chisq)
# 1 0.1887  1      0.664

# For large n, the t test and the z test above are practically the same. univariate() decides for us if n is sufficiently large and does the indicated test. 
# But we can get fancier, if we want to...

# The emmeans package can test marginal means, using a t test too. It can approximate df with two methods. 
# Kenward-Roger is not always indicated for use with mixed models, and this method is computationally intense.
emmeansgrid <- modelout |> emmeans::emmeans(~ 1+group, lmer.df = "satterthwaite")
emmeans::contrast(object = emmeansgrid, list("test"=contr))
# contrast estimate    SE  df t.ratio p.value
# test        0.177 0.406 8.8   0.434  0.6745
# 
# Degrees-of-freedom method: satterthwaite 

emmeansgrid <- modelout |> emmeans::emmeans(~ 1+group, lmer.df = "kenward-roger")
emmeans::contrast(object = emmeansgrid, list("test"=contr))
# contrast estimate    SE   df t.ratio p.value
# test        0.177 0.932 4.85   0.189  0.8574
# 
# Degrees-of-freedom method: kenward-roger 

# The lmerTest package can test differences in means as well, using a t test. This can do both Kenward-Roger or Satterthwaite approximations of the degrees of freedom.
# Setting joint = TRUE will do a global test for many contrasts. So, joint = TRUE will do an F test, like above. Again, if we have 1 contrast, it does not matter.
# If we have multiple contrasts in a matrix and we do joint = FALSE, a t test will be done for each contrast.
contest(model = modelout, L = contr, rhs = 0, joint = F, confint = F, ddf = "Satterthwaite")
#   Estimate Std. Error       df   t value  Pr(>|t|)
# 1 0.176508  0.4063628 8.798733 0.4343607 0.6744833

contest(model = modelout, L = contr, rhs = 0, joint = F, confint = F, ddf = "Kenward-Roger")
#   Estimate Std. Error       df   t value  Pr(>|t|)
# 1 0.176508  0.9315828 4.846627 0.1894711 0.8574018

# I think that it makes sense to fit a model for every gene using lme4, with REML, since we will be performing Wald tests. 
# If we were going to perform LRTs, then REML cannot be used. We can decide on the exact test used downstream, but we will 
# extract results for each gene in a similar fashion to this:
test <- modelout |> 
  glht(linfct = matrix(contr, nrow = 1), alternative = "two.sided", rhs = 0) |> 
  summary(test = univariate()) %$%test
out <- data.frame(
  coefficient = test$coefficients[[1]],
  sigma = test$sigma[[1]],
  tstat = test$tstat[[1]],
  pvalue = test$pvalues[[1]], 
  test_type = test$type
)
out
#   coefficient     sigma     tstat    pvalue  test_type
# 1    0.176508 0.4063628 0.4343607 0.6640266 univariate

## FUNCTIONS FOR MODELING AND TESTING ------------------------------------------
# Function for doing the modeling
lmem_de_model <- function(data_df) {
  lme4mod <- lme4::lmer(y ~ 1+group+(1+group|Patient_number), data_df, REML = T)
  return(lme4mod)
}

# Function for doing the testing
lmem_de_test <- function(outmod, .contr, .test_type = "univariate") {
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

# Trying it out and it works
mm <- model.matrix(~ 1+group, data = meta)
contrast_list <- list(
  "tumor_iic / tumor_iec" = (mm[meta$group == "tumor_iic", ] |> colMeans()) - (mm[meta$group == "tumor_iec", ] |> colMeans())
)

modeling_res <- lapply(X = d[1:3], FUN = lmem_de_model)
testing_res <- lapply(X = modeling_res, FUN = lmem_de_test, .contr = "tumor_iic / tumor_iec", .test_type = "univariate")
testing_res |> bind_rows(.id = "target")
#   target              contrast coefficient     sigma      tstat    pvalue  test_type
# 1   A1BG tumor_iic / tumor_iec   0.6324106 0.4288017  1.4748325 0.1402576 univariate
# 2   A1CF tumor_iic / tumor_iec  -0.4274537 0.3669365 -1.1649255 0.2440491 univariate
# 3    A2M tumor_iic / tumor_iec   0.1203431 0.4697354  0.2561933 0.7978016 univariate

testing_res <- lapply(X = modeling_res, FUN = lmem_de_test, .contr = "tumor_iic / tumor_iec", .test_type = "Chisqtest")
testing_res |> bind_rows(.id = "target")
#   target              contrast coefficient     sigma df    fstat      tstat    pvalue test_type
# 1   A1BG tumor_iic / tumor_iec   0.6324106 0.4288017  1 2.175131  1.4748325 0.1402576     Chisq
# 2   A1CF tumor_iic / tumor_iec  -0.4274537 0.3669365  1 1.357051 -1.1649255 0.2440491     Chisq
# 3    A2M tumor_iic / tumor_iec   0.1203431 0.4697354  1 0.065635  0.2561933 0.7978016     Chisq

testing_res <- lapply(X = modeling_res, FUN = lmem_de_test, .contr = "tumor_iic / tumor_iec", .test_type = "Ftest")
testing_res |> bind_rows(.id = "target")
#   target              contrast coefficient     sigma df1 df2    fstat      tstat    pvalue test_type
# 1   A1BG tumor_iic / tumor_iec   0.6324106 0.4288017   1  77 2.175131  1.4748325 0.1443348         F
# 2   A1CF tumor_iic / tumor_iec  -0.4274537 0.3669365   1  77 1.357051 -1.1649255 0.2476455         F
# 3    A2M tumor_iic / tumor_iec   0.1203431 0.4697354   1  77 0.065635  0.2561933 0.7984848         F

testing_res <- lapply(X = modeling_res, FUN = lmem_de_test, .contr = "tumor_iic / tumor_iec", .test_type = "Satterthwaite")
testing_res |> bind_rows(.id = "target")
#   target              contrast   Estimate Std. Error        df    t value  Pr(>|t|)
# 1   A1BG tumor_iic / tumor_iec  0.6324106  0.4288017  7.777033  1.4748325 0.1795521
# 2   A1CF tumor_iic / tumor_iec -0.4274537  0.3669365 28.962626 -1.1649255 0.2535546
# 3    A2M tumor_iic / tumor_iec  0.1203431  0.4697354 46.621060  0.2561933 0.7989287

## MODELING --------------------------------------------------------------------
# I will do the modeling for the entire dataset on eristwo, then I will load it back in
# I have a conda environment on eristwo called lme4-env that I activate in an interactive session

# Saving the data (it is already all set up)
saveRDS(d, file = "split_data_for_modeling.RDS")

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
# d <- readRDS(file = "split_data_for_modeling.RDS")
# 
# lmem_de_model <- function(data_df) {
#   lme4mod <- lme4::lmer(y ~ 1+group+(1+group|Patient_number), data_df, REML = T)
#   return(lme4mod)
# }
# 
# modeling_res <- parallel::mclapply(X = d, FUN = lmem_de_model, mc.cores = 35)
# saveRDS(modeling_res, file = "lmem_de_model_results_reml.RDS")
# 
# q()

# $ conda deactivate
# $ exit
# $ exit
# $ scp ccn22@eristwo.partners.org:/PHShome/ccn22/scratch/lmem_de_model_results_reml.RDS /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/analysis/ellisen-breast/lmem_de_model_results_reml.RDS 

# NOTE: I ran this twice: once using REML and once using ML.

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
