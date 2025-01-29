# GeoMx Linear Mixed Modeling with lme4 on an HPC
## Cole Nawrocki

# I uploaded the data to my scratch folder via scp.
# Then, I used an interactive job with 50 cpus and 50 GB of RAM.
# I had some trouble scheduling the job as non-interactive with the sbatch function.

# $ ssh -Y ccn22@eristwo.partners.org
# $ srun --pty -p interactive -c 50 --mem=50000MB -t 300 /bin/bash
# $ cd scratch
# $ conda activate lme4-env
# $ R

library(lme4)
library(lmerTest)
library(parallel)
library(magrittr)
options(mc.cores = 50)
set.seed(2001)

norm <- read.csv("ptnorm.csv", row.names = 1)
meta <- read.csv("meta_cleaned_v7.csv", row.names = 1)

norm <- t(norm) |> as.data.frame()
norm$aoi_id <- rownames(norm)
alldata <- dplyr::inner_join(x = meta, y = norm, by = "aoi_id")
alldata$Patient_number %<>% as.factor()
alldata <- tidyr::pivot_longer(data = alldata, cols = 20:ncol(alldata), names_to = "FeatureName", values_to = "y")
d <- alldata |> split(f = alldata$FeatureName)

lmem_de <- function(data_df) {
  lme4mod <- lmer(y ~ 1 + group + (1+group|Patient_number), data_df, REML = FALSE)
  covmat <- lme4mod@vcov_beta
  ests <- (summary(lme4mod) %$% coefficients)[,"Estimate"]
  return(list("est" = ests, "cov" = covmat))
}

demodres <- parallel::mclapply(X = d, FUN = lmem_de, mc.cores = 50)

saveRDS(demodres, file = "lmm_ptnorm.RDS")
