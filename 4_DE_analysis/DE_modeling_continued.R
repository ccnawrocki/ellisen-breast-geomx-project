## DE Modeling Continued ##
### Cole Nawrocki ###

# This is the same process as the other modeling script, just for a simpler model that does not 
# include the core/border information.

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
meta$group <- paste(meta$AOI_code, meta$infiltration_category, sep = "_")
alldata <- inner_join(x = meta, y = norm, by = "aoi_id")
alldata$Patient_number %<>% as.factor()
alldata$group %<>% as.factor()
alldata <- tidyr::pivot_longer(data = alldata, cols = 17:ncol(alldata), names_to = "FeatureName", values_to = "y")
d <- alldata |> split(f = alldata$FeatureName)

## MODELING --------------------------------------------------------------------
# I will do the modeling for the entire dataset on eristwo, then I will load it back in
# I have a conda environment on eristwo called lme4-env that I activate in an interactive session

# Saving the data (it is already all set up)
saveRDS(d, file = "split_data_for_simple_modeling.RDS")

# Now, switching over to my terminal to use the HPC.

# $ scp /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/analysis/ellisen-breast/split_data_for_simple_modeling.RDS ccn22@eristwo.partners.org:/PHShome/ccn22/scratch/split_data_for_simple_modeling.RDS 
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
# d <- readRDS(file = "split_data_for_simple_modeling.RDS")
# 
# lmem_de_model <- function(data_df) {
#   lme4mod <- lme4::lmer(y ~ 1+group+(1+group|Patient_number), data_df, REML = F)
#   return(lme4mod)
# }
# 
# modeling_res <- parallel::mclapply(X = d, FUN = lmem_de_model, mc.cores = 35)
# saveRDS(modeling_res, file = "lmem_de_simple_model_results_ml.RDS")
# 
# q()

# $ conda deactivate
# $ exit
# $ exit
# $ scp ccn22@eristwo.partners.org:/PHShome/ccn22/scratch/lmem_de_simple_model_results_ml.RDS /Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/analysis/ellisen-breast/lmem_de_simple_model_results_ml.RDS 

# NOTE: I ran this once: using ML.

# It worked. 
full_res <- readRDS(file = "lmem_de_simple_model_results_ml.RDS")
object.size(full_res)
# 870199896 bytes
870199896/1e9
# [1] 0.8701999

