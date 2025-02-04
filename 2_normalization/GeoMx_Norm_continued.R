## GeoMx Normalization Continued ##

# I Only foresaw needing the log2(q3 + 1) data. But immune deconvolution
# uses the non-log2-transformed data, so extracting that here. 

library(GeomxTools)
library(NanoStringNCTools)
library(GeoMxWorkflows)
library(dplyr)

# Data 
obj <- readRDS("geomx_obj_processed.RDS")

# q3 norm, from the GeomxTools package
q3 <- obj@assayData$q_norm
q3[1:4,1:4]
#       DSP-1001660011137-E-A02.dcc DSP-1001660011137-E-A03.dcc DSP-1001660011137-E-A05.dcc DSP-1001660011137-E-A06.dcc
# A2M                      1.562981                    2.344471                   39.856013                    2.344471
# ACADM                    4.688943                    7.033414                    7.033414                    4.688943
# ACADS                    4.688943                    4.688943                    7.033414                    0.000000
# ACAT1                    1.562981                    0.000000                    2.344471                    4.688943

# Making sure that I can recreate it, so I know what is going on
q3s <- obj@assayData$exprs |> apply(MARGIN = 2, FUN = quantile, 0.75)
nfs <- q3s / EnvStats::geoMean(q3s)
q3_homemade <- (obj@assayData$exprs |> sweep(MARGIN = 2, STATS = nfs, FUN = "/"))
q3_homemade[1:4,1:4]
#       DSP-1001660011137-E-A02.dcc DSP-1001660011137-E-A03.dcc DSP-1001660011137-E-A05.dcc DSP-1001660011137-E-A06.dcc
# A2M                      1.562981                    2.344471                   39.856013                    2.344471
# ACADM                    4.688943                    7.033414                    7.033414                    4.688943
# ACADS                    4.688943                    4.688943                    7.033414                    0.000000
# ACAT1                    1.562981                    0.000000                    2.344471                    4.688943

# Looks good
(q3 == q3_homemade) |> all()
# [1] TRUE

# Saving
write.csv(q3_homemade, "q3norm_nonlog.txt")




