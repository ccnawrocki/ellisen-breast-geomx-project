###### GeoMx Data QC Script ######
## Cole Nawrocki ##

# This script is adapted from the following vignette from NanoString: 
# https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html

# Check that we are in the correct environment
.libPaths()
# [1] "/Users/ccn22/micromamba/envs/geomx-env/lib/R/library"

# Packages
library(GeomxTools)
library(NanoStringNCTools)
library(GeoMxWorkflows)
library(dplyr)
library(stringr)
library(ggplot2)
library(knitr)
library(reshape2)
library(cowplot)

# GeoMx object
geomx_obj <- readRDS("geomx_obj_unprocessed.RDS")
geomx_obj@phenoData@data$Run <- str_split(string = geomx_obj@phenoData@data$`slide name`, pattern = "_", simplify = T)[,2]

### QC ###
# QC parameters
QC_params <-
  list(minSegmentReads = 5000, # Minimum number of reads (1000) --> default
       percentTrimmed = 90,    # Minimum % of reads trimmed (80%) --> default
       percentStitched = 90,   # Minimum % of reads stitched (80%) --> default
       percentAligned = 90,    # Minimum % of reads aligned (80%) --> default
       percentSaturation = 90, # Minimum sequencing saturation (90%) --> higher than default
       minNegativeCount = 0,   # Minimum negative control counts (0) --> lower than default
       maxNTCCount = 10000,     # Maximum counts observed in NTC well (10000) --> higher than default
       minNuclei = 10,         # Minimum # of nuclei estimated (10) --> much lower than default
       minArea = 500)         # Minimum segment area (500) --> much lower than default

# Flagging
geomx_obj <-
  setSegmentQCFlags(geomx_obj, 
                    qcCutoffs = QC_params)        

# Collate QC Results
QCResults <- protocolData(geomx_obj)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

# Graphical summaries of QC statistics
col_by <- "Run"
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = col_by)) +
    geom_histogram(bins = 50, alpha = 0.5, position = "identity") +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

QC_histogram(sData(geomx_obj), "Trimmed (%)", "segment", 90)
QC_histogram(sData(geomx_obj), "Stitched (%)", "segment", 90)
QC_histogram(sData(geomx_obj), "Aligned (%)", "segment", 90)
QC_histogram(sData(geomx_obj), "Saturated (%)", "segment", 90) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")
QC_histogram(sData(geomx_obj), "area", "segment", 500, scale_trans = "log10")
QC_histogram(sData(geomx_obj), "nuclei", "segment", 10)

# Calculate the negative geometric means for each module
pkcs <- annotation(geomx_obj)
modules <- gsub(".pkc", "", pkcs)
negativeGeoMeans <- 
  esBy(negativeControlSubset(geomx_obj), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(geomx_obj)[["NegGeoMean"]] <- negativeGeoMeans

# Explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(geomx_obj)[, negCols] <- sData(geomx_obj)[["NegGeoMean"]]

# Show the Negative geoMeans plot for each module
for(ann in negCols) {
  plt <- QC_histogram(pData(geomx_obj), ann, "segment", 2, scale_trans = "log10")
  print(plt)
}

# detatch neg_geomean columns ahead of aggregateCounts call
pData(geomx_obj) <- pData(geomx_obj)[, !colnames(pData(geomx_obj)) %in% negCols]

# show all NTC values, Freq = # of Segments with a given NTC count:
kable(table(NTC_Count = sData(geomx_obj)$NTC),
      col.names = c("NTC Count", "# of Segments"))

# Show the QC summary
kable(QC_Summary, caption = "QC Summary Table for each Segment")

# Filtering flagged segments out 
geomx_obj <- geomx_obj[, QCResults$QCStatus == "PASS"]
dim(geomx_obj)
# Features  Samples 
# 18815      177 

# Removing poor probes
# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
geomx_obj <- setBioProbeQCFlags(geomx_obj, 
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = TRUE)

ProbeQCResults <- fData(geomx_obj)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))
qc_df
# Passed Global Local
# 1  18797      0    18

# Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
  subset(geomx_obj, 
         fData(geomx_obj)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(geomx_obj)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)
# Features  Samples 
# 18815      177 
geomx_obj <- ProbeQCPassed 

# Check how many unique targets the object has
length(unique(featureData(geomx_obj)[["TargetName"]]))

# collapse to gene targets
bygene_geomx_obj <- aggregateCounts(object = geomx_obj)
dim(bygene_geomx_obj)

exprs(bygene_geomx_obj)[1:5, 1:2]

# Define LOQ SD threshold and minimum value
cutoff <- 1
minLOQ <- 1

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(bygene_geomx_obj))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(bygene_geomx_obj)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(bygene_geomx_obj)[, vars[1]] * 
             pData(bygene_geomx_obj)[, vars[2]] ^ cutoff)
  }
}
pData(bygene_geomx_obj)$LOQ <- LOQ

# Filtering based on LOQ
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(bygene_geomx_obj)$Module == module
  Mat_i <- t(esApply(bygene_geomx_obj[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(bygene_geomx_obj)$TargetName, ]

# Save detection rate information to pheno data
pData(bygene_geomx_obj)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(bygene_geomx_obj)$GeneDetectionRate <-
  pData(bygene_geomx_obj)$GenesDetected / nrow(bygene_geomx_obj)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(bygene_geomx_obj)$DetectionThreshold <- 
  cut(pData(bygene_geomx_obj)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(pData(bygene_geomx_obj),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = Run)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Run")
ggplot(pData(bygene_geomx_obj),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = segment)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segement Type")

# cut percent genes detected at <10%
kable(table(pData(bygene_geomx_obj)$DetectionThreshold,
            pData(bygene_geomx_obj)$segment))
bygene_geomx_obj <-
  bygene_geomx_obj[, pData(bygene_geomx_obj)$GeneDetectionRate >= 0.1]
dim(bygene_geomx_obj)

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(bygene_geomx_obj)]
fData(bygene_geomx_obj)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(bygene_geomx_obj)$DetectionRate <-
  fData(bygene_geomx_obj)$DetectedSegments / nrow(pData(bygene_geomx_obj))

# Gene of interest detection table
goi <- c("PDCD1", "CD274", "IFNG", "CD8A", "CD68", "EPCAM",
         "KRT18", "NPHS1", "NPHS2", "CALB1", "CLDN8")
goi_df <- data.frame(
  Gene = goi,
  Number = fData(bygene_geomx_obj)[goi, "DetectedSegments"],
  DetectionRate = scales::percent(fData(bygene_geomx_obj)[goi, "DetectionRate"]))
goi_df

# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(bygene_geomx_obj)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(bygene_geomx_obj))
rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")

# Subset to target genes detected in at least 10% of the samples.
# Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(bygene_geomx_obj), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
bygene_geomx_obj <- 
  bygene_geomx_obj[fData(bygene_geomx_obj)$DetectionRate >= 0.1 |
                    fData(bygene_geomx_obj)$TargetName %in% neg_probes, ]
dim(bygene_geomx_obj)

# retain only detected genes of interest
goi <- goi[goi %in% rownames(bygene_geomx_obj)]

# Add the quantile-3-normalization to the object
bygene_geomx_obj <- GeomxTools::normalize(object = bygene_geomx_obj, 
                                          norm_method = "quant", 
                                          desiredQuantile = .75,
                                          toElt = "q_norm")

### Saving ###
# Saving the counts, the metadata, and the "processed" object
write.csv(bygene_geomx_obj@assayData$exprs, file = "counts.csv")
meta <- cbind(bygene_geomx_obj@phenoData@data |> select(-LOQ), 
              sData(bygene_geomx_obj) |> dplyr::select(roi, aoi))
rownames(meta) <- gsub(pattern = "-", replacement = ".", x = rownames(meta))
write.csv(meta, "meta.csv")
saveRDS(bygene_geomx_obj, file = "geomx_obj_processed.RDS")

### Session ###
sessionInfo()
