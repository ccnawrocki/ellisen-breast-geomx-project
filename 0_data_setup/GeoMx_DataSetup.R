###### GeoMx Data Setup Script ######
## Cole Nawrocki ##

# Check that we are in the correct environment
.libPaths()
# [1] "/Users/ccn22/micromamba/envs/geomx-env/lib/R/library"

# Packages
library(GeomxTools)
library(NanoStringNCTools)
library(GeoMxWorkflows)
library(data.table)
library(openxlsx)
library(dplyr)
library(stringr)
library(magrittr)

# Defining where the data is
datadir <- "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/ellisen-breast-data"

# The directory should have the lab worksheets and the GeoMx_NGS_Pipeline_DCC folders files in it
list.files(datadir) 
# [1] "BreastCa_R1_073124_20240905T1040_LabWorksheet.txt" "BreastCa_R2_080224_20240905T1039_LabWorksheet.txt"
# [3] "GeoMx_NGS_Pipeline_DCC_R1"                         "GeoMx_NGS_Pipeline_DCC_R2"                        
# [5] "Hs_R_NGS_WTA_v1.0.pkc"                             "images"                                           

# Create the metadata excel file. 
labworksheets <- list.files(datadir)[list.files(datadir) |> grep(pattern = "LabWorksheet")]
metas <- list()
for (sheet in labworksheets) {
  
  # Getting the run name via the worksheet
  run_name <- str_split(string = sheet, pattern = "_[0-9]{3}", simplify = T)[,1]
  
  # Extracting the header and the NTC line from the worksheet
  header <- readLines(con = file.path(datadir, sheet), n = 17) 
  
  # Defining the column names for the metadata
  column_names <- header[16] |> 
    str_split(pattern = "\t") |> 
    unlist()
  
  # Reading the metadata, without the header and the NTC line
  meta <- fread(file.path(datadir, sheet), 
                skip = 17, sep = "\t", header = F) |> 
    as.data.frame()
  
  # Adding the column names
  colnames(meta) <- column_names
  
  # Reformatting the Roi column
  meta$Roi <- gsub(x = meta$Roi, pattern = "=|\"", replacement = "") 
  
  # Adding the NTC line at the top of the metadata
  ntcline <- header[17] |> str_split(pattern = "\t") |> unlist()
  
  # It needs to be the same number of columns as the other lines
  ntcline <- c(ntcline, rep("", length(column_names)-2)) |>
    matrix(nrow = 1, byrow = T, dimnames = list("0", colnames(meta)))
  
  # It needs to have the panel information added to it
  ntcline[which(colnames(meta) == "Panel")] <- meta$Panel[1]
  
  # Actually combining the NTC line and the other metadata
  meta <- rbind(ntcline, meta)
  
  # Adding to the metadata list
  metas[[run_name]] <- meta
}

# Combining all the metadata data frames
metadata <- bind_rows(metas)

# The columns must be in this order
metadata %<>% dplyr::select(Sample_ID, Panel, `Slide Name`, 
                            Roi, Segment, Aoi, Area, Nuclei, 
                            `Scan Width`, `Scan Height`, 
                            `ROI Coordinate X`, `ROI Coordinate Y`, 
                            `Scan Offset X`, `Scan Offset Y`, Tags)

# Some columns must be renamed
metadata %<>% rename(panel=Panel, `slide name`=`Slide Name`, roi=Roi, aoi=Aoi, segment=Segment, nuclei=Nuclei, area=Area)


write.xlsx(metadata, file = "metadata.xlsx")

# Move all the dcc files into one "dcc" folder. Have to do this in bash.
# $ cd ~/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/ellisen-breast-data
# $ mkdir dccs
# $ cp GeoMx_NGS_Pipeline_DCC*/* dccs

# Creating the GeoMx objects
PKCFiles <- list.files(datadir)[list.files(datadir) |> grep(pattern = ".pkc")]
DCCFiles <- dir(file.path(datadir, "dccs"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
ANNOFile <- "metadata.xlsx"
geomx_obj <- readNanoStringGeoMxSet(dccFiles = DCCFiles,
                                    pkcFiles = file.path(datadir, PKCFiles),
                                    phenoDataFile = ANNOFile,
                                    phenoDataSheet = "Sheet 1",
                                    phenoDataDccColName = "Sample_ID",
                                    protocolDataColNames = c("aoi", "roi"),
                                    experimentDataColNames = c("panel")
                                    )

# Checking the data looks right
geomx_obj@assayData$exprs[1:5,1:5]
geomx_obj@phenoData@data |> dplyr::glimpse()

# Exporting the geomx object as an RDS files
saveRDS(object = geomx_obj, file = "geomx_obj_unprocessed.RDS")

# Session
sessionInfo()
