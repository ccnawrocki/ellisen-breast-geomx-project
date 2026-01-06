# Project Environment Setup
# This should work on an ARM mac, but there may be issues on other devices.

# These are the typical packages that I install for GeoMx analysis, but not all 
# are necessarily needed.

renv::init()
.libPaths()
# [1] "/Users/ccn22/Library/CloudStorage/OneDrive-MassGeneralBrigham/geomx/analysis/ellisen-breast/11_stuff_for_presentation/renv/library/macos/R-4.4/aarch64-apple-darwin20"
# [2] "/Users/ccn22/Library/Caches/org.R-project.R/R/renv/sandbox/macos/R-4.4/aarch64-apple-darwin20/f7156815"              

# Should work for any device
renv::install(packages = list("BiocManager", "lme4", "lmerTest"))
renv::install(packages = list("Bioc::ComplexHeatmap", "Bioc::limma", "Bioc::edgeR", "Bioc::DESeq2", "Bioc::glmGamPoi", "ggdendro", "ggh4x", "bioc::sparseMatrixStats"))
install.packages("Matrix", type = "source")
install.packages("irlba", type = "source")
install.packages("remotes")
renv::install(packages = list("ggpubr", "jpeg", "openxlsx", "R.utils", "emmeans"))

# NanoString packages
renv::install("bioc::GeomxTools")
renv::install("bioc::NanoStringNCTools")
renv::install("bioc::GeoMxWorkflows")

# Some stuff may depend on your compiler. I used the homebrew gcc, evidenced here: 
system("which gfortran")
# /opt/homebrew/bin/gfortran

# Basically, I needed to edit some files behind the scenes:
Sys.setenv(
  LDFLAGS = "-L/opt/homebrew/lib -L/opt/homebrew/opt/gettext/lib",
  CPPFLAGS = "-I/opt/homebrew/include"
)

# Create ~/.R directory if it doesn't exist
# dir.create("~/.R", showWarnings = FALSE)

# Edit the Makevars: 
# file.create('~/.R/Makevars')

# Opened it and added these lines:
# FC = /opt/homebrew/Cellar/gcc/14.2.0_1/bin/gfortran
# F77 = /opt/homebrew/Cellar/gcc/14.2.0_1/bin/gfortran
# FLIBS = -L/opt/homebrew/Cellar/gcc/14.2.0_1/lib/gcc/14
# LDFLAGS=-L/opt/homebrew/lib -L/opt/homebrew/opt/gettext/lib
# CPPFLAGS=-I/opt/homebrew/opt/libpng/include

# See https://mac.r-project.org/tools/ for gcc install not from homebrew. 
# I believe that you can install it and specify the corresponding paths in the Makevars file.

# Now install
renv::install("gdtools")
renv::install("bioc::GeoDiff")

# Stuff for GSEA
renv::install("bioc::clusterProfiler")
renv::install("bioc::msigdbr")
renv::install("bioc::org.Hs.eg.db")
renv::install("BarryDigby/enrichplot")

# Stuff for deconvolution
renv::install("bioc::SpatialDecon")
renv::install("Seurat")

# Random stuff
renv::install("ggbiplot")
renv::install("multcomp")
renv::install("pheatmap")
renv::install("bioc::ComplexHeatmap")
renv::install("viridis")
renv::install("ggthemes", "ggprism")

# If all of these packages are installed, you should not experience any problems. 
# Installing other packages should be pretty easy.
renv::snapshot()

# Restart R as a precaution.
.rs.restartR(clean = T)

