# Project Information
Analysis completed on behalf of the MGH KFCCR Tumor Cartography Core for the Ellisen Lab. 
- Steps 0-8 require that all code is run in a micromamba environment made from the YML file provided in order to reproduce results exactly.
- Steps 9-10 were completed outside of this environment. Plots from Step 9 specifically will be near impossible to reproduce exactly, since I forgot to set the seed when I made the plots initially and sent them over for the paper.

## Steps
### 0_data_setup
Setting up the unprocessed GeoMx dataset object.

### 1_QC
Performing quality control and aggregating the counts among sets of probes that correspond to the same genes.

### 2_normalization
Q3-normalizing the data for downstream use with other NanoString tools, such as `SpatialDecon`, as well as exploring some alternative normalization methods.

### 3_segment_markers
Ensuring that the three tissue compartments (tumor, fibroblast, and immune) have marker genes that make sense as a sanity check.

### 4_DE_analysis
Performing differential expression analysis with `limma` and exploring other methods, such as mixed effects modeling with `lme4`.

### 5_immune_deconvolution
Deconvolving the immune AOIs to assess abundances of immune cell types, using `SpatialDecon`.

### 6_fibroblast_deconvolution
Deconvolving the fibroblast AOIs to assess abundances of CAF cell types, using `SpatialDecon`.

### 7_dim_reduction
Exploring some possible improvements in the gene selection and normalization processes for making PCA plots.

### 8_geodiff_exploration
Exploring the `GeoDiff` package, which offers an alternative QC and normalization workflow for this kind of data.

### 9_GSEA
Performing gene set enrichment analysis, using the `limma` results from Step 4. Publication figures: 
- TROP2 paper Figure S1C (as noted above, this figure may be difficult to reproduce exactly, due to failure to set the seed)

### 10_final_figures_and_methods
Creating some specific plots and writing up computational methods for publication. Publication figures: 
- TROP2 paper Figure 1C and 1E
- TROP2 paper Figure S1B

