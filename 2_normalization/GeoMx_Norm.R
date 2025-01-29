###### GeoMx Data Normalization Script ######
## Cole Nawrocki ##

# Check that we are in the correct environment
.libPaths()
# [1] "/Users/ccn22/micromamba/envs/geomx-env/lib/R/library"

# Packages
library(dplyr)
library(stringr)
library(ggplot2)
library(ComplexHeatmap)
library(limma)
library(edgeR)
library(ggbiplot)
library(magrittr)
library(patchwork)
library(lme4)
library(lmerTest)

# Colors
segment_cols <- c("tumor" = "green4", "immune" = "red", "fibroblast" = "dodgerblue")
run_cols <- c("R1"="yellow", "R2"="orange")
slide_cols <- grDevices::colors()[grep('(gr(a|e)y)|(white)', grDevices::colors(), invert = T)] |> sample(size = 8, replace = F)
names(slide_cols) <- unique(meta$slide.name)

# Counts and metadata
cts <- read.csv("counts.csv", row.names = 1)
meta <- read.csv("meta.csv", row.names = 1)
all(rownames(meta) == colnames(cts))

### Q3 Normalization ###
qs <- apply(X = cts, MARGIN = 2, FUN = quantile, 0.75)
nfs <- qs / EnvStats::geoMean(qs)
q3norm <- sweep(x = cts, MARGIN = 2, STATS = nfs, FUN = "/")
lq3 <- log2(q3norm + 1)

# Identifying the 5% of genes with the highest variance
variances <- apply(X = lq3, MARGIN = 1, FUN = var)
top_genes <- names(variances)[variances >= quantile(variances, 0.95)]

# PCA
mat <- lq3[top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = meta$segment) + 
  theme_classic() 
ggbiplot(pcaout, var.axes = F, groups = meta$slide.name) + 
  theme_classic()
ggbiplot(pcaout, var.axes = F, groups = meta$Run) + 
  theme_classic()

# Heatmap with these genes
ha <- HeatmapAnnotation(segment = meta$segment, 
                        slide =  meta$slide.name, 
                        run = meta$Run,
                        col = list("segment"=segment_cols, 
                                   "slide"=slide_cols, 
                                   "run"=run_cols))
Heatmap(matrix = mat |> t(), 
        top_annotation = ha, 
        row_names_gp = gpar(fontsize = 3), 
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled\nQ3-Normalized\nCounts", 
        show_column_names = F, 
        show_row_names = F
)

# Fitting a linear mixed effects model to the data, with segment as the fixed effect and slide as the random effect.
# Then, determining the 5% of genes that have the smallest MSE values (based on Pearson residuals). 
meta$segment %<>% as.factor()
meta$slide.name %<>% as.factor()
lmem_pear_mse <- function(gene) {
  d <- dplyr::mutate(meta, expr=(lq3[gene,] |> unlist()))
  testout <- lme4::lmer(formula = expr ~ segment + (1+segment|slide.name), data = d)
  return((resid(testout, type = "pear")**2 |> sum()) / nrow(d))
}

mses <- parallel::mclapply(rownames(lq3), lmem_pear_mse, mc.cores = 8)
names(mses) <- rownames(lq3)
mses %<>% unlist()

top_genes <- names(mses)[mses <= quantile(mses, 0.05)]

# PCA
mat <- lq3[top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = meta$segment) + 
  theme_classic() 
ggbiplot(pcaout, var.axes = F, groups = meta$slide.name) + 
  theme_classic()
ggbiplot(pcaout, var.axes = F, groups = meta$Run) + 
  theme_classic()

# Heatmap with these genes
Heatmap(matrix = mat |> t(), 
        top_annotation = ha, 
        row_names_gp = gpar(fontsize = 3), 
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled\nQ3-Normalized\nCounts", 
        show_column_names = F, 
        show_row_names = F
)

# Saving
write.csv(x = lq3, file = "q3norm.csv")

### Normalization with limma ###
# Creating the object
cts <- DGEList(counts = cts, samples = meta)

# Calculating normalization factors
cts <- calcNormFactors(object = cts)

# Normalizing with voom
norm <- voom(counts = cts, plot = T)

# Fitting the limma model
limfit <- lmFit(norm)

# Identifing the 5% most "variable" genes via dispersion
disps <- cbind(limfit$coefficients, "sigma"=limfit$sigma) |> as.data.frame() 
disp_top_genes <- rownames(disps)[disps$sigma >= quantile(disps$sigma, 0.95)]

# PCA
mat <- norm$E[disp_top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = meta$slide.name) + 
  theme_classic() + labs(title = "Dispersion")
ggbiplot(pcaout, var.axes = F, groups = meta$segment) + 
  theme_classic() + labs(title = "Dispersion")
ggbiplot(pcaout, var.axes = F, groups = meta$Run) + 
  theme_classic() + labs(title = "Dispersion")

# Heatmap with these "variable" genes
Heatmap(matrix = mat |> t(), 
        top_annotation = ha, 
        row_names_gp = gpar(fontsize = 3), 
        #column_split = meta$segment, 
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled log2CPM", 
        show_row_names = F
)

# Refitting the limma model with a similar mixed model setup as above. 
mm <- model.matrix(~segment, data = meta)
cts <- DGEList(counts = cts)
cts <- calcNormFactors(object = cts)
norm <- voom(cts, design = mm, plot = T)
corfit <- duplicateCorrelation(norm, design = mm, block = meta$slide.name)
limfit <- lmFit(object = norm, design = mm, block = meta$slide.name, correlation = corfit$consensus)

# Identifing the 5% most "variable" genes via dispersion 
disps <- cbind(limfit$coefficients, "sigma"=limfit$sigma) |> as.data.frame() 
disp_top_genes <- rownames(disps)[disps$sigma >= quantile(disps$sigma, 0.95)]

# PCA
mat <- norm$E[disp_top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = meta$slide.name) + 
  theme_classic()
ggbiplot(pcaout, var.axes = F, groups = meta$segment) + 
  theme_classic()
ggbiplot(pcaout, var.axes = F, groups = meta$Run) + 
  theme_classic()

# Heatmap with these "variable" genes
Heatmap(matrix = mat |> t(), 
        top_annotation = ha, 
        row_names_gp = gpar(fontsize = 3), 
        #column_split = meta$segment, 
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Scaled log2CPM", 
        show_row_names = F
)

# Saving
write.csv(x = norm$E, file = "vnorm.csv")

### Batch Correction ###
# We will use limma for batch correction of the q3 normalized data by run. We will find variable genes with variance.
norm_bc <- limma::removeBatchEffect(x = lq3, 
                                    batch = meta$Run, 
                                    #group = meta$segment, 
                                    #covariates = meta$area
)

variances <- apply(X = norm_bc, MARGIN = 1, FUN = var)
top_genes <- names(variances)[variances >= quantile(variances, 0.95)]

# PCA
mat <- norm_bc[top_genes,] |> apply(1, scale)
pcaout <- prcomp(mat, center = F, scale. = F)
plot(pcaout)
pca_var<-pcaout$sdev**2
pve <- pca_var/sum(pca_var)
plot(pve)
ggbiplot(pcaout, var.axes = F, groups = meta$slide.name) + 
  theme_classic() + labs(title = "Dispersion")
ggbiplot(pcaout, var.axes = F, groups = meta$segment) + 
  theme_classic() + labs(title = "Dispersion")
ggbiplot(pcaout, var.axes = F, groups = meta$Run) + 
  theme_classic() + labs(title = "Dispersion")

# Heatmap with the "variable" genes
Heatmap(matrix = mat |> t(), 
        top_annotation = ha, 
        row_names_gp = gpar(fontsize = 3), 
        cluster_rows = T, 
        cluster_columns = T, 
        name = "Batch\nCorrected\nlog2CPM", 
        show_row_names = F
        )

# Saving
write.csv(x = norm_bc, file = "q3norm_bc.csv")

### Comparison ###
boxplot(norm_bc,
        col = "#9EDAE5", 
        main = "Batch Corrected Q3 Normalized Counts", 
        names = 1:155, xlab = "AOI #")
boxplot(norm$E,
        col = "green4", 
        main = "Voom Normalized Counts", 
        names = 1:155, xlab = "AOI #")
boxplot(lq3,
        col = "orange", 
        main = "Q3 Normalized Counts", 
        names = 1:155, xlab = "AOI #")
boxplot(cts$counts,
        col = "red3", 
        main = "Counts", 
        names = 1:155, xlab = "AOI #")

### Thoughts/Conclusions ###
# I believe that Q3 normalized counts without batch correction will be good. If we begin to suspect that there is a batch 
# effect in play, then we can use the Q3 normalized counts, corrected for run number. 
