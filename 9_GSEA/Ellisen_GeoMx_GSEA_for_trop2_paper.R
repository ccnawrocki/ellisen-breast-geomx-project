### GSEA ###
## Cole Nawrocki ##

set.seed(2001)
.libPaths()
# [1] "/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library"

# Packages 
library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)
#devtools::install_github("BarryDigby/enrichplot") # Someone made some improvements that I like
library(enrichplot)
library(ggplot2)
library(msigdbr)
library(BiocParallel)

theme_set(theme_bw())
theme_update(axis.title.y = element_blank(), 
             plot.title = element_text(hjust = 0.5), 
             axis.text.y = element_blank(), 
             axis.ticks.y = element_blank(), 
             panel.grid.minor.y = element_blank(), 
             panel.grid.major.y = element_blank(),
             panel.grid.major.x = element_line(color = "black", linewidth = 0.05),
             panel.grid.minor.x = element_line(color = "black", linewidth = 0.05))

# Gene sets
hallmarkgenesets <- msigdbr(species = "Homo sapiens", category = "H") |> 
  dplyr::select(gs_name, gene_symbol)
c2genesets <- msigdbr(species = "Homo sapiens", category = "C2") |> 
  dplyr::select(gs_name, gene_symbol)
c5genesets <- msigdbr(species = "Homo sapiens", category = "C5") |> 
  dplyr::select(gs_name, gene_symbol)

# Data 
fib_comp <- read.csv("4_DE_analysis/DE_analysis_for_trop2_paper/fibroblast id vs ii+ie limma_for_trop2_paper.csv", row.names = 1)
tum_comp <- read.csv("4_DE_analysis/comparisons_first_cohort/tumor id vs ii+ie limma.csv", row.names = 1)

#pdf(file = "trop2_paper_gsea_results.pdf", width = 10, height = 8)

## GENE SETS H, C2, and C5 -----------------------------------------------------
# Tumor
tumgenelist <- arrange(tum_comp, desc(t)) |> pull(t)
names(tumgenelist) <- arrange(tum_comp, desc(t)) |> pull(target)

tum_hallmark <- GSEA(geneList = tumgenelist, TERM2GENE = hallmarkgenesets, eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
tum_hallmark@result$GeneRatio <- (substr(x = tum_hallmark@result$leading_edge, start = 6, stop = 7) |> as.numeric())/100
#write.csv(x = tum_hallmark |> as.data.frame(), file = "tumor_hot_vs_cold_gsea_hallmark.csv")
tum_hallmark@result$ID %<>% factor(levels = (arrange(tum_hallmark@result, desc(NES)) |> pull(ID)))

dat <- tum_hallmark@result |> mutate(updown = (sign(NES) > 0)) |> group_by(updown) |> top_n(n = -10, wt = p.adjust)
ggplot(data = dat) + 
  geom_point(mapping = aes(y = ID, x = NES, colour = p.adjust, size = GeneRatio)) + 
  geom_bar(mapping = aes(y = ID, x = NES), width = 0.05, color = NA, fill = "black", stat = "identity") +
  geom_text(data = dat |> filter(NES < 0),
            mapping = aes(x = 0.1, y = ID, label = ID),
            size = 2.5, hjust = 0) +
  geom_text(data = dat |> filter(NES > 0),
            mapping = aes(x = -0.1, y = ID, label = ID),
            size = 2.5, hjust = 1) +
  geom_vline(xintercept = 0) + 
  scale_color_viridis_c() +
  labs(title = "Tumor: Hallmark") 
ggsave(filename = "tumor_hallmark.pdf", width = 7, height = 5)

tum_c2 <- GSEA(geneList = tumgenelist, TERM2GENE = c2genesets, eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
tum_c2@result$GeneRatio <- (substr(x = tum_c2@result$leading_edge, start = 6, stop = 7) |> as.numeric())/100
#write.csv(x = tum_c2 |> as.data.frame(), file = "tumor_hot_vs_cold_gsea_c2.csv")
tum_c2@result$ID %<>% factor(levels = (arrange(tum_c2@result, desc(NES)) |> pull(ID)))
dat <- tum_c2@result |> mutate(updown = (sign(NES) > 0)) |> group_by(updown) |> top_n(n = -10, wt = p.adjust)
ggplot(data = dat) + 
  geom_point(mapping = aes(y = ID, x = NES, colour = p.adjust, size = GeneRatio)) + 
  geom_bar(mapping = aes(y = ID, x = NES), width = 0.05, color = NA, fill = "black", stat = "identity") +
  geom_text(data = dat |> filter(NES < 0),
            mapping = aes(x = 0.1, y = ID, label = ID),
            size = 2, hjust = 0) +
  geom_text(data = dat |> filter(NES > 0),
            mapping = aes(x = -0.1, y = ID, label = ID),
            size = 2, hjust = 1) +
  geom_vline(xintercept = 0) + 
  scale_color_viridis_c() + 
  theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = "Tumor: C2") 
ggsave(filename = "tumor_c2.pdf", width = 10, height = 5)

tum_c5 <- GSEA(geneList = tumgenelist, TERM2GENE = c5genesets, eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
tum_c5@result$GeneRatio <- (substr(x = tum_c5@result$leading_edge, start = 6, stop = 7) |> as.numeric())/100
#write.csv(x = tum_c5 |> as.data.frame(), file = "tumor_hot_vs_cold_gsea_c5.csv")
tum_c5@result$ID %<>% factor(levels = (arrange(tum_c5@result, desc(NES)) |> pull(ID)))
dat <- tum_c5@result |> mutate(updown = (sign(NES) > 0)) |> group_by(updown) |> top_n(n = -10, wt = p.adjust)
ggplot(data = dat) + 
  geom_point(mapping = aes(y = ID, x = NES, colour = p.adjust, size = GeneRatio)) + 
  geom_bar(mapping = aes(y = ID, x = NES), width = 0.05, color = NA, fill = "black", stat = "identity") +
  geom_text(data = dat |> filter(NES < 0),
            mapping = aes(x = 0.1, y = ID, label = ID),
            size = 2, hjust = 0) +
  geom_text(data = dat |> filter(NES > 0),
            mapping = aes(x = -0.1, y = ID, label = ID),
            size = 2, hjust = 1) +
  geom_vline(xintercept = 0) + 
  scale_color_viridis_c() + 
  theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = "Tumor: C5") 
ggsave(filename = "tumor_c5.pdf", width = 8, height = 5)

# Fibroblast
fibgenelist <- arrange(fib_comp, desc(t)) |> pull(t)
names(fibgenelist) <- arrange(fib_comp, desc(t)) |> pull(target)

fib_hallmark <- GSEA(geneList = fibgenelist, TERM2GENE = hallmarkgenesets, eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
fib_hallmark@result$GeneRatio <- (substr(x = fib_hallmark@result$leading_edge, start = 6, stop = 7) |> as.numeric())/100
#write.csv(x = fib_hallmark |> as.data.frame(), file = "fibroblast_hot_vs_cold_gsea_hallmark.csv")
fib_hallmark@result$ID %<>% factor(levels = (arrange(fib_hallmark@result, desc(NES)) |> pull(ID)))
dat <- fib_hallmark@result |> mutate(updown = (sign(NES) > 0)) |> group_by(updown) |> top_n(n = -10, wt = p.adjust)
ggplot(data = dat) + 
  geom_point(mapping = aes(y = ID, x = NES, colour = p.adjust, size = GeneRatio)) + 
  geom_bar(mapping = aes(y = ID, x = NES), width = 0.05, color = NA, fill = "black", stat = "identity") +
  geom_text(data = dat |> filter(NES < 0),
            mapping = aes(x = 0.1, y = ID, label = ID),
            size = 2.5, hjust = 0) +
  geom_text(data = dat |> filter(NES > 0),
            mapping = aes(x = -0.1, y = ID, label = ID),
            size = 2.5, hjust = 1) +
  geom_vline(xintercept = 0) + 
  scale_color_viridis_c() + 
  theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = "Fibroblast: Hallmark") 
ggsave(filename = "fibroblast_hallmark.pdf", width = 8, height = 5)

fib_c2 <- GSEA(geneList = fibgenelist, TERM2GENE = c2genesets, eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
fib_c2@result$GeneRatio <- (substr(x = fib_c2@result$leading_edge, start = 6, stop = 7) |> as.numeric())/100
#write.csv(x = fib_c2 |> as.data.frame(), file = "fibroblast_hot_vs_cold_gsea_c2.csv")
fib_c2@result$ID %<>% factor(levels = (arrange(fib_c2@result, desc(NES)) |> pull(ID)))
dat <- fib_c2@result |> mutate(updown = (sign(NES) > 0)) |> group_by(updown) |> top_n(n = -10, wt = p.adjust)
ggplot(data = dat) + 
  geom_point(mapping = aes(y = ID, x = NES, colour = p.adjust, size = GeneRatio)) + 
  geom_bar(mapping = aes(y = ID, x = NES), width = 0.05, color = NA, fill = "black", stat = "identity") +
  geom_text(data = dat |> filter(NES < 0),
            mapping = aes(x = 0.1, y = ID, label = ID),
            size = 2.5, hjust = 0) +
  geom_text(data = dat |> filter(NES > 0),
            mapping = aes(x = -0.1, y = ID, label = ID),
            size = 2.5, hjust = 1) +
  geom_vline(xintercept = 0) + 
  scale_color_viridis_c() + 
  theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = "Fibroblast: C2") 
ggsave(filename = "fibroblast_c2.pdf", width = 8.5, height = 5)

fib_c5 <- GSEA(geneList = fibgenelist, TERM2GENE = c5genesets, eps = 1e-100, nPermSimple=10000, pvalueCutoff = 0.05, pAdjustMethod = "BH")
fib_c5@result$GeneRatio <- (substr(x = fib_c5@result$leading_edge, start = 6, stop = 7) |> as.numeric())/100
#write.csv(x = fib_c5 |> as.data.frame(), file = "fibroblast_hot_vs_cold_gsea_c5.csv")
fib_c5@result$ID <- ifelse(test = fib_c5@result$ID == "GOBP_ADAPTIVE_IMMUNE_RESPONSE_BASED_ON_SOMATIC_RECOMBINATION_OF_IMMUNE_RECEPTORS_BUILT_FROM_IMMUNOGLOBULIN_SUPERFAMILY_DOMAINS", 
                                    yes = "GOBP_ADAPTIVE_IMMUNE_RESPONSE (shortened)", no = fib_c5@result$ID)
fib_c5@result$ID %<>% factor(levels = (arrange(fib_c5@result, desc(NES)) |> pull(ID)))
dat <- fib_c5@result |> mutate(updown = (sign(NES) > 0)) |> group_by(updown) |> top_n(n = -10, wt = p.adjust)
ggplot(data = dat) + 
  geom_point(mapping = aes(y = ID, x = NES, colour = p.adjust, size = GeneRatio)) + 
  geom_bar(mapping = aes(y = ID, x = NES), width = 0.05, color = NA, fill = "black", stat = "identity") +
  geom_text(data = dat |> filter(NES < 0),
            mapping = aes(x = 0.1, y = ID, label = ID),
            size = 2.5, hjust = 0) +
  geom_text(data = dat |> filter(NES > 0),
            mapping = aes(x = -0.1, y = ID, label = ID),
            size = 2.5, hjust = 1) +
  geom_vline(xintercept = 0) + 
  scale_color_viridis_c() + 
  theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = "Fibroblast: C5")
ggsave(filename = "fibroblast_c5.pdf", width = 8, height = 5)

## GENE SETS BP and MF ---------------------------------------------------------
# Tumor
tumbpout <- clusterProfiler::gseGO(keyType = "SYMBOL", 
                                OrgDb = "org.Hs.eg.db", 
                                ont = "BP", 
                                geneList = tumgenelist, 
                                pAdjustMethod = "BH", 
                                pvalueCutoff = 0.05,
                                nPermSimple = 10000, 
                                by = "fgsea", 
                                eps = 1e-100, 
                                verbose = T)
tumbpout@result$GeneRatio <- (substr(x = tumbpout@result$leading_edge, start = 6, stop = 7) |> as.numeric())/100
#write.csv(x = tumbpout |> as.data.frame(), file = "tumor_hot_vs_cold_gsea_bp.csv")
tumbpout@result$Description %<>% factor(levels = (arrange(tumbpout@result, desc(NES)) |> pull(Description)))
dat <- tumbpout@result |> mutate(updown = (sign(NES) > 0)) |> group_by(updown) |> top_n(n = -10, wt = p.adjust)
ggplot(data = dat) + 
  geom_point(mapping = aes(y = Description, x = NES, colour = p.adjust, size = GeneRatio)) + 
  geom_bar(mapping = aes(y = Description, x = NES), width = 0.05, color = NA, fill = "black", stat = "identity") +
  geom_text(data = dat |> filter(NES < 0),
            mapping = aes(x = 0.1, y = Description, label = Description),
            size = 3, hjust = 0) +
  geom_text(data = dat |> filter(NES > 0),
            mapping = aes(x = -0.1, y = Description, label = Description),
            size = 3, hjust = 1) +
  geom_vline(xintercept = 0) + 
  scale_color_viridis_c() + 
  theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = "Tumor: BP")
ggsave("tumor_bp.pdf", width = 8.5, height = 5)

tummfout <- clusterProfiler::gseGO(keyType = "SYMBOL", 
                                   OrgDb = "org.Hs.eg.db", 
                                   ont = "MF", 
                                   geneList = tumgenelist, 
                                   pAdjustMethod = "BH", 
                                   pvalueCutoff = 0.05,
                                   nPermSimple = 10000, 
                                   by = "fgsea", 
                                   eps = 1e-100, 
                                   verbose = T)
tummfout@result$GeneRatio <- (substr(x = tummfout@result$leading_edge, start = 6, stop = 7) |> as.numeric())/100
#write.csv(x = tummfout |> as.data.frame(), file = "tumor_hot_vs_cold_gsea_bp.csv")
tummfout@result$Description <- ifelse(test = tummfout@result$ID == "GO:0016709", yes = "oxidoreductase activity (shortened)", no = tummfout@result$Description)
tummfout@result$Description %<>% factor(levels = (arrange(tummfout@result, desc(NES)) |> pull(Description)))
dat <- tummfout@result |> mutate(updown = (sign(NES) > 0)) |> group_by(updown) |> top_n(n = -10, wt = p.adjust)
ggplot(data = dat) + 
  geom_point(mapping = aes(y = Description, x = NES, colour = p.adjust, size = GeneRatio)) + 
  geom_bar(mapping = aes(y = Description, x = NES), width = 0.05, color = NA, fill = "black", stat = "identity") +
  geom_text(data = dat |> filter(NES < 0),
            mapping = aes(x = 0.1, y = Description, label = Description),
            size = 3, hjust = 0) +
  geom_text(data = dat |> filter(NES > 0),
            mapping = aes(x = -0.1, y = Description, label = Description),
            size = 3, hjust = 1) +
  geom_vline(xintercept = 0) + 
  scale_color_viridis_c() + 
  theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = "Tumor: MF")
ggsave("tumor_mf.pdf", width = 7.5, height = 5)

# Fibroblast
fibbpout <- clusterProfiler::gseGO(keyType = "SYMBOL", 
                                   OrgDb = "org.Hs.eg.db", 
                                   ont = "BP", 
                                   geneList = fibgenelist, 
                                   pAdjustMethod = "BH", 
                                   pvalueCutoff = 0.05,
                                   nPermSimple = 10000, 
                                   by = "fgsea", 
                                   eps = 1e-100, 
                                   verbose = T)
fibbpout@result$GeneRatio <- (substr(x = fibbpout@result$leading_edge, start = 6, stop = 7) |> as.numeric())/100
#write.csv(x = fibbpout |> as.data.frame(), file = "fibroblast_hot_vs_cold_gsea_bp.csv")
fibbpout@result$Description <- ifelse(test = fibbpout@result$Description == "regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains", 
                                      yes = "regulation of adaptive immune response (shortened)", no = fibbpout@result$Description)
fibbpout@result$Description <- ifelse(test = fibbpout@result$Description == "adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains", 
                                      yes = "adaptive immune response (shortened)", no = fibbpout@result$Description)
fibbpout@result$Description %<>% factor(levels = (arrange(fibbpout@result, desc(NES)) |> pull(Description)))
dat <- fibbpout@result |> mutate(updown = (sign(NES) > 0)) |> group_by(updown) |> top_n(n = -10, wt = p.adjust)
ggplot(data = dat) + 
  geom_point(mapping = aes(y = Description, x = NES, colour = p.adjust, size = GeneRatio)) + 
  geom_bar(mapping = aes(y = Description, x = NES), width = 0.05, color = NA, fill = "black", stat = "identity") +
  geom_text(data = dat |> filter(NES < 0),
            mapping = aes(x = 0.1, y = Description, label = Description),
            size = 3, hjust = 0) +
  geom_text(data = dat |> filter(NES > 0),
            mapping = aes(x = -0.1, y = Description, label = Description),
            size = 3, hjust = 1) +
  geom_vline(xintercept = 0) + 
  scale_color_viridis_c() + 
  theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = "Fibroblast: BP")
ggsave("fibroblast_bp.pdf", width = 7.5, height = 5)

fibmfout <- clusterProfiler::gseGO(keyType = "SYMBOL", 
                                   OrgDb = "org.Hs.eg.db", 
                                   ont = "MF", 
                                   geneList = fibgenelist, 
                                   pAdjustMethod = "BH", 
                                   pvalueCutoff = 0.05,
                                   nPermSimple = 10000, 
                                   by = "fgsea", 
                                   eps = 1e-100, 
                                   verbose = T)
fibmfout@result$GeneRatio <- (substr(x = fibmfout@result$leading_edge, start = 6, stop = 7) |> as.numeric())/100
#write.csv(x = fibmfout |> as.data.frame(), file = "fibroblast_hot_vs_cold_gsea_bp.csv")
fibmfout@result$Description %<>% factor(levels = (arrange(fibmfout@result, desc(NES)) |> pull(Description)))
dat <- fibmfout@result |> mutate(updown = (sign(NES) > 0)) |> group_by(updown) |> top_n(n = -10, wt = p.adjust)
ggplot(data = dat) + 
  geom_point(mapping = aes(y = Description, x = NES, colour = p.adjust, size = GeneRatio)) + 
  geom_bar(mapping = aes(y = Description, x = NES), width = 0.025, color = NA, fill = "black", stat = "identity") +
  geom_text(data = dat |> filter(NES < 0),
            mapping = aes(x = 0.1, y = Description, label = Description),
            size = 3, hjust = 0) +
  geom_text(data = dat |> filter(NES > 0),
            mapping = aes(x = -0.1, y = Description, label = Description),
            size = 3, hjust = 1) +
  geom_vline(xintercept = 0) + 
  scale_color_viridis_c() + 
  theme(axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
  labs(title = "Fibroblast: MF")
ggsave("fibroblast_mf.pdf", width = 7, height = 5)

# One GSEA plot for example
gseaplot2(tum_hallmark, 
         geneSetID = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 
         title = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 
         subplots = 1:2, pvalue_table = T)

# All GSEA plots
tumbpout@result$Description %<>% as.character()

for (setgroup in c("tum_hallmark", "tum_c2", "tum_c5", "fib_hallmark", "fib_c2", "fib_c5")) {
  dat <- get(setgroup)
  dat@result$Description %<>% as.character()
  toplot <- dat@result |> mutate(updown = (sign(NES) > 0)) |> group_by(updown) |> top_n(n = -10, wt = p.adjust) |> pull(Description)
  cat(setgroup)
  pdf(file = paste(setgroup, "gsea_plots.pdf", sep = "_"), width = 7, height = 5)
  for (descr in toplot) {
    print(
      gseaplot2(dat, 
                geneSetID = descr, 
                title = descr, 
                subplots = 1:2, pvalue_table = T)
    )
  }
  dev.off()
  cat("... done\n")
}

for (setgroup in c("tumbpout", "tummfout", "fibbpout", "fibmfout")) {
  dat <- get(setgroup)
  dat@result$Description %<>% as.character()
  toplot <- dat@result |> mutate(updown = (sign(NES) > 0)) |> group_by(updown) |> top_n(n = -10, wt = p.adjust) |> pull(ID)
  cat(setgroup)
  pdf(file = paste(setgroup, "gsea_plots.pdf", sep = "_"), width = 7, height = 5)
  for (id in toplot) {
    descr <- (dat@result)[id, "Description"]
    print(
      gseaplot2(dat, 
                geneSetID = id, 
                title = descr, 
                subplots = 1:2, pvalue_table = T)
    )
  }
  dev.off()
  cat("... done\n")
}

# ~ 10 hours

sessionInfo()
# R version 4.4.2 (2024-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.3.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] BiocParallel_1.40.0    enrichplot_1.26.6      org.Hs.eg.db_3.20.0    AnnotationDbi_1.68.0   IRanges_2.40.1        
# [6] S4Vectors_0.44.0       Biobase_2.66.0         BiocGenerics_0.52.0    ggplot2_3.5.1          dplyr_1.1.4           
# [11] msigdbr_7.5.1          clusterProfiler_4.14.4
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1        viridisLite_0.4.2       farver_2.1.2            blob_1.2.4              R.utils_2.12.3         
# [6] Biostrings_2.74.1       lazyeval_0.2.2          fastmap_1.2.0           digest_0.6.37           lifecycle_1.0.4        
# [11] KEGGREST_1.46.0         tidytree_0.4.6          RSQLite_2.3.9           magrittr_2.0.3          compiler_4.4.2         
# [16] rlang_1.1.5             tools_4.4.2             igraph_2.1.4            data.table_1.16.4       ggtangle_0.0.6         
# [21] labeling_0.4.3          bit_4.5.0.1             gson_0.1.0              plyr_1.8.9              RColorBrewer_1.1-3     
# [26] aplot_0.2.4             babelgene_22.9          withr_3.0.2             purrr_1.0.4             R.oo_1.27.0            
# [31] grid_4.4.2              GOSemSim_2.32.0         colorspace_2.1-1        GO.db_3.20.0            scales_1.3.0           
# [36] cli_3.6.3               crayon_1.5.3            treeio_1.30.0           generics_0.1.3          rstudioapi_0.17.1      
# [41] ggtree_3.14.0           httr_1.4.7              reshape2_1.4.4          DBI_1.2.3               qvalue_2.38.0          
# [46] ape_5.8-1               cachem_1.1.0            DOSE_4.0.0              stringr_1.5.1           zlibbioc_1.52.0        
# [51] splines_4.4.2           parallel_4.4.2          ggplotify_0.1.2         XVector_0.46.0          yulab.utils_0.2.0      
# [56] vctrs_0.6.5             Matrix_1.7-1            jsonlite_1.8.9          gridGraphics_0.5-1      patchwork_1.3.0        
# [61] bit64_4.6.0-1           ggrepel_0.9.6           tidyr_1.3.1             glue_1.8.0              codetools_0.2-20       
# [66] cowplot_1.1.3           stringi_1.8.4           gtable_0.3.6            GenomeInfoDb_1.42.3     UCSC.utils_1.2.0       
# [71] munsell_0.5.1           tibble_3.2.1            pillar_1.10.1           fgsea_1.32.2            GenomeInfoDbData_1.2.13
# [76] R6_2.5.1                lattice_0.22-6          R.methodsS3_1.8.2       png_0.1-8               memoise_2.0.1          
# [81] ggfun_0.1.8             Rcpp_1.0.14             fastmatch_1.1-6         nlme_3.1-166            fs_1.6.5               
# [86] pkgconfig_2.0.3        
