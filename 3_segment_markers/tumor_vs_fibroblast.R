rm(list = ls())

meta <- read.csv("meta.csv", row.names = 1)
cts <- read.csv("counts.csv", row.names = 1)
all(rownames(meta) == colnames(cts))

idx <- meta$segment %in% c("tumor", "fibroblast")
meta <- meta[idx,]
cts <- cts[,idx]

segment_cols <- c("tumor" = "green4", "fibroblast" = "dodgerblue")
run_cols <- c("R1"="yellow", "R2"="orange")
slide_cols <- grDevices::colors()[grep('(gr(a|e)y)|(white)', grDevices::colors(), invert = T)] |> sample(size = 8, replace = F)
names(slide_cols) <- unique(meta$slide.name)

library(edgeR)
library(limma)

y <- DGEList(counts = cts) |> calcNormFactors()
plotMDS(x = y, col = plyr::mapvalues(x = meta$segment, from = names(segment_cols), to = segment_cols), pch = 1)
plotMDS(x = y, col = plyr::mapvalues(x = meta$Run, from = names(run_cols), to = run_cols), pch = 1)

mm <- model.matrix(~segment, data = meta)
#keep <- edgeR::filterByExpr(y = y, design = mm)
v <- edgeR::voomLmFit(counts = y,#[keep,,keep.lib.sizes = F], 
                      plot = T, 
                      design = mm, 
                      sample.weights = T,
                      block = meta$slide.name, 
                      normalize.method = "quantile")

fib_vs_tum <- (mm[meta$segment == "fibroblast", ] |> colMeans()) - (mm[meta$segment == "tumor", ] |> colMeans())
out <- contrasts.fit(fit = v, contrasts = fib_vs_tum) |> 
  eBayes() |> 
  topTable(n = Inf)
out$target <- rownames(out)

library(ggplot2)
ggplot() + 
  scattermore::geom_scattermore(data = out, mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = ifelse(out$adj.P.Val < 0.05 & abs(out$logFC) > 1, yes = "red", no = "black"), pointsize = 2) + 
  geom_vline(xintercept = c(-1, 1), color = "darkgrey", linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), color = "darkgrey", linetype = "dashed") + 
  ggrepel::geom_text_repel(data = out[abs(out$logFC) > 1.5 & out$adj.P.Val < 0.05,], mapping = aes(x = logFC, y = -log10(adj.P.Val), label = target), color = "red", min.segment.length = 0, box.padding = 0.25, size = 2, max.overlaps = 30, segment.size = 0.2) +
  ggthemes::theme_par()

# Random slopes model would be better...

lmem_de <- function(gene) {
  d <- dplyr::mutate(meta, expr=(v$EList$E[gene,] |> unlist()))
  modelout <- lmerTest::lmer(formula = expr ~ 1 + segment + (1 + segment|slide.name), data = d)
  out <- lmerTest::contest(model = modelout, L = matrix(fib_vs_tum, nrow = 1), joint = F, rhs = 0)
  out$gene <- gene
  return(out)
}

de <- parallel::mclapply(rownames(cts), lmem_de, mc.cores = 8)
dedf <- dplyr::bind_rows(de)
dedf$fdr <- p.adjust(p = dedf$`Pr(>|t|)`, method = "BH")

ggplot() + 
  scattermore::geom_scattermore(data = dedf, mapping = aes(x = Estimate, y = -log10(fdr)), 
                                color = dplyr::case_when(dedf$Estimate < -1 & dedf$fdr < 0.05 ~ "green4", 
                                                         dedf$Estimate > 1 & dedf$fdr < 0.05 ~ "dodgerblue", 
                                                         T ~ "grey"), pointsize = 2) + 
  geom_vline(xintercept = c(-1, 1), color = "black", linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed") + 
  ggrepel::geom_text_repel(data = dedf[abs(dedf$Estimate) > 1 & dedf$fdr < 0.05,], mapping = aes(x = Estimate, y = -log10(fdr), label = gene), color = "black", min.segment.length = 0, box.padding = 0.25, size = 2, max.overlaps = 15, segment.size = 0.2) +
  ggthemes::theme_par() + 
  labs(x = "log2FC", y = "-log10(adj. p-value)")
ggsave(filename = "fib_vs_tum.jpeg", width = 8, height = 8, device = "jpeg")







