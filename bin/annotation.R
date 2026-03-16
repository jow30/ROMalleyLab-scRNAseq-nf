library(dplyr)
library(tidyr)
library(purrr)
library(Seurat)
library(ggplot2)

# function for enrichment test
enrich_cluster_module <- function(cluster_genes, module_genes, background) {
  overlap <- length(intersect(cluster_genes, module_genes))
  only_cluster <- length(setdiff(cluster_genes, module_genes))
  only_module  <- length(setdiff(module_genes, cluster_genes))
  neither      <- length(setdiff(background, union(cluster_genes, module_genes)))
  
  mat <- matrix(c(overlap, only_cluster, only_module, neither), nrow = 2)
  test <- fisher.test(mat, alternative = "greater")
  data.frame(
    overlap = overlap,
    pvalue = test$p.value,
    odds_ratio = test$estimate
  )
}

marker_enrich_anno <- function(seur, marker_file){
  DefaultAssay(seur) <- "RNA"
  seur <- NormalizeData(seur)
  seur <- ScaleData(seur, features = rownames(seur))
  
  DefaultAssay(seur) <- "SCT"
  seur <- PrepSCTFindMarkers(seur)
  cluster_markers <- FindAllMarkers(seur, only.pos = T, assay = "SCT", logfc.threshold = log2(1.5), min.pct = 0.1)
  cluster_markers_sig <- cluster_markers %>% filter(p_val_adj < 0.01 & avg_log2FC > 1)

  cluster_markers_top50 <- cluster_markers_sig %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 50) %>%
    ungroup() 
  
  ref_marker <- read.csv(marker_file)
  ref_marker <- ref_marker[ref_marker$p_val_adj < 0.01 & ref_marker$avg_log2FC > 1,c("gene", "name","p_val_adj","avg_log2FC","clusterName")]
  
  ref_marker$gene <- gsub("$", ".Araport11.447", ref_marker$gene)
  
  ref_marker_top50 <- ref_marker %>%
    group_by(clusterName) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 50) %>%
    ungroup() %>% 
    select("gene","name","clusterName")
  
  module_list <- split(ref_marker_top50$gene, ref_marker_top50$clusterName)
  module_list <- lapply(module_list, function(x) x[!is.na(x)])
  # module_list <- module_list[sapply(module_list, function(g) sum(g %in% rownames(seur)))>=50]

  all_genes <- unique(rownames(seur))

  cluster_lists <- split(cluster_markers_top50$gene, cluster_markers_top50$cluster)

  # iterate over all cluster–module pairs
  enrichment_tb <- map_dfr(names(cluster_lists), function(cl) {
    cl_genes <- cluster_lists[[cl]]
    map_dfr(names(module_list), function(mod) {
      mod_genes <- module_list[[mod]]
      res <- enrich_cluster_module(cl_genes, mod_genes, all_genes)
      res$cluster <- cl
      res$module  <- mod
      res
    })
  })

  # adjust p-values (Benjamini–Hochberg)
  enrichment_tb <- enrichment_tb %>%
    mutate(p_adj = p.adjust(pvalue, method = "BH")) %>%
    arrange(p_adj)

  sig_tb <- enrichment_tb %>%
    filter(p_adj < 0.05 & overlap >= 3)

  ggplot(sig_tb, aes(x = cluster, y = module, fill = -log10(p_adj))) +
    geom_tile(color = "gray70") +
    scale_fill_viridis_c(option = "B") +
    theme_minimal() +
    labs(fill = "-log10(FDR)", x = "Cluster", y = "Module",
         title = "Module enrichment of cluster marker genes")

  # For each cluster, pick the module with the smallest FDR
  cluster_module_assign <- enrichment_tb %>%
    group_by(cluster) %>%
    slice_min(p_adj, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(cluster, module, overlap, p_adj, odds_ratio)
  
  # print(cluster_module_assign)

  # If p_adj >= 0.05, assign "Unassigned"
  cluster_module_assign$module[cluster_module_assign$p_adj >= 0.05] <- "Unknown/Contamination"

  cluster_module_assign$cluster <- as.character(cluster_module_assign$cluster)
  Idents(seur) <- as.character(Idents(seur))

  module_lookup <- setNames(cluster_module_assign$module, cluster_module_assign$cluster)
  
  # Map each cell’s cluster to its assigned module
  cell_modules <- module_lookup[as.character(Idents(seur))]
  names(cell_modules) <- colnames(seur)
  
  seur$celltype <- cell_modules
  
  return(seur)
}

anchor_transfer_anno <- function(seur, seur_ref){
  DefaultAssay(seur_ref) <- "RNA"
  seur_ref <- NormalizeData(seur_ref)
  seur_ref <- FindVariableFeatures(seur_ref)
  seur_ref <- ScaleData(seur_ref)

  DefaultAssay(seur) <- "RNA"
  seur <- NormalizeData(seur)
  seur <- FindVariableFeatures(seur)
  seur <- ScaleData(seur)

  # find anchors
  anchors <- FindTransferAnchors(reference = seur_ref, query = seur)

  # transfer labels
  predictions <- TransferData(
    anchorset = anchors,
    refdata = seur_ref$celltype
  )
  seur <- AddMetaData(object = seur, metadata = predictions)

  return(seur)
}

# # module score ------------------------------------------------------------
# 
# seur <- AddModuleScore(
#   seur,
#   module_list,
#   nbin = 30,
#   ctrl = 100,
#   assay = 'SCT',
#   name = 'module',
#   search = TRUE
# )
# # rename modules in metadata
# colnames(seur@meta.data)[grep("^module", colnames(seur@meta.data))] <- names(module_list) 
# 
# out_dir <- "featureplots_modules"
# dir.create(out_dir, showWarnings = FALSE)
# 
# for (i in seq_along(module_list)) {
#   module_name <- names(module_list)[i]
#   message("Plotting: ", module_name)
#   
#   p <- FeaturePlot(
#     seur,
#     features = names(module_list)[i],
#     cols = viridis(2),
#     keep.scale = "feature",
#     coord.fixed = TRUE,
#     order = TRUE
#   )
# 
#   ggsave(
#     file.path(out_dir, paste0("FeaturePlot_", gsub("[ /]", "_", module_name), ".png")), 
#     plot = p, width = 8, height = 6, dpi = 300
#   )
# 
#   vp <- VlnPlot(
#     seur,
#     features = names(module_list)[i],
#     alpha = 0.01
#   ) + ggtitle(module_name)
#   
#   ggsave(
#     filename = file.path(out_dir, paste0("VlnPlot_", gsub("[ /]", "_", module_name), ".png")),
#     plot = vp, width = 10, height = 8, dpi = 300
#   )
# }
# 
# 
# 
# top10 <- modules %>%
#   group_by(clusterName) %>%
#   arrange(desc(avg_log2FC), .by_group = TRUE) %>%
#   slice_head(n = 10) %>%
#   ungroup() %>% 
#   select("gene","clusterName")
# colnames(top10)<-c("Gene","geneType")
# 
# top10 <- top10[!duplicated(top10$Gene),]
# 
# gene_tmp <- rownames(seuratObj)[(apply(GetAssayData(seuratObj, layer = "counts", assay = "SCT") > 0.5, MARGIN = 1, FUN = 'sum') > 3)]
# top10 <- top10[top10$Gene %in% gene_tmp,]
# 
# top10$Gene <- factor(top10$Gene, levels = top10$Gene)
# top10$geneType <-
#   factor(
#     top10$geneType,
#     levels = unique(top10$geneType)[order(unique(top10$geneType))]
#   )
# p1 <- DotPlot(
#   seuratObj,
#   features = top10$Gene,
#   assay = "RNA",
#   cols = c('#D3D3D3', '#CC0000'),
#   col.min = 0.3,
#   col.max = 2.5,
#   scale.by = 'size',
#   dot.min = 0.01
# ) +
#   RotatedAxis() +
#   theme(
#     panel.grid.major = element_line(colour = "grey85"),
#     legend.text = element_text(color = "black", size = 12),
#     axis.text.x = element_text(
#       color = "black",
#       size = 12,
#       angle = 90,
#       vjust = 0.5
#     ),
#     axis.text.y = element_text(
#       color = "black",
#       size = 14,
#       angle = 0,
#       hjust = 1,
#       vjust = 0.5
#     ),
#     axis.title = element_blank()
#   )
# colors_50 <- c(
#   "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
#   "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
#   "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
#   "#3182bd", "#e6550d", "#31a354", "#756bb1", "#636363",
#   "#6baed6", "#fd8d3c", "#74c476", "#9e9ac8", "#bdbdbd",
#   "#9ecae1", "#fdae6b", "#a1d99b", "#bcbddc", "#d9d9d9",
#   "#c6dbef", "#fdd0a2", "#c7e9c0", "#dadaeb", "#f0f0f0",
#   "#084594", "#7f2704", "#00441b", "#3f007d", "#252525",
#   "#2171b5", "#cb181d", "#238b45", "#88419d", "#525252",
#   "#4292c6", "#ef3b2c", "#41ab5d", "#8c6bb1", "#737373",
#   "yellow","yellow2","yellow3","yellow4"
# )
# p2 <- ggplot(top10) +
#   geom_bar(
#     mapping = aes(x = Gene, y = 1, fill = geneType),
#     stat = "identity",
#     width = 1
#   ) +
#   scale_fill_manual(values = colors_50) +
#   theme_void() +
#   theme(panel.spacing.x = grid::unit(1, "mm")) +
#   theme(legend.text = element_text(color = "black", size = 12),
#         legend.title = element_blank())
# legend <-
#   plot_grid(
#     get_legend(p2),
#     get_legend(p1),
#     ncol = 2,
#     align = 'v',
#     axis = 't'
#   )
# p1 <- p1 + theme(legend.position = "none")
# p2 <- p2 + theme(legend.position = "none")
# p3 <- plot_grid(
#   p2,
#   p1,
#   align = "v",
#   ncol = 1,
#   axis = 'lr',
#   rel_heights = c(0.3, 10)
# )
# p <- plot_grid(
#   p3,
#   legend,
#   nrow = 2,
#   align = 'h',
#   axis = 'none',
#   rel_widths = c(26, 30)
# )
# ggsave(
#   paste0(resDir, "/markers.png"),
#   p,
#   width = 30,
#   height = 12,
#   units = "in"
# )  
# 
# 
# DefaultAssay(seuratObj) <- "SCT"
# markers_all <- FindAllMarkers(
#   seuratObj,
#   only.pos        = TRUE,       # only positive markers
#   min.pct         = 0.1,        # expressed in ≥10% of cells in either group
#   logfc.threshold = 0.25,       # minimum log2FC
#   test.use        = "wilcox",   # "wilcox" (default) | "MAST" | "DESeq2" | "LR" | ...
#   return.thresh   = 0.05        # adjusted p-value cutoff
# )
# top50_per_cluster <- markers_all %>%
#   group_by(cluster) %>%
#   slice_max(order_by = avg_log2FC, n = 50, with_ties = FALSE) %>%
#   ungroup()
# write.csv(top50_per_cluster, "seurat_tmp/cluster_markers_top50.csv")
# 
