rm(list = ls())

library(dplyr)
library(Seurat)
library(glmGamPoi)
library(patchwork)
library(ggplot2)
library(reshape2)
library(tidyr)
library(RColorBrewer)

suppressWarnings(suppressPackageStartupMessages(library(optparse)))

option_list <- list(
  make_option(c("-i", "--inputRds"),
              type = "character",
              help = "Input RDS files containing Seurat objects (comma-separated) [REQUIRED]"),
  
  make_option(c("-o", "--resDir"),
              type = "character",
              help = "Output directory [REQUIRED]"),
              
  make_option(c("-n", "--nFeatures"),
              type = "integer",
              default = 3000,
              help = "Number of integration features [default = %default]"),
  
  make_option(c("--markers"),
              type = "character",
              default = NULL,
              help = "Path to marker CSV file for enrichment-based annotation [default = NULL]"),
  
  make_option(c("--seurat_reference"),
              type = "character",
              default = NULL,
              help = "Path to Seurat reference RDS file for anchor-transfer annotation [default = NULL]")
)

opt <- parse_args(OptionParser(option_list=option_list))

required_args <- c("inputRds", "resDir")
missing <- required_args[ sapply(required_args, function(x) is.null(opt[[x]])) ]

if (length(missing) > 0) {
  cat("###ERROR### Missing required arguments:", paste(missing, collapse = ", "), "\n\n")
  print_help(OptionParser(option_list = option_list))
  quit(status = 1)
}

resDir <- opt$resDir
inputRds <- opt$inputRds
nFeatures <- opt$nFeatures

if (!dir.exists(resDir)) {dir.create(resDir, recursive = TRUE)}

inputRds_files <- strsplit(inputRds, ",")[[1]]
seuratObjs <- lapply(inputRds_files, function(r) {
  if(!file.exists(r)) {
    stop(paste("Input file not found:", r))
  }
  seur <- readRDS(r)
  DefaultAssay(seur) <- "RNA"
  seur <- SCTransform(seur, verbose = FALSE)
  seur
})
names(seuratObjs) <- sapply(seuratObjs, function(x) unique(x$orig.ident)[1])

features <- SelectIntegrationFeatures(seuratObjs, nfeatures = nFeatures)
seuratObjs <- PrepSCTIntegration(seuratObjs, anchor.features = features)
anchors <- FindIntegrationAnchors(seuratObjs, normalization.method = "SCT", anchor.features = features, reduction = "rpca")

# Seurat SCT integration can fail when the number of anchor cells is
# smaller than k.weight (default: 100). In such cases, retry with a
# lowered k.weight parsed from Seurat's error message.
k_weight_initial <- 100
seuratObjs.integrated <- tryCatch(
  {
    IntegrateData(anchors, normalization.method = "SCT", k.weight = k_weight_initial)
  },
  error = function(e) {
    msg <- conditionMessage(e)
    # Example: "Number of anchor cells is less than k.weight. Consider lowering k.weight to less than 74 ..."
    if (grepl("Number of anchor cells is less than k\\.weight", msg) || grepl("Number of anchor cells is less than k.weight", msg)) {
      m <- regmatches(msg, regexec("less than ([0-9]+)", msg))
      if (length(m) >= 2) {
        suggested_upper <- as.numeric(m[2])
        retry_k_weight <- max(suggested_upper - 1, 1)
      } else {
        # Fallback when we cannot parse the suggested bound.
        retry_k_weight <- 50
      }
      cat("###WARNING### Seurat::IntegrateData failed with default k.weight =", k_weight_initial, "\n")
      cat("###WARNING### Retrying with k.weight =", retry_k_weight, " (parsed from error)\n")
      return(
        IntegrateData(
          anchors,
          normalization.method = "SCT",
          k.weight = retry_k_weight
        )
      )
    }
    stop(e)
  }
)

seuratObjs.integrated <- RunPCA(seuratObjs.integrated)
seuratObjs.integrated <- RunUMAP(seuratObjs.integrated, dims = 1:30)
seuratObjs.integrated <- FindNeighbors(seuratObjs.integrated, dims = 1:30)
seuratObjs.integrated <- FindClusters(seuratObjs.integrated, resolution = 0.8)

dp <- DimPlot(seuratObjs.integrated, reduction = "umap", group.by = c("orig.ident"), shuffle = T)
ggsave(file.path(resDir, "UMAP_integration.png"), plot = dp, width = 6, height = 5)
dp <- DimPlot(seuratObjs.integrated, reduction = "umap", group.by = c("seurat_clusters"), shuffle = T, label = T)
ggsave(file.path(resDir, "UMAP_integration_seurat_clusters.png"), plot = dp, width = 6, height = 5)
vp <- VlnPlot(seuratObjs.integrated, features = "nFeature_RNA", group.by = "seurat_clusters", pt.size = 0.1) + NoLegend() + theme(legend.position = "none")
ggsave(file.path(resDir, "VlnPlot_nFeature_RNA_integration_seurat_clusters.png"), plot = vp, width = 10, height = 5)
vp <- VlnPlot(seuratObjs.integrated, features = "nCount_RNA", group.by = "seurat_clusters", pt.size = 0.1) + NoLegend() + theme(legend.position = "none")
ggsave(file.path(resDir, "VlnPlot_nCount_RNA_integration_seurat_clusters.png"), plot = vp, width = 10, height = 5)
vp <- VlnPlot(seuratObjs.integrated, features = "pct.mt", group.by = "seurat_clusters", pt.size = 0.1) + NoLegend() + theme(legend.position = "none")
ggsave(file.path(resDir, "VlnPlot_percentMT_RNA_integration_seurat_clusters.png"), plot = vp, width = 10, height = 5)
vp <- VlnPlot(seuratObjs.integrated, features = "pct.cp", group.by = "seurat_clusters", pt.size = 0.1) + NoLegend() + theme(legend.position = "none")
ggsave(file.path(resDir, "VlnPlot_percentCP_RNA_integration_seurat_clusters.png"), plot = vp, width = 10, height = 5)
vp <- VlnPlot(seuratObjs.integrated, features = "pct.rb", group.by = "seurat_clusters", pt.size = 0.1) + NoLegend() + theme(legend.position = "none")
ggsave(file.path(resDir, "VlnPlot_percentRP_RNA_integration_seurat_clusters.png"), plot = vp, width = 10, height = 5)
vp <- VlnPlot(seuratObjs.integrated, features = "score.debris", group.by = "seurat_clusters", pt.size = 0.1) + NoLegend() + theme(legend.position = "none")
ggsave(file.path(resDir, "VlnPlot_debris_score_RNA_integration_seurat_clusters.png"), plot = vp, width = 10, height = 5)
vp <- VlnPlot(seuratObjs.integrated, features = "score.doublet", group.by = "seurat_clusters", pt.size = 0.1) + NoLegend() + theme(legend.position = "none")
ggsave(file.path(resDir, "VlnPlot_doublet_score_RNA_integration_seurat_clusters.png"), plot = vp, width = 10, height = 5)

saveRDS(seuratObjs.integrated, file.path(resDir, "seuratObjs_integrated.rds"))

# ── Cell annotation (optional, parameterized) ──────────────────────────────

source(Sys.which("annotation.R"))

# Marker-enrichment annotation
if (!is.null(opt$markers) && file.exists(opt$markers)) {
  if (exists("marker_enrich_anno")) {
    cat("###CheckPoint### Running marker-enrichment annotation with:", opt$markers, "\n")
    seuratObjs.integrated <- marker_enrich_anno(seuratObjs.integrated, opt$markers)
    dp <- DimPlot(seuratObjs.integrated, group.by = "celltype", label = T)
    ggsave(file.path(resDir, "UMAP_integration_celltype.png"), plot = dp, width = 8, height = 6)
    dp <- DimPlot(seuratObjs.integrated, group.by = "celltype", label = T, split.by = "orig.ident", ncol = 2)
    ggsave(file.path(resDir, "UMAP_integration_celltype_by_sample.png"), plot = dp, width = 14, height = 6)
  } else {
    cat("###WARNING### marker_enrich_anno function not found. Skipping marker annotation.\n")
  }
} else {
  cat("###CheckPoint### No marker file provided or file not found. Skipping marker annotation.\n")
}

# Anchor-transfer annotation
if (!is.null(opt$seurat_reference) && file.exists(opt$seurat_reference)) {
  if (exists("anchor_transfer_anno")) {
    cat("###CheckPoint### Running anchor-transfer annotation with:", opt$seurat_reference, "\n")
    seur_ref <- readRDS(opt$seurat_reference)
    seuratObjs.integrated <- UpdateSeuratObject(seuratObjs.integrated)
    rownames(seuratObjs.integrated) <- gsub("\\.Araport.*$", "", rownames(seuratObjs.integrated))
    rownames(seuratObjs.integrated) <- gsub("\\.v.*$", "", rownames(seuratObjs.integrated))
    seuratObjs.integrated <- anchor_transfer_anno(seuratObjs.integrated, seur_ref)
    dp <- DimPlot(seuratObjs.integrated, group.by = "predicted.id", label = T, label.size = 3)
    ggsave(file.path(resDir, "UMAP_integration_predicted_id.png"), plot = dp, width = 8, height = 6)
    dp <- DimPlot(seuratObjs.integrated, group.by = "predicted.id", split.by = "orig.ident", label = T, ncol = 2, label.size = 2)
    ggsave(file.path(resDir, "UMAP_integration_predicted_id_by_sample.png"), plot = dp, width = 14, height = 6)
    dp <- DimPlot(seuratObjs.integrated, split.by = "predicted.id", group.by = "predicted.id", label = F, label.size = 3, ncol = 6)
    ggsave(file.path(resDir, "UMAP_integration_predicted_id_split.png"), plot = dp, width = 20, height = 10)
  } else {
    cat("###WARNING### anchor_transfer_anno function not found. Skipping anchor-transfer annotation.\n")
  }
} else {
  cat("###CheckPoint### No Seurat reference provided or file not found. Skipping anchor-transfer annotation.\n")
}

saveRDS(seuratObjs.integrated, file.path(resDir, "seuratObjs_integrated.rds"))

# seuratObjs_merged <- merge(x = seuratObjs[[1]], y = c(seuratObjs[[2]], seuratObjs[[3]]))

# options(future.globals.maxSize = 3e+10)
# seuratObjs_merged <- IntegrateLayers(object = seuratObjs_merged, method = CCAIntegration, new.reduction = "integrated.cca", normalization.method = "SCT", verbose = F)
# seuratObjs_merged <- IntegrateLayers(object = seuratObjs_merged, method = RPCAIntegration, new.reduction = "integrated.rpca", normalization.method = "SCT", verbose = F)
# seuratObjs_merged <- IntegrateLayers(object = seuratObjs_merged, method = HarmonyIntegration, new.reduction = "harmony", normalization.method = "SCT", verbose = F)
# 
# seuratObjs_merged <- FindNeighbors(seuratObjs_merged, dims = 1:30, reduction = "integrated.cca")
# seuratObjs_merged <- FindClusters(seuratObjs_merged, resolution = 0.8, cluster.name = "cca_clusters")
# seuratObjs_merged <- RunUMAP(seuratObjs_merged, dims = 1:30, reduction = "integrated.cca", reduction.name = "umap.cca")
# 
# seuratObjs_merged <- FindNeighbors(seuratObjs_merged, dims = 1:30, reduction = "integrated.rpca")
# seuratObjs_merged <- FindClusters(seuratObjs_merged, resolution = 0.8, cluster.name = "rpca_clusters")
# seuratObjs_merged <- RunUMAP(seuratObjs_merged, dims = 1:30, reduction = "integrated.rpca", reduction.name = "umap.rpca")
# 
# seuratObjs_merged <- FindNeighbors(seuratObjs_merged, dims = 1:30, reduction = "harmony")
# seuratObjs_merged <- FindClusters(seuratObjs_merged, resolution = 0.8, cluster.name = "harmony_clusters") # latest clustering results will be stored in object metadata under 'seurat_clusters'. Idents will be updated, too. Note that 'seurat_clusters' will be overwritten everytime FindClusters is run
# seuratObjs_merged <- RunUMAP(seuratObjs_merged, dims = 1:30, reduction = "harmony", reduction.name = "umap.harmony")
# 
# dp <- DimPlot(seuratObjs_merged, reduction = "umap.cca", group.by = c("orig.ident"))
# ggsave(file.path(resDir, "UMAP_after_integration_all_CCA.png"), plot = dp, width = 6, height = 5)
# dp <- DimPlot(seuratObjs_merged, reduction = "umap.cca", group.by = c("cca_clusters"))
# ggsave(file.path(resDir, "UMAP_after_integration_all_CCA_clusters.png"), plot = dp, width = 6, height = 5)
# 
# dp <- DimPlot(seuratObjs_merged, reduction = "umap.rpca", group.by = c("orig.ident"))
# ggsave(file.path(resDir, "UMAP_after_integration_all_RPCA.png"), plot = dp, width = 6, height = 5)
# dp <- DimPlot(seuratObjs_merged, reduction = "umap.rpca", group.by = c("rpca_clusters"))
# ggsave(file.path(resDir, "UMAP_after_integration_all_RPCA_clusters.png"), plot = dp, width = 6, height = 5)
# 
# dp <- DimPlot(seuratObjs_merged, reduction = "umap.harmony", group.by = c("orig.ident"))
# ggsave(file.path(resDir, "UMAP_after_integration_all_Harmony.png"), plot = dp, width = 6, height = 5)
# dp <- DimPlot(seuratObjs_merged, reduction = "umap.harmony", group.by = c("harmony_clusters"))
# ggsave(file.path(resDir, "UMAP_after_integration_all_Harmony_clusters.png"), plot = dp, width = 6, height = 5)
# 
# vp <- VlnPlot(seuratObjs_merged, features = "nFeature_RNA", group.by = "cca_clusters", pt.size = 0.1) + NoLegend() + theme(legend.position = "none")
# ggsave(file.path(resDir, "VlnPlot_nFeature_RNA_after_integration_all_CCA_clusters.png"), plot = vp, width = 10, height = 5)
# vp <- VlnPlot(seuratObjs_merged, features = "nFeature_RNA", group.by = "rpca_clusters", pt.size = 0.1) + NoLegend() + theme(legend.position = "none")
# ggsave(file.path(resDir, "VlnPlot_nFeature_RNA_after_integration_all_RPCA_clusters.png"), plot = vp, width = 10, height = 5)
# vp <- VlnPlot(seuratObjs_merged, features = "nFeature_RNA", group.by = "harmony_clusters", pt.size = 0.1) + NoLegend() + theme(legend.position = "none")
# ggsave(file.path(resDir, "VlnPlot_nFeature_RNA_after_integration_all_Harmony_clusters.png"), plot = vp, width = 10, height = 5)



# percentage bar plots of samples in each cluster --------------------------------

if ("celltype" %in% colnames(seuratObjs.integrated@meta.data)) {
  init_df <- data.frame(cluster = seuratObjs.integrated$celltype, sample = seuratObjs.integrated$orig.ident)
  trans_df <- dcast(init_df, cluster ~ sample, fun.aggregate = length)
  long_df <- pivot_longer(trans_df, cols = -cluster, names_to = "sample", values_to = "Count")
  long_df <- long_df %>% group_by(cluster) %>% mutate(Percentage = (Count / sum(Count)) * 100)
  wide_df <- long_df %>% select(-Count) %>% pivot_wider(names_from = sample, values_from = Percentage)

  write.csv(trans_df, file = paste0(resDir, "/cell_counts_samples_per_harmony_cluster.csv"))
  write.csv(wide_df, file = paste0(resDir, "/cell_percentage_samples_per_harmony_cluster.csv"))

  long_df$sample <- factor(long_df$sample, levels = unique(long_df$sample))
  sample_color <- brewer.pal(n = 7, name = 'Set2')
  names(sample_color) <- levels(long_df$sample)

  p <- ggplot(long_df, aes(x = cluster, y = Percentage, fill = sample)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = NULL, y = NULL, fill = "Sample") +
    scale_fill_manual(values = sample_color) +
    theme(text = element_text(size = 20)) +
    theme(axis.text = element_text(colour = "black")) +
    coord_flip()
  ggsave(
    file.path(resDir, "percentage_barplot_samples_per_harmony_cluster.png"),
    p,
    width = 8,
    height = 8,
    units = "in"
  )
} else {
  cat("###CheckPoint### 'celltype' column not found. Skipping celltype barplot.\n")
}

# percentage bar plots of clusters for each sample --------------------------------

if ("predicted.id" %in% colnames(seuratObjs.integrated@meta.data)) {
  init_df <- data.frame(cluster = seuratObjs.integrated$predicted.id, sample = seuratObjs.integrated$orig.ident)
  trans_df <- dcast(init_df, sample ~ cluster, fun.aggregate = length)
  long_df <- pivot_longer(trans_df, cols = -sample, names_to = "cluster", values_to = "Count")
  long_df <- long_df %>% group_by(sample) %>% mutate(Percentage = (Count / sum(Count)) * 100)
  wide_df <- long_df %>% select(-Count) %>% pivot_wider(names_from = cluster, values_from = Percentage)

  write.csv(trans_df, file = paste0(resDir, "/cell_counts_harmony_clusters_per_sample.csv"))
  write.csv(wide_df, file = paste0(resDir, "/cell_percentage_harmony_clusters_per_sample.csv"))

  long_df$sample <- factor(long_df$sample, levels = rev(unique(long_df$sample)))
  long_df$cluster <- factor(long_df$cluster, levels = unique(long_df$cluster))
  cluster_color <- c(brewer.pal(n = 13, name = 'Set2'), brewer.pal(n = 13, name = 'Set1'), brewer.pal(n = 13, name = 'Set3'), brewer.pal(n = 13, name = 'Paired'))[1:length(levels(long_df$cluster))]
  names(cluster_color) <- levels(long_df$cluster)

  p <- ggplot(long_df, aes(x = sample, y = Percentage, fill = cluster)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = NULL, y = NULL, fill = "Celltype") +
    scale_fill_manual(values = cluster_color) +
    theme(text = element_text(size = 20)) +
    theme(axis.text = element_text(colour = "black")) +
    coord_flip()
  ggsave(
    file.path(resDir, "percentage_barplot_harmony_clusters_per_sample.png"),
    p,
    width = 10,
    height = 6,
    units = "in"
  )
} else {
  cat("###CheckPoint### 'predicted.id' column not found. Skipping predicted.id barplot.\n")
}


