#!/usr/bin/env Rscript
rm(list = ls())

suppressWarnings(suppressPackageStartupMessages(library(optparse)))

option_list <- list(
  make_option(c("-i", "--inputRds"),
              type = "character",
              default = NULL,
              help = "Input RDS files containing Seurat objects (comma-separated). If omitted, all .rds/.RDS files in current directory are used [default = auto]"),
  
  make_option(c("-o", "--resDir"),
              type = "character",
              default = ".",
              help = "Output directory [default = %default]"),
              
  make_option(c("-n", "--nFeatures"),
              type = "integer",
              default = 3000,
              help = "Number of integration features [default = %default]"),
  
  make_option(c("-m", "--integration_method"),
              type = "character",
              default = "RPCA",
              help = "Integration method to apply (RPCA, CCA, or harmony) [default = %default]",
              ),
  
  make_option(c("-s", "--species"),
              type = "character",
              default = "Arabidopsis thaliana",
              help = "Species full name [default = %default]"),
  
  make_option(c("--ref_yaml"),
              type = "character",
              default = "/project/gzy8899/qiaoshan/scRNAseq/nextflow/refs/scQC.yaml",
              help = "The organelle gene lists and annotation references of all available species curated in yaml [default = %default]"),

  make_option(c("--memory"),
              type = "integer",
              default = 16,
              help = "Memory limit in GB [default = %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# options(
#   inputRds = "/project/gzy8899/qiaoshan/scRNAseq/nextflow/test/single/results/preprocess/seur_clean_control.rds,/project/gzy8899/qiaoshan/scRNAseq/nextflow/test/single/results/preprocess/seur_clean_ABA.rds",
#   resDir = ".",
#   nFeatures = 3000,
#   integration_method = "RPCA",
#   species = "Arabidopsis thaliana",
#   ref_yaml = "/project/gzy8899/qiaoshan/scRNAseq/nextflow/refs/scQC.yaml"
# )
# opt <- options()

suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(Seurat)))
suppressWarnings(suppressPackageStartupMessages(library(glmGamPoi)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(reshape2)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))
suppressWarnings(suppressPackageStartupMessages(library(RColorBrewer)))
suppressWarnings(suppressPackageStartupMessages(library(yaml)))

required_args <- c("resDir", "species", "ref_yaml")
missing <- required_args[ sapply(required_args, function(x) is.null(opt[[x]])) ]

if (length(missing) > 0) {
  cat("###ERROR### Missing required arguments:", paste(missing, collapse = ", "), "\n\n")
  print_help(OptionParser(option_list = option_list))
  quit(status = 1)
}

qc_patterns <- read_yaml(opt$ref_yaml)
if (opt$species %in% names(qc_patterns)) {
  if(!is.null(qc_patterns[[opt$species]]$annotation_ref_seurat_obj)) {
    seurat_reference = qc_patterns[[opt$species]]$annotation_ref_seurat_obj
  }else{
    seurat_reference = NULL
  }
  if(!is.null(qc_patterns[[opt$species]]$celltype_markers)) {
    markers = qc_patterns[[opt$species]]$celltype_markers
  }else{
    markers = NULL
  }
} else {
  stop("###ERROR### This species is absent from the reference yaml. Please add annotation references to ref.yaml\n")
}

resDir <- opt$resDir
inputRds <- opt$inputRds
nFeatures <- opt$nFeatures
integration_method <- opt$integration_method
mem <- opt$memory

if (!dir.exists(resDir)) {dir.create(resDir, recursive = TRUE)}

if (is.null(inputRds) || inputRds == "") {
  inputRds_files <- list.files(resDir, pattern = "\\.[Rr][Dd][Ss]$", full.names = TRUE)
} else {
  inputRds_files <- strsplit(inputRds, ",")[[1]]
}

if (length(inputRds_files) == 0) {
  stop("###ERROR### No .rds/.RDS files found in current directory '.', and --inputRds was not provided.")
}

seuratObjs <- lapply(inputRds_files, function(r) {
  if(!file.exists(r)) {
    stop(paste("Input file not found:", r))
  }
  seur <- readRDS(r)
  seur
})
names(seuratObjs) <- sapply(seuratObjs, function(x) unique(x$orig.ident)[1])

seuratObjs.integrated <- merge(seuratObjs[[1]], y = seuratObjs[-1])
options(future.globals.maxSize = mem * 1024^3) # set memory limit to 8G

if(integration_method == "RPCA") {
  seuratObjs.integrated <- NormalizeData(seuratObjs.integrated)
  seuratObjs.integrated <- FindVariableFeatures(seuratObjs.integrated)
  seuratObjs.integrated <- ScaleData(seuratObjs.integrated)
  seuratObjs.integrated <- RunPCA(seuratObjs.integrated)
  
  seuratObjs.integrated <- IntegrateLayers(
    object = seuratObjs.integrated, method = RPCAIntegration, 
    orig.reduction = "pca", new.reduction = "integrated.rpca")
  
  umap_reduction_name <- "umap.rpca"
  cluster_name <- "rpca_clusters"
  seuratObjs.integrated <- FindNeighbors(seuratObjs.integrated, reduction = "integrated.rpca", dims = 1:30)
  seuratObjs.integrated <- FindClusters(seuratObjs.integrated, graph.name = "RNA_snn", resolution = 0.8, cluster.name = cluster_name)
  seuratObjs.integrated <- RunUMAP(seuratObjs.integrated, reduction = "integrated.rpca", dims = 1:30, reduction.name = umap_reduction_name)
}else if (integration_method == "CCA") {
  seuratObjs.integrated <- NormalizeData(seuratObjs.integrated)
  seuratObjs.integrated <- FindVariableFeatures(seuratObjs.integrated)
  seuratObjs.integrated <- ScaleData(seuratObjs.integrated)
  seuratObjs.integrated <- RunPCA(seuratObjs.integrated)
  
  seuratObjs.integrated <- IntegrateLayers(
    object = seuratObjs.integrated, method = CCAIntegration, 
    orig.reduction = "pca", new.reduction = "integrated.cca")
  
  umap_reduction_name <- "umap.cca"
  cluster_name <- "cca_clusters"
  seuratObjs.integrated <- FindNeighbors(seuratObjs.integrated, reduction = "integrated.cca", dims = 1:30)
  seuratObjs.integrated <- FindClusters(seuratObjs.integrated, graph.name = "RNA_snn", resolution = 0.8, cluster.name = cluster_name)
  seuratObjs.integrated <- RunUMAP(seuratObjs.integrated, reduction = "integrated.cca", dims = 1:30, reduction.name = umap_reduction_name)
}else if (integration_method == "harmony") {
  seuratObjs.integrated <- SCTransform(seuratObjs.integrated)
  seuratObjs.integrated <- RunPCA(seuratObjs.integrated, reduction.name = "sctpca")
  
  seuratObjs.integrated <- IntegrateLayers(
    object = seuratObjs.integrated, method = HarmonyIntegration, assay = "SCT",
    orig.reduction = "sctpca", new.reduction = "integrated.harmony")
  
  umap_reduction_name <- "umap.harmony"
  cluster_name <- "harmony_clusters"
  seuratObjs.integrated <- FindNeighbors(seuratObjs.integrated, reduction = "integrated.harmony", dims = 1:30)
  seuratObjs.integrated <- FindClusters(seuratObjs.integrated, graph.name = "SCT_snn", resolution = 0.8, cluster.name = cluster_name)
  seuratObjs.integrated <- RunUMAP(seuratObjs.integrated, reduction = "integrated.harmony", dims = 1:30, reduction.name = umap_reduction_name)
}

dp <- DimPlot(seuratObjs.integrated, reduction = umap_reduction_name, group.by = c("orig.ident"), shuffle = T)
ggsave(file.path(resDir, "UMAP_bySample.png"), plot = dp, width = 6, height = 5)
dp <- DimPlot(seuratObjs.integrated, reduction = umap_reduction_name, group.by = cluster_name, shuffle = T, label = T)
ggsave(file.path(resDir, "UMAP_byCluster.png"), plot = dp, width = 6, height = 5)
vp <- VlnPlot(seuratObjs.integrated, features = "nFeature_RNA", group.by = cluster_name, pt.size = 0.1) + NoLegend() + theme(legend.position = "none")
ggsave(file.path(resDir, "VlnPlot_nFeatureRNA_byClusters.png"), plot = vp, width = 10, height = 5)
vp <- VlnPlot(seuratObjs.integrated, features = "nCount_RNA", group.by = cluster_name, pt.size = 0.1) + NoLegend() + theme(legend.position = "none")
ggsave(file.path(resDir, "VlnPlot_nCountRNA_byClusters.png"), plot = vp, width = 10, height = 5)
# vp <- VlnPlot(seuratObjs.integrated, features = "pct.mt", group.by = cluster_name, pt.size = 0.1) + NoLegend() + theme(legend.position = "none")
# ggsave(file.path(resDir, "VlnPlot_percentMT_RNA_integration_seurat_clusters.png"), plot = vp, width = 10, height = 5)
# vp <- VlnPlot(seuratObjs.integrated, features = "pct.cp", group.by = cluster_name, pt.size = 0.1) + NoLegend() + theme(legend.position = "none")
# ggsave(file.path(resDir, "VlnPlot_percentCP_RNA_integration_seurat_clusters.png"), plot = vp, width = 10, height = 5)
# vp <- VlnPlot(seuratObjs.integrated, features = "pct.rb", group.by = cluster_name, pt.size = 0.1) + NoLegend() + theme(legend.position = "none")
# ggsave(file.path(resDir, "VlnPlot_percentRP_RNA_integration_seurat_clusters.png"), plot = vp, width = 10, height = 5)
# vp <- VlnPlot(seuratObjs.integrated, features = "score.debris", group.by = cluster_name, pt.size = 0.1) + NoLegend() + theme(legend.position = "none")
# ggsave(file.path(resDir, "VlnPlot_debris_score_RNA_integration_seurat_clusters.png"), plot = vp, width = 10, height = 5)
# vp <- VlnPlot(seuratObjs.integrated, features = "score.doublet", group.by = cluster_name, pt.size = 0.1) + NoLegend() + theme(legend.position = "none")
# ggsave(file.path(resDir, "VlnPlot_doublet_score_RNA_integration_seurat_clusters.png"), plot = vp, width = 10, height = 5)

source(Sys.which("annotation.R"))
# source("/project/gzy8899/qiaoshan/scRNAseq/nextflow/bin/annotation.R")

# Marker-enrichment annotation
if (!is.null(markers) && file.exists(markers)) {
  if (exists("marker_enrich_anno")) {
    cat("###CheckPoint### Running marker-enrichment annotation with:", markers, "\n")
    seuratObjs.integrated <- marker_enrich_anno(seuratObjs.integrated, markers)
    dp <- DimPlot(seuratObjs.integrated, group.by = "marker_anno", label = T)
    ggsave(file.path(resDir, "UMAP_byMarkerBasedCelltype.png"), plot = dp, width = 8, height = 6)
    dp <- DimPlot(seuratObjs.integrated, group.by = "marker_anno", label = T, split.by = "orig.ident", ncol = 2)
    ggsave(file.path(resDir, "UMAP_byMarkerBasedCelltype_perSample.png"), plot = dp, width = 14, height = 6)
  } else {
    cat("###WARNING### marker_enrich_anno function not found. Skipping marker annotation.\n")
  }
} else {
  cat("###CheckPoint### No marker file provided or file not found. Skipping marker annotation.\n")
}

# Anchor-transfer annotation
if (!is.null(seurat_reference) && file.exists(seurat_reference)) {
  if (exists("anchor_transfer_anno")) {
    cat("###CheckPoint### Running anchor-transfer annotation with:", seurat_reference, "\n")
    seur_ref <- readRDS(seurat_reference)
    seuratObjs.integrated <- JoinLayers(seuratObjs.integrated, assay = "RNA")
    if (!("anchor_anno" %in% colnames(seuratObjs.integrated@meta.data))) {
      seuratObjs.integrated$anchor_anno <- seuratObjs.integrated$predicted.id
    }
    seuratObjs.integrated$anchor_anno_solo <- seuratObjs.integrated$anchor_anno
    seuratObjs.integrated <- anchor_transfer_anno(seuratObjs.integrated, seur_ref)
    
    dp <- DimPlot(seuratObjs.integrated, group.by = "anchor_anno", label = T, label.size = 3) + theme(legend.position = "right")
    ggsave(file.path(resDir, "UMAP_byCelltype.png"), plot = dp, width = 12, height = 6)
    
    dp <- DimPlot(seuratObjs.integrated, group.by = "anchor_anno", split.by = "orig.ident", label = T, label.size = 2, ncol = 3) + theme(legend.position = "bottom")
    h <- 8+3*length(unique(seuratObjs.integrated$orig.ident))/3
    ggsave(file.path(resDir, "UMAP_byCelltype_perSample.png"), plot = dp, width = 14, height = h)
    
    dp <- DimPlot(seuratObjs.integrated, split.by = "anchor_anno", group.by = "anchor_anno", label = F, label.size = 3, ncol = 5) + theme(legend.position = "none")
    h <- 3+3*length(unique(seuratObjs.integrated$anchor_anno))/6
    ggsave(file.path(resDir, "UMAP_perCelltype.png"), plot = dp, width = 20, height = h)
    
    dp <- DimPlot(seuratObjs.integrated, split.by = "anchor_anno", group.by = "orig.ident", label = F, label.size = 3, ncol = 5) + theme(legend.position = "none")
    h <- 3+3*length(unique(seuratObjs.integrated$anchor_anno))/6
    ggsave(file.path(resDir, "UMAP_bySample_perCelltype.png"), plot = dp, width = 20, height = h)
    
  } else {
    cat("###WARNING### anchor_transfer_anno function not found. Skipping anchor-transfer annotation.\n")
  }
} else {
  cat("###CheckPoint### No Seurat reference provided or file not found. Skipping anchor-transfer annotation.\n")
}

saveRDS(seuratObjs.integrated, paste0(resDir, "/seuratObjs_", integration_method, "_integrated.rds"))

# percentage bar plots of samples in each cluster --------------------------------

if ("anchor_anno" %in% colnames(seuratObjs.integrated@meta.data)) {
  init_df <- data.frame(cluster = seuratObjs.integrated$anchor_anno, sample = seuratObjs.integrated$orig.ident)
  trans_df <- dcast(init_df, cluster ~ sample, fun.aggregate = length)
  long_df <- pivot_longer(trans_df, cols = -cluster, names_to = "sample", values_to = "Count")
  long_df <- long_df %>% group_by(cluster) %>% mutate(Percentage = (Count / sum(Count)) * 100)
  wide_df <- long_df %>% select(-Count) %>% pivot_wider(names_from = sample, values_from = Percentage)

  write.csv(trans_df, file = paste0(resDir, "/cell_counts_", integration_method, "_samples_perCelltype.csv"))
  write.csv(wide_df, file = paste0(resDir, "/cell_percentage_", integration_method, "_samples_perCelltype.csv"))

  long_df$sample <- factor(long_df$sample, levels = unique(long_df$sample))
  sample_color <- c(brewer.pal(n = 8, name = 'Set2'), brewer.pal(n = 9, name = 'Set1'), brewer.pal(n = 12, name = 'Set3'), brewer.pal(n = 12, name = 'Paired'))[1:length(levels(long_df$sample))]
  names(sample_color) <- levels(long_df$sample)

  p <- ggplot(long_df, aes(x = cluster, y = Percentage, fill = sample)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = NULL, y = NULL, fill = "Sample") +
    scale_fill_manual(values = sample_color) +
    theme(text = element_text(size = 20)) +
    theme(axis.text = element_text(colour = "black")) +
    coord_flip()
  ggsave(
    paste0(resDir, "/percentage_barplot_", integration_method, "_samples_perCelltype.png"),
    p,
    width = 8,
    height = 8,
    units = "in"
  )
} else {
  cat("###CheckPoint### 'anchor_anno' column not found. Skipping celltype barplot.\n")
}

# percentage bar plots of clusters for each sample --------------------------------

if ("anchor_anno" %in% colnames(seuratObjs.integrated@meta.data)) {
  init_df <- data.frame(cluster = seuratObjs.integrated$anchor_anno, sample = seuratObjs.integrated$orig.ident)
  trans_df <- dcast(init_df, sample ~ cluster, fun.aggregate = length)
  long_df <- pivot_longer(trans_df, cols = -sample, names_to = "cluster", values_to = "Count")
  long_df <- long_df %>% group_by(sample) %>% mutate(Percentage = (Count / sum(Count)) * 100)
  wide_df <- long_df %>% select(-Count) %>% pivot_wider(names_from = cluster, values_from = Percentage)

  write.csv(trans_df, file = paste0(resDir, "/cell_counts_", integration_method, "_celltype_perSample.csv"))
  write.csv(wide_df, file = paste0(resDir, "/cell_percentage_", integration_method, "_celltype_perSample.csv"))

  long_df$sample <- factor(long_df$sample, levels = rev(unique(long_df$sample)))
  long_df$cluster <- factor(long_df$cluster, levels = unique(long_df$cluster))
  cluster_color <- c(brewer.pal(n = 8, name = 'Set2'), brewer.pal(n = 9, name = 'Set1'), brewer.pal(n = 12, name = 'Set3'), brewer.pal(n = 12, name = 'Paired'))[1:length(levels(long_df$cluster))]
  names(cluster_color) <- levels(long_df$cluster)

  p <- ggplot(long_df, aes(x = sample, y = Percentage, fill = cluster)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = NULL, y = NULL, fill = "Celltype") +
    scale_fill_manual(values = cluster_color) +
    theme(text = element_text(size = 20)) +
    theme(axis.text = element_text(colour = "black")) +
    coord_flip()
  ggsave(
    paste0(resDir, "/percentage_barplot_", integration_method, "_celltype_per_sample.png"),
    p,
    width = 20,
    height = 6,
    units = "in"
  )
} else {
  cat("###CheckPoint### 'anchor_anno' column not found. Skipping anchor_anno barplot.\n")
}


