rm(list = ls())

suppressWarnings(suppressPackageStartupMessages(library(optparse)))

option_list <- list(
  make_option(c("-s", "--sample"),
              type = "character",
              help = "Sample name [REQUIRED]"),

  make_option(c("-S", "--species"),
              type = "character",
              help = "Species full name [REQUIRED]"),

  make_option(c("-o", "--outDir"),
              type = "character",
              help = "Output directory [REQUIRED]"),

  make_option(c("-v", "--velocyto_txt"),
              type = "character",
              help = "Velocyto txt output file path [REQUIRED]"),

  # Optional arguments with defaults

  make_option(c("--ref_yaml"),
              type = "character",
              default = "/project/gzy8899/qiaoshan/scRNAseq/references/scQC.yaml",
              help = "The reference organelle gene lists of all available species curated in yaml [default = %default]"),

  make_option(c("--min_unsplice_ratio"),
              type = "numeric",
              default = 0.1,
              help = "Minimum Velocyto-derived unsplice ratio [default = %default]"),

  make_option(c("--min_nCount_RNA"),
              type = "integer",
              default = 400,
              help = "Minimum UMI per cell [default = %default]"),

  make_option(c("--min_nFeature_RNA"),
              type = "integer",
              default = 300,
              help = "Minimum genes per cell [default = %default]"),

  make_option(c("--max_mt"),
              type = "double",
              default = 10,
              help = "Maximum percent mitochondrial gene expression [default = %default]"),

  make_option(c("--max_cp"),
              type = "double",
              default = 15,
              help = "Maximum percent chloroplast gene expression [default = %default]"),

  make_option(c("--nHVG"),
              type = "integer",
              default = 3000,
              help = "Number of highly variable genes (HVGs) [default = %default]"),

  make_option(c("--min_ncell_expr"),
              type = "integer",
              default = 5,
              help = "Minimum number of cells expressing a gene [default = %default]"),

  make_option(c("--remove_doublet"),
              type = "logical",
              default = TRUE,
              help = "Whether to remove doublets [default = %default]"),

  make_option(c("--max_doublet_score"),
              type = "double",
              default = 0.4,
              help = "Maximum doubletFinder pANN [default = %default]"),

  make_option(c("--min_nClusterMarker"),
              type = "integer",
              default = 5,
              help = "Minimum number of markers detected for a valid cluster [default = %default]"),

  make_option(c("--min_cells"),
              type = "integer",
              default = 500,
              help = "Minimum number of cells required to continue; samples with fewer cells are skipped [default = %default]"),

  make_option(c("--threads"),
              type = "integer",
              default = 1,
              help = "Threads [default = %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# init --------------------------------------------------------------------

required_args <- c("sample", "species", "outDir", "velocyto_txt")

missing <- required_args[ sapply(required_args, function(x) is.null(opt[[x]])) ]

if (length(missing) > 0) {
  cat("###ERROR### Missing required arguments:", paste(missing, collapse = ", "), "\n\n")
  print_help(OptionParser(option_list = option_list))
  quit(status = 1)
}

suppressWarnings(suppressPackageStartupMessages(library(yaml)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(Seurat)))
suppressWarnings(suppressPackageStartupMessages(library(glmGamPoi)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))
suppressWarnings(suppressPackageStartupMessages(library(SingleCellExperiment)))
suppressWarnings(suppressPackageStartupMessages(library(cowplot)))
suppressWarnings(suppressPackageStartupMessages(library(viridis)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(Matrix)))
suppressWarnings(suppressPackageStartupMessages(library(DoubletFinder)))
suppressWarnings(suppressPackageStartupMessages(library(ggpubr)))
suppressWarnings(suppressPackageStartupMessages(library(scDblFinder)))

cat("###CheckPoint### Arguments\n")
cat("###CheckPoint### ---------\n")
cat("###CheckPoint### Sample     :", opt$sample, "\n")
cat("###CheckPoint### Species    :", opt$species, "\n")
cat("###CheckPoint### Output dir :", opt$outDir, "\n")
cat("###CheckPoint### Velocyto txt", opt$velocyto_txt, "\n")
cat("###CheckPoint### ref_yaml         :", opt$ref_yaml, "\n")
cat("###CheckPoint### min_velocyto_unsplice_ratio  :", opt$min_unsplice_ratio, "\n")
cat("###CheckPoint### min_nCount_RNA   :", opt$min_nCount_RNA, "\n")
cat("###CheckPoint### min_nFeature_RNA :", opt$min_nFeature_RNA, "\n")
cat("###CheckPoint### max_mt :", opt$max_mt, "\n")
cat("###CheckPoint### max_cp :", opt$max_cp, "\n")
cat("###CheckPoint### nHVG   :", opt$nHVG, "\n")
cat("###CheckPoint### min_ncell_expr   :", opt$min_ncell_expr, "\n")
cat("###CheckPoint### remove_doublet   :", opt$remove_doublet, "\n")
cat("###CheckPoint### max_doublet_score:", opt$max_doublet_score, "\n")
cat("###CheckPoint### min_nClusterMarker:", opt$min_nClusterMarker, "\n")
cat("###CheckPoint### threads:", opt$threads, "\n\n")

qc_patterns <- read_yaml(opt$ref_yaml)
if (opt$species %in% names(qc_patterns)) {
  if(!is.null(qc_patterns[[opt$species]]$mitochondrial_pattern)) {
    mitochondrial_pattern = qc_patterns[[opt$species]]$mitochondrial_pattern
  }else{
    mitochondrial_pattern = NULL
    mt_genes = readLines(qc_patterns[[opt$species]]$mitochondrial_genes)
  }
  if(!is.null(qc_patterns[[opt$species]]$chloroplast_pattern)) {
    chloroplast_pattern = qc_patterns[[opt$species]]$chloroplast_pattern
  }else{
    chloroplast_pattern = NULL
    cp_genes = readLines(qc_patterns[[opt$species]]$chloroplast_genes)
  }
  ribosomal_genes = readLines(qc_patterns[[opt$species]]$ribosomal_genes)
} else {
  stop("###ERROR### This species is absent from the reference yaml. Please add mitochondrial/chloroplast/ribosomal genes to ref.yaml\n")
}

if (!dir.exists(opt$outDir)) {dir.create(opt$outDir, recursive = T)}

setwd(opt$outDir)

saveRDS(opt, paste0("summary_opts_", opt$sample, ".rds"))

set.seed(123)

# Exit gracefully when too few cells remain; Nextflow treats exit 0 as success
# so the process produces no outputs and downstream channels drop this sample.
check_cells <- function(obj, step, min_cells = 500) {
  n <- ncol(obj)
  if (n < min_cells) {
    cat("###CheckPoint### ERROR: Only", n, "cells remaining after", step,
        "(minimum required:", min_cells, "). Skipping this sample.\n")
    quit(save = "no", status = 0)
  }
}

# Load RDS files produced by preprocess_initial.R (staged into work dir by Nextflow)
summary_dims <- readRDS(paste0("summary_dims_", opt$sample, ".rds"))
summary_tbs  <- readRDS(paste0("summary_tbs_",  opt$sample, ".rds"))
summary_plts <- readRDS(paste0("summary_plts_", opt$sample, ".rds"))
seur_objs    <- readRDS(paste0("seur_objs_",    opt$sample, ".rds"))
seur_diem    <- seur_objs[["After DIEM Debris Removal"]]

# Un-spliced ratio filtering -----------------------------------------------

velocyto_tb <- read.table(opt$velocyto_txt, header = T)

rownames(velocyto_tb) <- gsub(".*:(.*)x", "\\1", velocyto_tb$cell_barcode)
velocyto_tb <- velocyto_tb[,-1]

velocyto_tb <- velocyto_tb[colnames(seur_diem), ]

seur_diem <- AddMetaData(seur_diem, metadata = velocyto_tb)

vp <- VlnPlot(seur_diem, assay = "RNA", layer = "counts", features = c("total_counts", "n_genes", "pct.mt", "pct.cp", "pct.rb", "score.debris", "unsplice_ratio"), pt.size = 0.01, alpha = 0.1, ncol = 7)
# ggsave(paste0("QC_feature_violin_diem_", opt$sample, ".png"), plot = vp, width = 18, height = 7, units = "in")
summary_plts[["After DIEM Debris Removal"]] <- c("Violin"=vp)
summary_tbs[["After DIEM Debris Removal"]] <- do.call(rbind, lapply(seur_diem@meta.data[c("total_counts","n_genes","pct.mt","pct.cp","pct.rb","score.debris", "unsplice_ratio")], summary))

seur_diem_velocyto <- subset(seur_diem, subset = unsplice_ratio > opt$min_unsplice_ratio)

vp <- VlnPlot(seur_diem_velocyto, assay = "RNA", layer = "counts", features = c("total_counts", "n_genes", "pct.mt", "pct.cp", "pct.rb", "score.debris", "unsplice_ratio"), pt.size = 0.01, alpha = 0.1, ncol = 7)
# ggsave(paste0("QC_feature_violin_diem_velocyto_", opt$sample, ".png"), plot = vp, width = 18, height = 7, units = "in")
summary_plts[["After Unspliced Ratio Filtering"]] <- c("Violin"=vp)
summary_dims[["After Unspliced Ratio Filtering"]] <- c("nCell"=ncol(seur_diem_velocyto), "nGene"=nrow(seur_diem_velocyto))
summary_tbs[["After Unspliced Ratio Filtering"]] <- do.call(rbind, lapply(seur_diem_velocyto@meta.data[c("total_counts","n_genes","pct.mt","pct.cp","pct.rb","score.debris", "unsplice_ratio")], summary))
seur_objs[["After Unspliced Ratio Filtering"]] <- seur_diem_velocyto

cat("###CheckPoint### Number of remaining cell barcodes after removing low un-spliced ratio:", ncol(seur_diem_velocyto), "\n")
check_cells(seur_diem_velocyto, "unspliced ratio filtering", opt$min_cells)

# Remove ribosomal genes --------------------------------------------------
# Here the ribosomal genes in the list are actually ribosomal-protein-coding genes. They can be quantified as usual.
# They don't need to be removed unless they are highly expressed and dominate the cell clustering.
# seur_diem_velocyto <- seur_diem_velocyto[ !(rownames(seur_diem_velocyto) %in% ribosomal_genes), ]
# cat("###CheckPoint###", length(ribosomal_genes), "ribosomal genes removed.\n")

# Count/Feature/CP/MT filtering ---------------------------------------------------------

# automatic max cutoff setting
max_nCount_RNA  <- as.integer(quantile(seur_diem_velocyto$nCount_RNA, probs = 0.999, na.rm = TRUE))
cat("###CheckPoint### The maximum UMI count is", max_nCount_RNA, "\n")

seur_diem_velocyto_flts <- subset(seur_diem_velocyto, subset = nCount_RNA >= opt$min_nCount_RNA & nCount_RNA <= max_nCount_RNA & nFeature_RNA >= opt$min_nFeature_RNA & pct.mt < opt$max_mt & pct.cp < opt$max_cp)

# Filter out features/genes detected in less than 5 cells
mat <- GetAssayData(seur_diem_velocyto_flts, layer = "counts")
keep <- Matrix::rowSums(mat > 0) > opt$min_ncell_expr
seur_diem_velocyto_flts <- seur_diem_velocyto_flts[keep, ]

cat("###CheckPoint### Number of cell barcodes left after Count/Feature/percentCP/percentMT filtering:", ncol(seur_diem_velocyto_flts), "\n")
cat("###CheckPoint### Number of genes left after Count/Feature/percentCP/percentMT filtering:", nrow(seur_diem_velocyto_flts), "\n\n")

if(!is.null(mitochondrial_pattern)) {mt_genes <- grep(mitochondrial_pattern, rownames(seur_diem_velocyto_flts), ignore.case=TRUE, value=TRUE)}
if(!is.null(chloroplast_pattern)) {cp_genes <- grep(chloroplast_pattern, rownames(seur_diem_velocyto_flts), ignore.case=TRUE, value=TRUE)}

# remove CP/MT genes before clustering
seur_diem_velocyto_flts <- seur_diem_velocyto_flts[ !(rownames(seur_diem_velocyto_flts) %in% mt_genes), ]
seur_diem_velocyto_flts <- seur_diem_velocyto_flts[ !(rownames(seur_diem_velocyto_flts) %in% cp_genes), ]

cat("###CheckPoint###", length(mt_genes), "MT genes removed.\n")
cat("###CheckPoint###", length(cp_genes), "CP genes removed.\n")
cat("###CheckPoint### Number of genes left after CP/MT gene filtering:", nrow(seur_diem_velocyto_flts), "\n\n")
check_cells(seur_diem_velocyto_flts, "UMI/CP/MT/Gene filtering", opt$min_cells)

summary_dims[["After UMI/CP/MT/Gene Filtering"]] <- c("nCell"=ncol(seur_diem_velocyto_flts), "nGene"=nrow(seur_diem_velocyto_flts))

seur_diem_velocyto_flts$total_counts <- Matrix::colSums(GetAssayData(seur_diem_velocyto_flts, assay="RNA", layer="count"))
seur_diem_velocyto_flts$n_genes <- Matrix::colSums(GetAssayData(seur_diem_velocyto_flts, assay="RNA", layer="count")>0)
seur_diem_velocyto_flts$pct.mt <- PercentageFeatureSet(seur_diem_velocyto_flts, features = mt_genes, assay = "RNA")
seur_diem_velocyto_flts$pct.cp <- PercentageFeatureSet(seur_diem_velocyto_flts, features = cp_genes, assay = "RNA")
seur_diem_velocyto_flts$pct.rb <- PercentageFeatureSet(seur_diem_velocyto_flts, features = ribosomal_genes, assay = "RNA")

# Initial clustering ------------------------------------------------------

seur_diem_velocyto_flts <- SCTransform(seur_diem_velocyto_flts, variable.features.n = opt$nHVG, verbose = FALSE)
seur_diem_velocyto_flts <- RunPCA(seur_diem_velocyto_flts, verbose = FALSE)
seur_diem_velocyto_flts <- RunUMAP(seur_diem_velocyto_flts, dims = 1:30, verbose = FALSE)
seur_diem_velocyto_flts <- FindNeighbors(seur_diem_velocyto_flts, dims = 1:30, verbose = FALSE) # PCA reduction was used
seur_diem_velocyto_flts <- FindClusters(seur_diem_velocyto_flts, verbose = FALSE)

dp <- DimPlot(seur_diem_velocyto_flts, label = TRUE)
# ggsave(paste0("UMAP_diem_velocyto_flts_", opt$sample, ".png"), plot = dp, width = 6, height = 5, units = "in")
summary_plts[["After UMI/CP/MT/Gene Filtering"]] <- c("UMAP"=dp)

plot1 <- FeatureScatter(seur_diem_velocyto_flts, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- FeatureScatter(seur_diem_velocyto_flts, feature1 = "nCount_RNA", feature2 = "pct.mt")
plot3 <- FeatureScatter(seur_diem_velocyto_flts, feature1 = "nCount_RNA", feature2 = "pct.cp")
plot4 <- FeatureScatter(seur_diem_velocyto_flts, feature1 = "nCount_RNA", feature2 = "pct.rb")
plot5 <- FeatureScatter(seur_diem_velocyto_flts, feature1 = "nCount_RNA", feature2 = "score.debris")
plot6 <- FeatureScatter(seur_diem_velocyto_flts, feature1 = "nCount_RNA", feature2 = "unsplice_ratio")
sp <- (plot1 + plot2 + plot3)/(plot4 + plot5 + plot6)
# ggsave(paste0("QC_feature_scatter_diem_velocyto_flts_", opt$sample, ".png"), plot = sp, width = 20, height = 6, units = "in")
summary_plts[["After UMI/CP/MT/Gene Filtering"]] <- c(summary_plts[["After UMI/CP/MT/Gene Filtering"]], "Scatter"=sp)

# saveRDS(seur_diem_velocyto_flts, paste0("seur_diem_velocyto_flts_", opt$sample, ".rds"))

# Doublet detection -------------------------------------------------------

if(opt$remove_doublet){
  invisible(
    capture.output({
      sweep.res <- paramSweep(seur_diem_velocyto_flts, PCs = 1:30, sct = TRUE, num.cores = 1)
    })
  )
  sweep.stats <- summarizeSweep(sweep.res)
  best.pK <- as.numeric(as.character(find.pK(sweep.stats)$pK[which.max(find.pK(sweep.stats)$BCmetric)]))

  annotations <- seur_diem_velocyto_flts$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075*ncol(seur_diem_velocyto_flts))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  seur_diem_velocyto_flts <- doubletFinder(seur_diem_velocyto_flts, PCs = 1:30, pN = 0.25, pK = best.pK, nExp = nExp_poi, sct=TRUE)

  pANN_colname <- colnames(seur_diem_velocyto_flts@meta.data)[grep("pANN",colnames(seur_diem_velocyto_flts@meta.data))]
  seur_diem_velocyto_flts[["score.doublet"]] <- seur_diem_velocyto_flts[[pANN_colname]]
  seur_diem_velocyto_flts[[pANN_colname]] <- NULL

  seur_diem_velocyto_flts <- doubletFinder(seur_diem_velocyto_flts, PCs = 1:30, pN = 0.25, pK = best.pK, nExp = nExp_poi.adj, sct=TRUE, reuse.pANN = "score.doublet")
  seur_diem_velocyto_flts <- AddMetaData(seur_diem_velocyto_flts, metadata = list("doubletFinder"=seur_diem_velocyto_flts@meta.data[,ncol(seur_diem_velocyto_flts@meta.data)]))
  df_colname <- colnames(seur_diem_velocyto_flts@meta.data)[grep("DF\\.",colnames(seur_diem_velocyto_flts@meta.data))]
  for (c in df_colname) { seur_diem_velocyto_flts[[c]] <- NULL }

  doublet_frac <- round(table(seur_diem_velocyto_flts[["doubletFinder"]])[["Doublet"]]/ncol(seur_diem_velocyto_flts) * 100, digits = 1)
  cat("###CheckPoint### Final doublet fraction:", doublet_frac, "\n")
}else{
  seur_diem_velocyto_flts[["doubletFinder"]] <- "Singlet"
  seur_diem_velocyto_flts[["score.doublet"]] <- 0
}

summary_tbs[["After UMI/CP/MT/Gene Filtering"]] <- do.call(rbind, lapply(seur_diem_velocyto_flts@meta.data[c("total_counts","n_genes","pct.mt","pct.cp","pct.rb","score.debris", "unsplice_ratio", "score.doublet", "scDblFinder.score")], summary))
seur_objs[["After UMI/CP/MT/Gene Filtering"]] <- seur_diem_velocyto_flts

dp1 <- DimPlot(seur_diem_velocyto_flts, group.by = "doubletFinder")
dp2 <- FeaturePlot(seur_diem_velocyto_flts, features = "score.doublet")
dp <- dp1 + dp2
# ggsave(paste0("UMAP_doublets_diem_velocyto_flts_", opt$sample, ".png"), plot = dp, width = 10, height = 5, units = "in")
summary_plts[["After UMI/CP/MT/Gene Filtering"]] <- c(summary_plts[["After UMI/CP/MT/Gene Filtering"]], "UMAP-doubletFinderCall"=dp)

dp1 <- DimPlot(seur_diem_velocyto_flts, group.by = "scDblFinder.class")
dp2 <- FeaturePlot(seur_diem_velocyto_flts, features = "scDblFinder.score")
dp <- dp1 + dp2
# ggsave(paste0("UMAP_doublets_diem_velocyto_flts_", opt$sample, ".png"), plot = dp, width = 10, height = 5, units = "in")
summary_plts[["After UMI/CP/MT/Gene Filtering"]] <- c(summary_plts[["After UMI/CP/MT/Gene Filtering"]], "UMAP-scDblFinderCall"=dp)

sp <- FeatureScatter(seur_diem_velocyto_flts, feature1 = "score.doublet", feature2 = "scDblFinder.score")
summary_plts[["After UMI/CP/MT/Gene Filtering"]] <- c(summary_plts[["After UMI/CP/MT/Gene Filtering"]], "Scatter-doublets"=sp)

vp <- VlnPlot(seur_diem_velocyto_flts, assay = "RNA", layer = "counts", group.by = "orig.ident", features = c("total_counts", "n_genes", "pct.mt", "pct.cp", "pct.rb", "score.debris", "unsplice_ratio", "score.doublet", "scDblFinder.score"), pt.size = 0.01, alpha = 0.1, ncol = 9)
# ggsave(paste0("QC_feature_violin_diem_velocyto_flts_", opt$sample, ".png"), plot = vp, width = 15, height = 7, units = "in")
summary_plts[["After UMI/CP/MT/Gene Filtering"]] <- c(summary_plts[["After UMI/CP/MT/Gene Filtering"]], "Violin"=vp)

# check the nCount and nFeature between Singlets and Doublets
vp <- VlnPlot(seur_diem_velocyto_flts, assay = "RNA", layer = "counts", group.by = "doubletFinder", features = c("total_counts", "n_genes", "pct.mt", "pct.cp", "pct.rb", "score.debris", "unsplice_ratio", "score.doublet", "scDblFinder.score"), pt.size = 0.01, alpha = 0.1, ncol = 9)
# ggsave(paste0("QC_feature_violin_perDF_diem_velocyto_flts_", opt$sample, ".png"), plot = vp, width = 18, height = 7, units = "in")
summary_plts[["After UMI/CP/MT/Gene Filtering"]] <- c(summary_plts[["After UMI/CP/MT/Gene Filtering"]], "Violin-doubletFinderCall"=vp)

vp <- VlnPlot(seur_diem_velocyto_flts, assay = "RNA", layer = "counts", group.by = "scDblFinder.class", features = c("total_counts", "n_genes", "pct.mt", "pct.cp", "pct.rb", "score.debris", "unsplice_ratio", "score.doublet", "scDblFinder.score"), pt.size = 0.01, alpha = 0.1, ncol = 9)
# ggsave(paste0("QC_feature_violin_perDF_diem_velocyto_flts_", opt$sample, ".png"), plot = vp, width = 18, height = 7, units = "in")
summary_plts[["After UMI/CP/MT/Gene Filtering"]] <- c(summary_plts[["After UMI/CP/MT/Gene Filtering"]], "Violin-scDblFinderCall"=vp)

vp <- VlnPlot(seur_diem_velocyto_flts, assay = "RNA", layer = "counts", features = c("total_counts", "n_genes", "pct.mt", "pct.cp", "pct.rb", "score.debris", "unsplice_ratio", "score.doublet", "scDblFinder.score"), pt.size = 0.01, alpha = 0.1, ncol = 9)
# ggsave(paste0("QC_feature_violin_perCls_diem_velocyto_flts_", opt$sample, ".png"), plot = vp, width = 25, height = 7, units = "in")
summary_plts[["After UMI/CP/MT/Gene Filtering"]] <- c(summary_plts[["After UMI/CP/MT/Gene Filtering"]], "Violin-clusters"=vp)

### Remove doublets

if(opt$remove_doublet){
  cat("###CheckPoint###", "The number of doublets assigned by doubletFinder is", table(seur_diem_velocyto_flts[["doubletFinder"]])[["Doublet"]], "\n")
  if ("doublet" %in% names(table(seur_diem_velocyto_flts[["scDblFinder.class"]]))) {
    cat("###CheckPoint###", "The number of doublets assigned by scDblFinder is", table(seur_diem_velocyto_flts[["scDblFinder.class"]])[["doublet"]], "\n")
  } else {
    cat("###CheckPoint###", "The number of doublets assigned by scDblFinder is 0\n")
  }

  seur_diem_velocyto_flts_dblt <- subset(seur_diem_velocyto_flts, subset = doubletFinder == "Singlet")
  cat("###CheckPoint###", ncol(seur_diem_velocyto_flts_dblt@assays$RNA), " cells left after removing doubletFinder-assigned doublets. \n")
  seur_diem_velocyto_flts_dblt <- subset(seur_diem_velocyto_flts_dblt, subset = scDblFinder.class == "singlet")
  cat("###CheckPoint###", ncol(seur_diem_velocyto_flts_dblt@assays$RNA), " cells left after removing scDblFinder-assigned doublets. \n")

  seur_diem_velocyto_flts_dblt <- subset(seur_diem_velocyto_flts_dblt, subset = score.doublet <= opt$max_doublet_score)
  cat("###CheckPoint###", ncol(seur_diem_velocyto_flts_dblt@assays$RNA), " cells left after filtering by doublet scores. \n")
}else{
  seur_diem_velocyto_flts_dblt <- seur_diem_velocyto_flts
  cat("###CheckPoint###", ncol(seur_diem_velocyto_flts_dblt@assays$RNA), " cells left. Doublet removal was disabled. \n")
}
check_cells(seur_diem_velocyto_flts_dblt, "doublet removal", opt$min_cells)

# reclustering
seur_diem_velocyto_flts_dblt <- SCTransform(seur_diem_velocyto_flts_dblt, variable.features.n = opt$nHVG, verbose = FALSE)
seur_diem_velocyto_flts_dblt <- RunPCA(seur_diem_velocyto_flts_dblt, verbose = FALSE)
seur_diem_velocyto_flts_dblt <- RunUMAP(seur_diem_velocyto_flts_dblt, dims = 1:30, verbose = FALSE)
seur_diem_velocyto_flts_dblt <- FindNeighbors(seur_diem_velocyto_flts_dblt, dims = 1:30, verbose = FALSE) # PCA reduction was used
seur_diem_velocyto_flts_dblt <- FindClusters(seur_diem_velocyto_flts_dblt, verbose = FALSE)

dp <- DimPlot(seur_diem_velocyto_flts_dblt, label = TRUE)
# ggsave(paste0("UMAP_diem_velocyto_flts_dblt_", opt$sample, ".png"), plot = dp, width = 6, height = 5, units = "in")
summary_plts[["After Doublet Removal"]] <- c("UMAP"=dp)

vp <- VlnPlot(seur_diem_velocyto_flts_dblt, assay = "RNA", layer = "counts", group.by = "orig.ident", features = c("total_counts", "n_genes", "pct.mt", "pct.cp", "pct.rb", "score.debris", "unsplice_ratio", "score.doublet", "scDblFinder.score"), pt.size = 0.01, alpha = 0.1, ncol = 9)
summary_plts[["After Doublet Removal"]] <- c(summary_plts[["After Doublet Removal"]], "Violin"=vp)
vp <- VlnPlot(seur_diem_velocyto_flts_dblt, assay = "RNA", layer = "counts", features = c("total_counts", "n_genes", "pct.mt", "pct.cp", "pct.rb", "score.debris", "unsplice_ratio", "score.doublet", "scDblFinder.score"), pt.size = 0.01, alpha = 0.1, ncol = 9)
summary_plts[["After Doublet Removal"]] <- c(summary_plts[["After Doublet Removal"]], "Violin-clusters"=vp)

summary_dims[["After Doublet Removal"]] <- c("nCell"=ncol(seur_diem_velocyto_flts_dblt@assays$RNA), "nGene"=nrow(seur_diem_velocyto_flts_dblt@assays$RNA))
summary_tbs[["After Doublet Removal"]] <- do.call(rbind, lapply(seur_diem_velocyto_flts_dblt@meta.data[c("total_counts","n_genes","pct.mt","pct.cp","pct.rb","score.debris", "unsplice_ratio", "score.doublet", "scDblFinder.score")], summary))

# saveRDS(seur_diem_velocyto_flts_dblt, paste0("seur_diem_velocyto_flts_dblt_", opt$sample, ".rds"))
seur_objs[["After Doublet Removal"]] <- seur_diem_velocyto_flts_dblt

# Potential Seurat cluster of doublets ---------------------------

if(opt$remove_doublet){
  n_clusters <- length(unique(seur_diem_velocyto_flts_dblt$seurat_clusters))
  if (n_clusters >= 3) {
    dbl <- scDblFinder::findDoubletClusters(GetAssayData(seur_diem_velocyto_flts_dblt, layer = "counts"), seur_diem_velocyto_flts_dblt$seurat_clusters)
    write.csv(dbl, paste0("findDoubletClusters_", opt$sample, ".csv"))
  } else {
    cat("###CheckPoint### Only", n_clusters, "cluster(s) detected; skipping findDoubletClusters (requires >= 3).\n")
  }
}

# Filtering by number of cluster markers ----------------------------------

seur_diem_velocyto_flts_dblt <- PrepSCTFindMarkers(seur_diem_velocyto_flts_dblt)
cluster_markers <- FindAllMarkers(seur_diem_velocyto_flts_dblt, only.pos = T, assay = "SCT", logfc.threshold = log2(1.5), min.pct = 0.1)

if(nrow(cluster_markers)<opt$min_nClusterMarker){
  seur_diem_velocyto_flts_dblt_cls <- seur_diem_velocyto_flts_dblt
  cat("###CheckPoint###", "Too few markers detected. Skip filtering by number of cluster markers. \n")
  summary_plts[["After Min Cluster Marker Filtering"]] <- summary_plts[["After Doublet Removal"]]
  summary_dims[["After Min Cluster Marker Filtering"]] <- summary_dims[["After Doublet Removal"]]
  summary_tbs[["After Min Cluster Marker Filtering"]] <- summary_tbs[["After Doublet Removal"]]

  seur_objs[["After Min Cluster Marker Filtering"]] <- seur_diem_velocyto_flts_dblt_cls
}else{
  cat("###CheckPoint###", "The number of markers for each cluster:\n")
  print(table(cluster_markers$cluster))

  cluster_markers_flt <- cluster_markers %>%
    group_by(cluster) %>%
    filter(n() >= opt$min_nClusterMarker) %>%
    ungroup()

  removed_clusters <- setdiff(unique(cluster_markers$cluster), unique(cluster_markers_flt$cluster))

  if (length(removed_clusters)>0) {
    seur_diem_velocyto_flts_dblt_cls <- subset(seur_diem_velocyto_flts_dblt, subset = seurat_clusters %in% removed_clusters, invert = T)
    cat("###CheckPoint###", "Clusters removed due to too few markers:", paste(removed_clusters, collapse = ", "), "\n")
    cat("###CheckPoint### Number of cell barcodes left after low-marker cluster filtering:", ncol(seur_diem_velocyto_flts_dblt_cls@assays$RNA), "\n")
    cat("###CheckPoint### Number of genes left after low-marker cluster filtering:", nrow(seur_diem_velocyto_flts_dblt_cls@assays$RNA), "\n\n")
    check_cells(seur_diem_velocyto_flts_dblt_cls, "low-marker cluster filtering", opt$min_cells)

    # reclustering
    seur_diem_velocyto_flts_dblt_cls <- SCTransform(seur_diem_velocyto_flts_dblt_cls, variable.features.n = opt$nHVG, verbose = FALSE)
    seur_diem_velocyto_flts_dblt_cls <- RunPCA(seur_diem_velocyto_flts_dblt_cls, verbose = FALSE)
    seur_diem_velocyto_flts_dblt_cls <- RunUMAP(seur_diem_velocyto_flts_dblt_cls, dims = 1:30, verbose = FALSE)
    seur_diem_velocyto_flts_dblt_cls <- FindNeighbors(seur_diem_velocyto_flts_dblt_cls, dims = 1:30, verbose = FALSE) # PCA reduction was used
    seur_diem_velocyto_flts_dblt_cls <- FindClusters(seur_diem_velocyto_flts_dblt_cls, verbose = FALSE)

    dp <- DimPlot(seur_diem_velocyto_flts_dblt_cls, label = TRUE)
    summary_plts[["After Min Cluster Marker Filtering"]] <- c("UMAP"=dp)
    vp <- VlnPlot(seur_diem_velocyto_flts_dblt_cls, assay = "RNA", layer = "counts", group.by = "orig.ident", features = c("total_counts", "n_genes", "pct.mt", "pct.cp", "pct.rb", "score.debris", "unsplice_ratio", "score.doublet", "scDblFinder.score"), pt.size = 0.01, alpha = 0.1, ncol = 9)
    summary_plts[["After Min Cluster Marker Filtering"]] <- c(summary_plts[["After Min Cluster Marker Filtering"]], "Violin"=vp)
    vp <- VlnPlot(seur_diem_velocyto_flts_dblt_cls, assay = "RNA", layer = "counts", features = c("total_counts", "n_genes", "pct.mt", "pct.cp", "pct.rb", "score.debris", "unsplice_ratio", "score.doublet", "scDblFinder.score"), pt.size = 0.01, alpha = 0.1, ncol = 9)
    summary_plts[["After Min Cluster Marker Filtering"]] <- c(summary_plts[["After Min Cluster Marker Filtering"]], "Violin-clusters"=vp)

    summary_dims[["After Min Cluster Marker Filtering"]] <- c("nCell"=ncol(seur_diem_velocyto_flts_dblt_cls@assays$RNA), "nGene"=nrow(seur_diem_velocyto_flts_dblt_cls@assays$RNA))
    summary_tbs[["After Min Cluster Marker Filtering"]] <- do.call(rbind, lapply(seur_diem_velocyto_flts_dblt_cls@meta.data[c("total_counts","n_genes","pct.mt","pct.cp","pct.rb","score.debris", "unsplice_ratio", "score.doublet", "scDblFinder.score")], summary))

    seur_objs[["After Min Cluster Marker Filtering"]] <- seur_diem_velocyto_flts_dblt_cls
  } else {
    seur_diem_velocyto_flts_dblt_cls <- seur_diem_velocyto_flts_dblt
    cat("###CheckPoint###", "No Cluster removed due to too few markers.\n")
    summary_plts[["After Min Cluster Marker Filtering"]] <- summary_plts[["After Doublet Removal"]]
    summary_dims[["After Min Cluster Marker Filtering"]] <- summary_dims[["After Doublet Removal"]]
    summary_tbs[["After Min Cluster Marker Filtering"]] <- summary_tbs[["After Doublet Removal"]]

    seur_objs[["After Min Cluster Marker Filtering"]] <- seur_diem_velocyto_flts_dblt_cls
  }
}

# Filtering by low-UMI clusters -------------------------------------------

meta <- seur_diem_velocyto_flts_dblt_cls@meta.data

# pairwise Wilcoxon tests: compare each cluster vs all others
clusters <- unique(meta$seurat_clusters)

if(length(clusters)>1){
  test_results <- lapply(clusters, function(cl) {
    meta$group <- ifelse(meta$seurat_clusters == cl, cl, "others")
    w_nCount <- wilcox.test(nCount_RNA ~ group, data = meta)
    w_nFeature <- wilcox.test(nFeature_RNA ~ group, data = meta)
    data.frame(
      cluster = cl,
      p_nCount = w_nCount$p.value,
      p_nFeature = w_nFeature$p.value
    )
  }) %>% bind_rows()

  # adjust for multiple testing (FDR)
  test_results <- test_results %>%
    mutate(
      p_adj_nCount = p.adjust(p_nCount, method = "BH"),
      p_adj_nFeature = p.adjust(p_nFeature, method = "BH")
    )

  summary_stats <- meta %>%
    group_by(seurat_clusters) %>%
    summarise(
      median_nCount = median(nCount_RNA),
      median_nFeature = median(nFeature_RNA)
    )

  combined <- left_join(summary_stats, test_results, by = c("seurat_clusters" = "cluster"))

  cat("###CheckPoint###", "The Wilcoxon tests for each cluster:", "\n")
  print(combined)

  low_clusters <- combined %>%
    filter(
      p_adj_nCount < 0.001 | p_adj_nFeature < 0.001,          # significant
      median_nCount < median(median_nCount)*0.6 | median_nFeature < median(median_nFeature)*0.6   # below global median
    )

  cat("###CheckPoint###", "The Wilcoxon tests for low-median clusters:", "\n")
  if(nrow(low_clusters)>0) {
    seur_diem_velocyto_flts_dblt_cls <- subset(seur_diem_velocyto_flts_dblt_cls, subset = seurat_clusters %in% low_clusters$seurat_clusters, invert = T)
    cat("###CheckPoint###", "Clusters removed due to low median UMI/nGene:", paste(low_clusters$seurat_clusters, collapse = ", "), "\n")
    cat("###CheckPoint### Number of cell barcodes left after low-marker cluster filtering:", ncol(seur_diem_velocyto_flts_dblt_cls@assays$RNA), "\n")
    cat("###CheckPoint### Number of genes left after low-marker cluster filtering:", nrow(seur_diem_velocyto_flts_dblt_cls@assays$RNA), "\n\n")
    check_cells(seur_diem_velocyto_flts_dblt_cls, "low-median-UMI/nGene cluster filtering", opt$min_cells)

    # reclustering
    seur_diem_velocyto_flts_dblt_cls <- SCTransform(seur_diem_velocyto_flts_dblt_cls, variable.features.n = opt$nHVG, verbose = FALSE)
    seur_diem_velocyto_flts_dblt_cls <- RunPCA(seur_diem_velocyto_flts_dblt_cls, verbose = FALSE)
    seur_diem_velocyto_flts_dblt_cls <- RunUMAP(seur_diem_velocyto_flts_dblt_cls, dims = 1:30, verbose = FALSE)
    seur_diem_velocyto_flts_dblt_cls <- FindNeighbors(seur_diem_velocyto_flts_dblt_cls, dims = 1:30, verbose = FALSE) # PCA reduction was used
    seur_diem_velocyto_flts_dblt_cls <- FindClusters(seur_diem_velocyto_flts_dblt_cls, verbose = FALSE)

    dp <- DimPlot(seur_diem_velocyto_flts_dblt_cls, label = TRUE)
    summary_plts[["After low-median-UMI/nGene Cluster Filtering"]] <- c("UMAP"=dp)
    vp <- VlnPlot(seur_diem_velocyto_flts_dblt_cls, assay = "RNA", layer = "counts", group.by = "orig.ident", features = c("total_counts", "n_genes", "pct.mt", "pct.cp", "pct.rb", "score.debris", "unsplice_ratio", "score.doublet", "scDblFinder.score"), pt.size = 0.01, alpha = 0.1, ncol = 9)
    summary_plts[["After low-median-UMI/nGene Cluster Filtering"]] <- c(summary_plts[["After low-median-UMI/nGene Cluster Filtering"]], "Violin"=vp)
    vp <- VlnPlot(seur_diem_velocyto_flts_dblt_cls, assay = "RNA", layer = "counts", features = c("total_counts", "n_genes", "pct.mt", "pct.cp", "pct.rb", "score.debris", "unsplice_ratio", "score.doublet", "scDblFinder.score"), pt.size = 0.01, alpha = 0.1, ncol = 9)
    summary_plts[["After low-median-UMI/nGene Cluster Filtering"]] <- c(summary_plts[["After low-median-UMI/nGene Cluster Filtering"]], "Violin-clusters"=vp)

    summary_dims[["After low-median-UMI/nGene Cluster Filtering"]] <- c("nCell"=ncol(seur_diem_velocyto_flts_dblt_cls@assays$RNA), "nGene"=nrow(seur_diem_velocyto_flts_dblt_cls@assays$RNA))
    summary_tbs[["After low-median-UMI/nGene Cluster Filtering"]] <- do.call(rbind, lapply(seur_diem_velocyto_flts_dblt_cls@meta.data[c("total_counts","n_genes","pct.mt","pct.cp","pct.rb","score.debris", "unsplice_ratio", "score.doublet", "scDblFinder.score")], summary))
    seur_objs[["After low-median-UMI/nGene Cluster Filtering"]] <- seur_diem_velocyto_flts_dblt_cls

  }else{
    cat("###CheckPoint### No low-median-UMI/nGene cluster detected. \n\n")
    summary_plts[["After low-median-UMI/nGene Cluster Filtering"]] <- summary_plts[["After Min Cluster Marker Filtering"]]
    summary_dims[["After low-median-UMI/nGene Cluster Filtering"]] <- summary_dims[["After Min Cluster Marker Filtering"]]
    summary_tbs[["After low-median-UMI/nGene Cluster Filtering"]] <- summary_tbs[["After Min Cluster Marker Filtering"]]

    seur_objs[["After low-median-UMI/nGene Cluster Filtering"]] <- seur_diem_velocyto_flts_dblt_cls
  }
}else{
  cat("###CheckPoint### Only one cluster. Skip filtering by low-UMI clusters. \n\n")
  summary_plts[["After low-median-UMI/nGene Cluster Filtering"]] <- summary_plts[["After Min Cluster Marker Filtering"]]
  summary_dims[["After low-median-UMI/nGene Cluster Filtering"]] <- summary_dims[["After Min Cluster Marker Filtering"]]
  summary_tbs[["After low-median-UMI/nGene Cluster Filtering"]] <- summary_tbs[["After Min Cluster Marker Filtering"]]

  seur_objs[["After low-median-UMI/nGene Cluster Filtering"]] <- seur_diem_velocyto_flts_dblt_cls
}


# saveRDS(seur_diem_velocyto_flts_dblt_cls, paste0("seur_diem_velocyto_flts_dblt_cls_", opt$sample, ".rds"))
saveRDS(seur_diem_velocyto_flts_dblt_cls, paste0("seur_clean_", opt$sample, ".rds"))

saveRDS(summary_dims, paste0("summary_dims_", opt$sample, ".rds"))
saveRDS(summary_tbs, paste0("summary_tbs_", opt$sample, ".rds"))
saveRDS(summary_plts, paste0("summary_plts_", opt$sample, ".rds"))
saveRDS(seur_objs, paste0("seur_objs_", opt$sample, ".rds"))
