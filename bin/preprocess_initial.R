rm(list = ls())

suppressWarnings(suppressPackageStartupMessages(library(optparse)))

option_list <- list(
  make_option(c("-i", "--inputDir"),
              type = "character",
              help = "Input directory which contains the raw_feature_bc_matrix directory [REQUIRED]"),

  make_option(c("-s", "--sample"),
              type = "character",
              help = "Sample name [REQUIRED]"),

  make_option(c("-S", "--species"),
              type = "character",
              help = "Species full name [REQUIRED]"),

  make_option(c("-o", "--outDir"),
              type = "character",
              help = "Output directory [REQUIRED]"),

  make_option(c("--ref_yaml"),
              type = "character",
              default = "/project/gzy8899/qiaoshan/scRNAseq/references/scQC.yaml",
              help = "The reference organelle gene lists of all available species curated in yaml [default = %default]"),

  make_option(c("--min_diem_debris_score"),
              type = "numeric",
              default = 1,
              help = "Minimum DIEM-derived debris score [default = %default]"),

  make_option(c("--threads"),
              type = "integer",
              default = 1,
              help = "Threads [default = %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# init --------------------------------------------------------------------

required_args <- c("inputDir", "sample", "species", "outDir")

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
suppressWarnings(suppressPackageStartupMessages(library(diem)))
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(Matrix)))
suppressWarnings(suppressPackageStartupMessages(library(inflection)))
suppressWarnings(suppressPackageStartupMessages(library(ggpubr)))
suppressWarnings(suppressPackageStartupMessages(library(scDblFinder)))

cat("###CheckPoint### Arguments\n")
cat("###CheckPoint### ---------\n")
cat("###CheckPoint### Input dir  :", opt$inputDir, "\n")
cat("###CheckPoint### Sample     :", opt$sample, "\n")
cat("###CheckPoint### Species    :", opt$species, "\n")
cat("###CheckPoint### Output dir :", opt$outDir, "\n")
cat("###CheckPoint### ref_yaml         :", opt$ref_yaml, "\n")
cat("###CheckPoint### min_diem_debris_score        :", opt$min_diem_debris_score, "\n")
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

summary_tbs <- list()
summary_dims <- list()
summary_plts <- list()
seur_objs <- list()

set.seed(123)

# Diem debris filtering -------------------------------------------------------------------

auto_knee <- function(sce, plot = TRUE) {
  # Get total counts from DIEM object
  d <- droplet_data(sce)$total_counts
  d <- sort(d, decreasing = TRUE)

  # x = rank; y = counts
  x <- seq_along(d)
  y <- d

  # Detect knee (2nd derivative detection)
  kne <- inflection::uik(x, y)

  if (plot) {
    df <- data.frame(x = x, y = y)
    p <- ggplot(df, aes(x, y)) +
      geom_line() +
      geom_vline(xintercept = kne, color = "red") +
      scale_y_log10() +
      scale_x_log10() +
      ggtitle(paste("Knee index =", kne)) +
      xlab("Rank") +
      ylab("UMI count per cell")
    # print(p)
  }

  return(kne)
}

countsDir <- paste(opt$inputDir, "raw_feature_bc_matrix", sep = "/")

counts <- read_10x(countsDir) # Read 10X data into sparse matrix
sce <- create_SCE(counts) # Create SCE object from counts

cat("###CheckPoint### Raw Input\n")
cat("###CheckPoint### ---------\n")
cat("###CheckPoint### Number of cell barcodes in raw data:", ncol(sce), "\n")
cat("###CheckPoint### Number of genes in raw data:", nrow(sce), "\n\n")

summary_dims[["Cellranger Raw"]] <- c("nCell"=ncol(sce), "nGene"=nrow(sce))

if(!is.null(mitochondrial_pattern)) {mt_genes <- grep(pattern=mitochondrial_pattern, x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)}
if(!is.null(chloroplast_pattern)) {cp_genes <- grep(pattern=chloroplast_pattern, x=rownames(sce@gene_data), ignore.case=TRUE, value=TRUE)}
sce <- get_gene_pct(x = sce, genes=mt_genes, name="pct.mt")
sce <- get_gene_pct(x = sce, genes=cp_genes, name="pct.cp")
sce <- get_gene_pct(x = sce, genes=ribosomal_genes, name="pct.rb")

summary_tbs[["Cellranger Raw"]] <- do.call(rbind, lapply(sce@droplet_data, summary))

top_n <- auto_knee(sce)
cat("###CheckPoint### Number of test droplets assigned to DIEM:", top_n, "\n")

brp <- barcode_rank_plot(sce, ret = T) + geom_vline(xintercept = top_n, color = "red") +
  ggtitle(paste("Knee index =", top_n))
# ggsave(paste0("diem_barcode_rank_plot_", opt$sample, ".png"), brp, width = 6, height = 6, units = "in")

summary_plts[["Cellranger Raw"]] <- c("Barcode Rank"=brp)

diem_min_counts <- 200
if (sum(sce@droplet_data$total_counts>=diem_min_counts) < 2000) {
  diem_min_counts <- quantile(sce@droplet_data$total_counts, 1-top_n/ncol(sce))
}
cat("###CheckPoint### Min UMI counts for DIEM set to:", diem_min_counts, "\n")

sce <- set_debris_test_set(sce, min_counts = diem_min_counts, top_n = top_n)
sce <- filter_genes(sce)

diem_min_genes <- 200

warn_msg <- NULL
error_msg <- NULL
sce_tmp <- NULL

sce_tmp <- tryCatch(
  withCallingHandlers(
    get_pcs(sce, min_genes = diem_min_genes),
    warning = function(w) {
      warn_msg <<- conditionMessage(w)
      invokeRestart("muffleWarning")   # prevents printing the warning
    }
  ),
  error = function(e) {
    error_msg <<- conditionMessage(e)
    return(NULL)   # or your fallback behavior
  }
)

if (!is.null(warn_msg)|!is.null(error_msg)) {
  while (!is.null(warn_msg)|!is.null(error_msg)) {
    diem_min_genes <- diem_min_genes-5
    warn_msg <- NULL
    error_msg <- NULL
    sce_tmp <- tryCatch(
      withCallingHandlers(
        get_pcs(sce, min_genes = diem_min_genes),
        warning = function(w) {
          warn_msg <<- conditionMessage(w)
          invokeRestart("muffleWarning")   # prevents printing the warning
        }
      ),
      error = function(e) {
        error_msg <<- conditionMessage(e)
        return(NULL)   # or your fallback behavior
      }
    )
  }
}
sce <- sce_tmp

cat("###CheckPoint### Min gene counts for DIEM set to:", diem_min_genes, "\n\n")

sce <- init(sce)
sce <- run_em(sce, threads = opt$threads)
sce <- assign_clusters(sce)
sce <- estimate_dbr_score(sce, thresh_genes = diem_min_genes)

# Evaluate debris scores
# sm <- summarize_clusters(sce)

# saveRDS(sce, file = paste0("sce_", opt$sample, ".rds"))

p <- plot_clust(sce, feat_x = "total_counts", feat_y = "score.debris", log_x = TRUE, log_y = FALSE, ret = T)
# ggsave(paste0("diem_clust_plot_total_counts_", opt$sample, ".png"), p, width = 5, height = 5, units = "in")
summary_plts[["Cellranger Raw"]] <- c(summary_plts[["Cellranger Raw"]], "Scatter nUMI x debris.score"=p)

p <- plot_clust(sce, feat_x = "n_genes", feat_y = "score.debris", log_x = TRUE, log_y = FALSE, ret = T)
# ggsave(paste0("diem_clust_plot_n_genes_", opt$sample, ".png"), p, width = 5, height = 5, units = "in")
summary_plts[["Cellranger Raw"]] <- c(summary_plts[["Cellranger Raw"]], "Scatter nGene x debris.score"=p)

p <- plot_clust(sce, feat_x = "pct.mt", feat_y = "score.debris", log_x = TRUE, log_y = FALSE, ret = T)
# ggsave(paste0("diem_clust_plot_pct_mt_", opt$sample, ".png"), p, width = 5, height = 5, units = "in")
summary_plts[["Cellranger Raw"]] <- c(summary_plts[["Cellranger Raw"]], "Scatter pctMT x debris.score"=p)

p <- plot_clust(sce, feat_x = "pct.cp", feat_y = "score.debris", log_x = TRUE, log_y = FALSE, ret = T)
# ggsave(paste0("diem_clust_plot_pct_cp_", opt$sample, ".png"), p, width = 5, height = 5, units = "in")
summary_plts[["Cellranger Raw"]] <- c(summary_plts[["Cellranger Raw"]], "Scatter pctCP x debris.score"=p)

p <- plot_clust(sce, feat_x = "pct.rb", feat_y = "score.debris", log_x = TRUE, log_y = FALSE, ret = T)
# ggsave(paste0("diem_clust_plot_pct_rb_", opt$sample, ".png"), p, width = 5, height = 5, units = "in")
summary_plts[["Cellranger Raw"]] <- c(summary_plts[["Cellranger Raw"]], "Scatter pctRB x debris.score"=p)

# Call targets using debris score for single-nucleus data
sce <- call_targets(sce, clusters = "debris", thresh = NULL, min_genes = quantile(sce@droplet_data$n_genes, 0.8))
sce_metadata <- droplet_data(sce)

cat("###CheckPoint### The number of clean cells and debris called by DIEM: \n")
print(table(sce_metadata$Call))

cat("###CheckPoint### Among all cells/nuclei candidates (", nrow(sce_metadata), "), ", sum(sce_metadata$score.debris>opt$min_diem_debris_score), " have a debris score above ", opt$min_diem_debris_score, ". Among them, the number of cells/nuclei and debris called by DIEM are: \n", sep = "")
print(table(sce_metadata$Call[sce_metadata$score.debris>opt$min_diem_debris_score]))

# saveRDS(sce_metadata, file = paste0("sce_metadata_", opt$sample, ".rds"))

seur_diem <- convert_to_seurat(sce)

seur_diem <- subset(seur_diem, subset = score.debris <= opt$min_diem_debris_score)

cat("###CheckPoint### Number of cell barcodes left after removing debris:", ncol(seur_diem), "\n")
cat("###CheckPoint### Number of genes left after removing debris:", nrow(seur_diem), "\n\n")

summary_dims[["After DIEM Debris Removal"]] <- c("nCell"=ncol(seur_diem), "nGene"=nrow(seur_diem))
summary_tbs[["After DIEM Debris Removal"]] <- do.call(rbind, lapply(seur_diem@meta.data[c("total_counts","n_genes","pct.mt","pct.cp","pct.rb","score.debris")], summary))

Project(seur_diem) <- opt$sample
seur_diem$orig.ident <- opt$sample

sce_diem <- scDblFinder(GetAssayData(seur_diem, assay = "RNA", layer = "counts"))
seur_diem <- AddMetaData(seur_diem, metadata = as.data.frame(colData(sce_diem)))

seur_objs[["After DIEM Debris Removal"]] <- seur_diem

# saveRDS(seur_diem, paste0("seur_diem_", opt$sample, ".rds"))
saveRDS(summary_dims, paste0("summary_dims_", opt$sample, ".rds"))
saveRDS(summary_tbs, paste0("summary_tbs_", opt$sample, ".rds"))
saveRDS(summary_plts, paste0("summary_plts_", opt$sample, ".rds"))
saveRDS(seur_objs, paste0("seur_objs_", opt$sample, ".rds"))

seur_diem_bc <- gsub("$", "-1", Cells(seur_diem))

con <- gzfile(paste0("seur_diem_barcodes_", opt$sample, ".tsv.gz"), "w")
writeLines(seur_diem_bc, con)
close(con)

cat("###CheckPoint### DIEM debris filtering is completed. \n\n")
