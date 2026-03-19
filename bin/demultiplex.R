#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library(optparse)))

option_list <- list(
  make_option(c("-i", "--inputDir"),
              type = "character",
              help = "Input directory which contains the raw_feature_bc_matrix directory [REQUIRED]"),
  
  make_option(c("-o", "--outDir"),
              type = "character",
              help = "Output directory [REQUIRED]"),
  
  make_option(c("-s", "--sample"),
              type = "character",
              help = "Sample name [REQUIRED]"),
  
  make_option(c("-S", "--species"),
              type = "character",
              help = "Species full name (comma-separated) [REQUIRED]", metavar = "character"),
  
  # Optional arguments with defaults
  
  make_option(c("--ref_yaml"),
              type = "character",
              default = "/project/gzy8899/qiaoshan/scRNAseq/references/scQC.yaml",
              help = "The reference organelle gene lists of all available species curated in yaml [default = %default]"),
  
  make_option(c("--min_UMI_per_cell_barcode"),
              type = "numeric",
              default = 400,
              help = "Minimum UMI per cell barcode [default = %default]"),
  
  make_option(c("--chisq_pvalues_max"),
              type = "numeric",
              default = 0.01,
              help = "Maximum pvalue in chi-squared test [default = %default]"),
  
  make_option(c("--ambient_rate_max"),
              type = "numeric",
              default = 0.5,
              help = "Maximum ambient rate [default = %default]"),

  make_option(c("--multiple_species_per_droplet"),
              type = "logical",
              default = FALSE,
              help = "Allow droplets to contribute to multiple species outputs [default = %default]"),
  
  make_option(c("--demux_matrix"),
              type = "logical",
              default = TRUE,
              help = "Write demultiplexed matrices to outDir/[species]/raw_feature_bc_matrix [default = %default]")

)

args <- parse_args(OptionParser(option_list=option_list))

# for testing
#options(
#  inputDir = "/project/gzy8899/qiaoshan/scRNAseq/experiments/control_exp_brassica-mix_shoot/cellranger/brassica-mix/brassica_mix_rep1/outs",
#  sample = "brassica_mix_rep1",
#  species = c("Arabidopsis thaliana","Arabidopsis lyrata","Brassica oleracea","Capsella rubella"),
#  outDir = "/project/gzy8899/qiaoshan/scRNAseq/experiments/control_exp_brassica-mix_shoot/demultiplex",
#  min_UMI_per_cell_barcode = 400,
#  chisq_pvalues_max = 0.01,
#  ambient_rate_max = 0.5
#)
#args <- options()


# init --------------------------------------------------------------------

required_args <- c("inputDir", "sample", "species", "outDir")

missing <- required_args[ sapply(required_args, function(x) is.null(args[[x]])) ]

if (length(missing) > 0) {
  cat("###ERROR### Missing required arguments:", paste(missing, collapse = ", "), "\n\n")
  print_help(OptionParser(option_list = option_list))
  quit(status = 1)
}

species_names <- strsplit(args$species, ",")[[1]]
outdir <- args$outDir
sample <- args$sample
min_UMI <- args$min_UMI_per_cell_barcode
chisq_pvalues_max <- args$chisq_pvalues_max
ambient_rate_max <- args$ambient_rate_max
multiple_species <- isTRUE(args$multiple_species_per_droplet)
demux_matrix <- isTRUE(args$demux_matrix)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(DropletUtils)
  library(yaml)
  library(scDblFinder)
})

matrix_dir <- paste(args$inputDir, "raw_feature_bc_matrix", sep = "/")

if (!dir.exists(matrix_dir)) {
  stop(sprintf("Matrix directory not found: %s", matrix_dir), call. = FALSE)
}

if(!dir.exists(outdir)) {dir.create(outdir, recursive = TRUE, showWarnings = FALSE)}
setwd(outdir)

# Load the raw feature-barcode matrix --------------------------------------
cat("Loading matrix from: ", matrix_dir, "\n", sep = "", file = stdout())
counts <- Read10X(data.dir = matrix_dir)

if (is.list(counts)) {
  counts <- counts[[1]]
  cat("Detected multiple assay entries; using the first element returned by Read10X().", "\n", sep = "", file = stdout())
}

if (!inherits(counts, "dgCMatrix")) {
  stop("Read10X did not return a sparse matrix. Please verify the input directory.", call. = FALSE)
}

feature_names <- rownames(counts)
cells <- colnames(counts)

# Plan step 1: gather per-species gene sets and background counts -----------
species_gene_list <- lapply(species_names, function(sp) {
  pattern <- paste0("^", gsub(" ", "_", sp), "_")
  genes <- feature_names[grepl(pattern, feature_names)]
  if (length(genes) == 0) {
    warning(sprintf("No features matched for species '%s'", sp), call. = FALSE)
  }
  genes
})
names(species_gene_list) <- species_names
cat("Step 1/6 complete: prepared per-species gene lists.", "\n", sep = "", file = stdout())

if(demux_matrix) {
  for (sp in species_names) {
    sp_dir_name <- gsub(" ", "_", sp)
    sp_dir <- file.path(outdir, sp_dir_name, sample, "raw_feature_bc_matrix")

    sp_counts <- counts[species_gene_list[[sp]],]
    rownames(sp_counts) <- gsub(sp_dir_name, "", rownames(sp_counts), fixed = TRUE)
    rownames(sp_counts) <- gsub("^_+", "", rownames(sp_counts))
    
    sp_seur_obj <- Seurat::CreateSeuratObject(
      counts = sp_counts,
      project = sample,
      assay = "RNA"
    )
    if(dir.exists(sp_dir)) { unlink(sp_dir, recursive = TRUE) }
    dir.create(dirname(sp_dir), recursive = TRUE, showWarnings = FALSE)
    write10xCounts(
      sp_dir,
      GetAssayData(sp_seur_obj, assay="RNA", layer="counts"),
      gene.id = rownames(sp_seur_obj),
      gene.symbol = rownames(sp_seur_obj),
      barcodes = colnames(sp_seur_obj),
      type = "sparse",
      version = "3"
    )
  }
}

species_counts <- matrix(
  0,
  nrow = length(cells),
  ncol = length(species_names),
  dimnames = list(cells, species_names)
)

# Plan step 2: per-cell UMI totals for each species -------------------------
for (sp in species_names) {
  species_counts[, sp] <- as.numeric(Matrix::colSums(counts[species_gene_list[[sp]], , drop = FALSE]))
}
cat("Step 2/6 complete: tallied UMIs per cell for each species.", "\n", sep = "", file = stdout())

# Background totals from plan step 1 ----------------------------------------
background_totals <- colSums(species_counts)
active_species <- background_totals > 0

if (!any(active_species)) {
  stop("No counts detected for the supplied species list.", call. = FALSE)
}

zero_species <- names(background_totals)[background_totals == 0]

if (length(zero_species) > 0) {
  warning(
    paste0(
      "Dropping species with zero detected counts:\n  - ",
      paste(zero_species, collapse = "\n  - ")
    ),
    call. = FALSE
  )
  
  species_counts <- species_counts[, background_totals > 0, drop = FALSE]
  species_names  <- colnames(species_counts)
  background_totals <- background_totals[species_names]
}

if (length(species_names) < 2) {
  stop(
    paste(
      "Only one species with detected counts remains after filtering;",
      "this script requires multiple species to perform chi-squared testing."
    ),
    call. = FALSE
  )
}

background_props <- background_totals / sum(background_totals)
cat("Proportion total UMIs in each species: ", "\n", sep = "", file = stdout())
print(background_props)
cat("Step 3/6 setup complete: computed background species proportions.", "\n", sep = "", file = stdout())

# Plan step 3: per-barcode chi-squared test against background --------------
chisq_stat <- rep(NA_real_, nrow(species_counts))
chisq_pvalues <- rep(NA_real_, nrow(species_counts))

for (i in seq_len(nrow(species_counts))) { # for each cell
  observed <- species_counts[i,]
  total_observed <- sum(observed)
  if (total_observed < min_UMI) {
    next
  }
  expected <- total_observed * background_props
  if (any(expected == 0)) {
    stop("Expected UMIs in cell i have zero.")
  }
  chisq_stat[i] <- sum((observed - expected)^2 / expected)
  df <- length(observed) - 1
  if (df <= 0) {
    next
  }
  chisq_pvalues[i] <- stats::pchisq(chisq_stat[i], df = df, lower.tail = FALSE)
}

keep_cells <- !is.na(chisq_pvalues) & chisq_pvalues < chisq_pvalues_max

cat("The number of chi-squared tested and kept droplets:", "\n", sep = "", file = stdout())
print(table(
  is_tested = !is.na(chisq_pvalues),
  is_kept   = keep_cells
))
cat("Step 3/6 complete: applied chi-squared filter.", "\n", sep = "", file = stdout())

# Plan step 4: identify clean nuclei via 99th percentile thresholds ----------
species_thresholds <- vapply(species_names, function(sp) {
  counts_sp <- species_counts[, sp]
  positive_counts <- counts_sp[counts_sp > 0]
  if (length(positive_counts) == 0) {
    return(Inf)
  }
  stats::quantile(positive_counts, probs = 0.99, type = 7, names = FALSE)
}, numeric(1))

threshold_matrix <- matrix(
  species_thresholds,
  nrow = nrow(species_counts),
  ncol = ncol(species_counts),
  byrow = TRUE,
  dimnames = dimnames(species_counts)
)

true_mask <- species_counts >= threshold_matrix & species_counts > 0
true_counts <- species_counts
true_counts[!true_mask] <- 0
cat(
  "Step 4/6 complete: derived clean nuclei masks via 99th percentile thresholds.",
  "\n",
  sep = "",
  file = stdout()
)
cat(
  paste0(
    "Thresholds (UMIs) per species: ",
    paste(
      paste0(names(species_thresholds), ": ", signif(species_thresholds, 4)),
      collapse = "; "
    )
  ),
  "\n",
  sep = "",
  file = stdout()
)

# Plan step 5: estimate ambient RNA contribution per barcode ----------------
sum_true_counts <- rowSums(true_counts)
total_umi <- as.numeric(Matrix::colSums(counts))
ambient_counts <- total_umi - sum_true_counts
ambient_rate <- ifelse(total_umi > 0, ambient_counts / total_umi, NA_real_)
cat("Step 5/6 complete: calculated ambient RNA rates per barcode.", "\n", sep = "", file = stdout())

metadata <- data.frame(
  cell = cells,
  total_umi = total_umi,
  ambient_counts = ambient_counts,
  ambient_rate = ambient_rate,
  chisq_stat = chisq_stat,
  chisq_pvalue = chisq_pvalues,
  stringsAsFactors = FALSE
)

for (sp in species_names) {
  metadata[[paste0("umi_", gsub(" ", "_", sp))]] <- species_counts[, sp]
  metadata[[paste0("clean_", gsub(" ", "_", sp))]] <- true_mask[, sp]
  metadata[[paste0("threshold_", gsub(" ", "_", sp))]] <- species_thresholds[[sp]]
}

multi_species_counts <- rowSums(true_mask[, species_names, drop = FALSE])
species_passed_label <- vapply(
  seq_len(nrow(true_mask)),
  function(i) {
    mask_row <- true_mask[i, species_names, drop = FALSE]
    species_i <- species_names[mask_row[1, ]]
    if (length(species_i) == 0) "" else paste(species_i, collapse = ";")
  },
  character(1)
)
metadata$species_passed <- species_passed_label

cat("Clean nuclei counts per barcode (number of species meeting threshold):", "\n", sep = "", file = stdout())
print(table(num_species = multi_species_counts))

if (multiple_species) {
  metadata$dominant_species <- NA_character_
  cat("Multiple-species mode enabled: retaining all species passing thresholds per barcode.", "\n", sep = "", file = stdout())
} else {
  normalized_counts <- sweep(species_counts, 2, background_props[colnames(species_counts)], FUN = "/")
  normalized_counts[!is.finite(normalized_counts)] <- 0
  dominant_idx <- max.col(normalized_counts, ties.method = "first")
  dominant_species <- species_names[dominant_idx]
  dominant_species[rowSums(species_counts) < min_UMI] <- NA_character_
  metadata$dominant_species <- dominant_species
  cat("Assigned dominant species based on background-normalized UMIs.", "\n", sep = "", file = stdout())
}

rownames(metadata) <- metadata$cell

# Plan step 6: remove high ambient barcodes and keep clean nuclei ------------
keep_cells2 <- !is.na(metadata$ambient_rate) & metadata$ambient_rate <= ambient_rate_max
keep_both <- keep_cells & keep_cells2
cat("Step 6/6 initialized: filtered barcodes by ambient RNA rate threshold.", "\n", sep = "", file = stdout())

cat("Post-filter clean nuclei counts (kept barcodes per number of species):", "\n", sep = "", file = stdout())
print(table(num_species = multi_species_counts[keep_both]))

# Output: build Seurat objects per species after filtering ------------------
for (sp in species_names) {
  genes <- species_gene_list[[sp]]
  if (length(genes) == 0) {
    cat(sprintf("Skipping species '%s' because no matching features were found.", sp), "\n", sep = "", file = stdout())
    next
  }
  
  if (multiple_species) {
    cells_for_species <- rownames(metadata)[keep_both & true_mask[, sp]]
  } else {
    cells_for_species <- rownames(metadata)[keep_both & metadata$dominant_species == sp]
  }
  if (length(cells_for_species) == 0) {
    cat(sprintf("No cells passed filters for species '%s'; no RDS written.", sp), "\n", sep = "", file = stdout())
    next
  }

  counts_subset <- counts[genes, cells_for_species, drop = FALSE]
  species_prefix <- gsub(" ", "_", sp)
  rownames(counts_subset) <- gsub(species_prefix, "", rownames(counts_subset), fixed = TRUE)
  rownames(counts_subset) <- gsub("^_+", "", rownames(counts_subset))
  metadata_subset <- metadata[cells_for_species, , drop = FALSE]
  metadata_subset$assigned_species <- sp

  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts_subset,
    project = sample,
    assay = "RNA",
    meta.data = metadata_subset,
    min.cells = 0,
    min.features = 0
  )

  qc_patterns <- read_yaml(args$ref_yaml)
  if (sp %in% names(qc_patterns)) {
    if(!is.null(qc_patterns[[sp]]$mitochondrial_pattern)) {
      mitochondrial_pattern = qc_patterns[[sp]]$mitochondrial_pattern
    }else{
      mitochondrial_pattern = NULL
      mt_genes = readLines(qc_patterns[[sp]]$mitochondrial_genes)
    }
    if(!is.null(qc_patterns[[sp]]$chloroplast_pattern)) {
      chloroplast_pattern = qc_patterns[[sp]]$chloroplast_pattern
    }else{
      chloroplast_pattern = NULL
      cp_genes = readLines(qc_patterns[[sp]]$chloroplast_genes)
    }
    ribosomal_genes = readLines(qc_patterns[[sp]]$ribosomal_genes)
  } else {
    stop("###ERROR### This species is absent from the reference yaml. Please add mitochondrial/chloroplast/ribosomal genes to ref.yaml\n")
  }
  
  saveRDS(args, paste0(outdir, "/", species_prefix, "/summary_opts_", args$sample, ".rds"))
  
  summary_tbs <- list()
  summary_dims <- list()
  summary_plts <- list()
  seur_objs <- list()
  
  set.seed(123)
  
  seur_diem<-seurat_obj
  summary_dims[["After DIEM Debris Removal"]] <- c("nCell"=ncol(seur_diem), "nGene"=nrow(seur_diem))
  
  seur_diem$total_counts <- seur_diem[[paste0("umi_", gsub(" ", "_", sp))]]
  seur_diem$n_genes <- Matrix::colSums(Seurat::GetAssayData(seur_diem, assay = "RNA", layer = "counts") > 0)
  if(!is.null(mitochondrial_pattern)) {mt_genes <- grep(pattern=mitochondrial_pattern, x=rownames(seur_diem), ignore.case=TRUE, value=TRUE)}
  if(!is.null(chloroplast_pattern)) {cp_genes <- grep(pattern=chloroplast_pattern, x=rownames(seur_diem), ignore.case=TRUE, value=TRUE)}
  
  seur_diem$pct.mt <- PercentageFeatureSet(seur_diem, features = mt_genes, assay = "RNA")
  seur_diem$pct.cp <- PercentageFeatureSet(seur_diem, features = cp_genes, assay = "RNA")
  seur_diem$pct.rb <- PercentageFeatureSet(seur_diem, features = ribosomal_genes, assay = "RNA")
  seur_diem$score.debris <- seur_diem$ambient_rate
  
  summary_tbs[["After DIEM Debris Removal"]] <- do.call(rbind, lapply(seur_diem@meta.data[c("total_counts","n_genes","pct.mt","pct.cp","pct.rb","score.debris")], summary))
  
  Project(seur_diem) <- args$sample
  seur_diem$orig.ident <- args$sample
  
  sce_diem <- scDblFinder(GetAssayData(seur_diem, layer="counts"))
  seur_diem <- AddMetaData(seur_diem, metadata = as.data.frame(colData(sce_diem)))
  
  saveRDS(summary_dims, paste0(outdir, "/", species_prefix, "/summary_dims_", args$sample, ".rds"))
  saveRDS(summary_tbs, paste0(outdir, "/", species_prefix, "/summary_tbs_", args$sample, ".rds"))
  saveRDS(summary_plts, paste0(outdir, "/", species_prefix, "/summary_plts_", args$sample, ".rds"))
  
  seur_diem_bc <- Cells(seur_diem)
  seur_diem <- RenameCells(seur_diem, new.names = gsub("-1$", "", colnames(seur_diem)))
  seur_objs[["After DIEM Debris Removal"]] <- seur_diem
  saveRDS(seur_objs, paste0(outdir, "/", species_prefix, "/seur_objs_", args$sample, ".rds"))
  
  con <- gzfile(paste0(outdir, "/", species_prefix, "/seur_diem_barcodes_", args$sample, ".tsv.gz"), "w")
  writeLines(seur_diem_bc, con)
  close(con)
  
  cat(sprintf(
    "Saved %s Seurat object with %d cells to %s",
    sp,
    length(cells_for_species),
    outdir
  ), "\n", sep = "", file = stdout())
}

cat("Demultiplexing finished: all species processed.", "\n", sep = "", file = stdout())
