# snRNAseq Processing Pipeline — User Guide

A Nextflow pipeline for processing single-nucleus RNA sequencing data from plants.

## Table of Contents
- [Bash Environment Preset](#bash-environment-preset)
- [Quick Start](#quick-start)
- [Parameters](#parameters)
- [Output Structure](#output-structure)
- [Pipeline Steps](#pipeline-steps)
- [Integration](#integration)
- [Q&A](#qa)
  - [What if I want to use different cutoff values for different samples?](#what-if-i-want-to-use-different-cutoff-values-for-different-samples)
  - [How can I build one summary report across samples from different experiments?](#how-can-i-build-one-summary-report-across-samples-from-different-experiments)

---

## Bash Environment Preset

Use the following command to add the module environment to your `.bashrc` file:

```bash
echo "module use /project/gzy8899/modules" >> ~/.bashrc
```

Then source the `.bashrc` file:

```bash
source ~/.bashrc
```

Check if the module environment is added:

```bash
module avail
```

If you see the `/project/gzy8899/modules` listed, the module environment is added successfully.

Require modules for this pipeline: 
- `nextflow/25.10.4`
- `cellranger/9.0.1`
- `agat/1.6.1`
- `scrnaseq/1.0`
- `velocyto/0.17.17`
- `samtools/1.22.1`
- `deeptools/3.5.6`

*This setup is needed to be done only once.*

---

## Quick Start

For single-species:

```bash
module load nextflow/25.10.4
nextflow run jow30/ROMalleyLab-scRNAseq-nf \
  --input samplesheet.csv \
  --out exp1-nf-results \
  --species "Arabidopsis thaliana"
```

For multi-species:

```bash
module load nextflow/25.10.4
nextflow run jow30/ROMalleyLab-scRNAseq-nf \
  --input samplesheet.csv \
  --out exp2-nf-results \
  --species "Arabidopsis thaliana,Capsella rubella"
```

Resume a run after an interruption:

```bash
nextflow run jow30/ROMalleyLab-scRNAseq-nf \
  --input samplesheet.csv \
  --out exp1-nf-results \
  {$otherParams} \
  -resume
```

---

Print help message:

```bash
nextflow run jow30/ROMalleyLab-scRNAseq-nf --help
```

---

## Parameters

### Mandatory Parameters

| Parameter  | Description                  |
|------------|------------------------------|
| `--input`  | Path to samplesheet CSV      |
| `--out`    | Output directory path        |
| `--species`| Species name(s), comma-separated (default: `Arabidopsis thaliana`) |

Samplesheet should be a CSV file with the following columns:

| Column      | Description                                          |
|-------------|------------------------------------------------------|
| `sample`    | Sample name (letters, numbers, `-`, `_` only)        |
| `fastq_dir` | Path to the directory containing FASTQ files         |
| `fastq_R1`  | Filename of the R1 FASTQ file              |
| `fastq_R2`  | Filename of the R2 FASTQ file              |

Example `samplesheet.csv`:

```
sample,fastq_dir,fastq_R1,fastq_R2
exp,/project/gzy8899/data/202510_DAP_10X,RO1_S1_L001_R1_001.fastq.gz,RO1_S1_L001_R2_001.fastq.gz
ctrl,/project/gzy8899/data/202510_DAP_10X,RO2_S2_L001_R1_001.fastq.gz,RO2_S2_L001_R2_001.fastq.gz
```

FASTQ filenames should follow standard Illumina naming: `*_S*_L*_R*_*.fastq.gz`.

### Species Configuration

For the given species in `--species`, the pipeline will search for existing reference files in [`refs`](refs) and use them for read mapping, QC filtering and cell-type annotation. 

The following species are pre-configured in [`refs`](refs):
- `Arabidopsis thaliana`
- `Arabidopsis lyrata`
- `Capsella rubella`
- `Brassica oleracea`

The configured species references include:
- genome assemblies in FASTA format
- genome annotations in GTF format
- CellRanger index folders
- **A yaml file** which contains organellar and ribosomal gene lists for QC filtering and cell-type annotation reference atlas/markers

By default, the pipeline reads the file [`refs/scQC.yaml`](refs/scQC.yaml) to identify organellar and ribosomal genes for each species given in `--species` during QC filtering. Seurat reference atlas and celltype marker references are also included in this file for automated cell-type annotation. ***Edit this file to update organelle gene lists, add new celltype references, or when you [add a new species](#add-a-new-species).***

#### Config the yaml file

For each species in [`refs/scQC.yaml`](refs/scQC.yaml), use **either** the pattern fields **or** the gene-list fields for mitochondria/chloroplast — not both. Ribosomal genes always require a file. Provide plain-text files with one gene ID per line. The `annotation_ref_seurat_obj` field and `celltype_markers` field are optional; if omitted, cell-type annotation is skipped.

| Field | Type | Description |
|-------|------|-------------|
| `mitochondrial_pattern` | regex | Grep pattern to match mitochondrial gene IDs (e.g. `"^ATMG"`) |
| `chloroplast_pattern` | regex | Grep pattern to match chloroplast gene IDs (e.g. `"^ATCG"`) |
| `mitochondrial_genes` | file path | One mitochondrial gene ID per line |
| `chloroplast_genes` | file path | One chloroplast gene ID per line |
| `ribosomal_genes` | file path | One ribosomal gene ID per line |
| `annotation_ref_seurat_obj` | file path | Seurat reference RDS with a `celltype` column for anchor-based label transfer |
| `celltype_markers` | file path | CSV file with `gene`, `name`, `clusterName`, `p_val_adj`, `avg_log2FC` columns for celltype marker genes |

Please make sure the gene IDs in the gene-list files are consistent with the gene IDs in the genome annotation file.

Example 1:
```yaml
Arabidopsis thaliana:
    mitochondrial_pattern: "^ATMG"
    chloroplast_pattern: "^ATCG"
    ribosomal_genes: "/path/to/ribosomal_genes.txt"
    annotation_ref_seurat_obj: "/path/to/reference_atlas.RDS"   # optional
    celltype_markers: "/path/to/celltype_markers.csv"           # optional
```

Example 2:
```yaml
Brassica oleracea:
    mitochondrial_genes: "/path/to/mitochondrial_genes.txt"
    chloroplast_genes: "/path/to/chloroplast_genes.txt"
    ribosomal_genes: "/path/to/ribosomal_genes.txt"
```

#### Add a new species

To configure a new species, add an entry under `params.species_map` in [nextflow.config](nextflow.config):

```groovy
'My new species': [
    genome:      '/path/to/genome.fa',
    gtf:         '/path/to/annotation.gtf',
    cellranger:  "${projectDir}/refs/My_new_species/cellranger"  // or null to build on-the-fly
]
```

If `cellranger` points to an existing directory, it will be used directly. Otherwise, the pipeline builds it from `genome` and `gtf`.

#### Configure multi-species references

Pre-built multi-species references exist for:
- `Arabidopsis thaliana,Capsella rubella`
- `Arabidopsis thaliana,Arabidopsis lyrata,Capsella rubella,Brassica oleracea`

To configure a new multi-species reference, add an entry under `params.combined_ref_map` in [nextflow.config](nextflow.config):

```groovy
'My new species A,My new species B': "${projectDir}/refs/My_new_species_A_B/cellranger"
```

The pipeline will build the combined reference from the given species references in `params.species_map` and save it to the path specified in `params.combined_ref_map`. If the path already exists, the pipeline will use it directly.

### Optional Parameters

| Parameter  | Description                          |
|------------|--------------------------------------|
| `--genome` | Path to genome FASTA file            |
| `--gtf`    | Path to annotation file (`.gtf`, `.gff3`, or `.gff`) |
| `--ref_yaml` | Path to scQC YAML which contains information about organelles and annotation references (default: `refs/scQC.yaml`) |

**Library QC:**

| Parameter                | Default | Description                                                  |
|--------------------------|---------|--------------------------------------------------------------|
| `--bin_size`             | `500`   | Window size to count reads for mapping efficiency analysis   |

**Cell Filtering:**

These thresholds are applied during preprocessing. Defaults work well for most plant snRNA-seq datasets.

| Parameter                | Default | Description                          |
|--------------------------|---------|--------------------------------------|
| `--min_diem_debris_score`| `1`     | Minimum DIEM debris score            |
| `--min_unsplice_ratio`   | `0.1`   | Minimum unspliced/total ratio        |
| `--min_nCount_RNA`       | `400`   | Minimum UMIs per cell                |
| `--min_nFeature_RNA`     | `300`   | Minimum genes per cell               |
| `--max_mt`               | `10`    | Maximum mitochondrial % per cell     |
| `--max_cp`               | `15`    | Maximum chloroplast % per cell       |
| `--nHVG`                 | `3000`  | Number of highly variable genes      |
| `--min_ncell_expr`       | `5`     | Minimum cells expressing a gene      |
| `--remove_doublet`       | `true`  | Remove doublets                      |
| `--max_doublet_score`    | `0.4`   | Maximum doublet pANN score           |
| `--min_nClusterMarker`   | `5`     | Minimum cluster marker genes         |
| `--min_cells`            | `500`   | Minimum cells required after each filtering step; samples with fewer cells are skipped |

**Multi-species Demultiplexing and Filtering:**

These parameters only apply when running with multiple species (`--clean chi` is auto-selected for multi-species).

| Parameter                        | Default | Description                                    |
|----------------------------------|---------|------------------------------------------------|
| `--clean`                        | auto    | Cleaning method: `diem` (single-species) or `chi` (multi-species) |
| `--min_UMI_per_cell_barcode`     | `400`   | Minimum UMIs per barcode for chi demultiplexing |
| `--chisq_pvalues_max`            | `0.01`  | Maximum chi-squared p-value                    |
| `--ambient_rate_max`             | `0.5`   | Maximum ambient RNA rate                       |
| `--multiple_species_per_droplet` | `true`  | Allow multi-species droplets                   |

**Cleanup:**

| Parameter   | Default | Description                          |
|-------------|---------|--------------------------------------|
| `--cleanup` | `true`  | Cleanup large intermediate files     |

---

## Output Structure

```
<out>/
├── cellranger/          # CellRanger count outputs per sample
└── preprocess/          # Processed seurat objects, HTML summary report, and other results
```

seur_diem_*.rds: Seurat object after DIEM filtering with full metadata (including debris score, unsplicing ratio, doublet score, etc.) -> good for exploring different filtering cutoffs
seur_clean_*.rds: Seurat object after all filtering steps with full metadata -> good for downstream analysis like integration and differential analysis

---

## Pipeline Steps

1. **Reference preparation** — builds a CellRanger reference from FASTA + GTF (skipped if a pre-built reference exists)
2. **CellRanger count** — aligns FASTQs and generates cell × gene count matrices
3. **Preprocessing / Demultiplexing** — DIEM-based debris removal (single-species) or chi-squared demultiplexing (multi-species)
4. **Read counting** — separates uniquely/multi-mapped reads; counts reads in 2 kb windows
5. **Velocyto** — generates spliced/unspliced RNA velocity loom files
6. **Full preprocessing** — Seurat-based QC, normalization, clustering, and annotation using velocyto ratios
7. **Summary report** — HTML report aggregating QC metrics across samples

---

## Integration

Use `integrate` (entrypoint for `bin/integration.R`) on the per-sample Seurat objects produced by this pipeline (typically `preprocess/seur_clean_*.rds`).

### Usage 

By default, `integrate` automatically uses all `.rds/.RDS` files in the current directory (`.`) and produces all output files in the same directory.

```bash
mkdir integration && cd integration
ln -s /path/to/exp1/preprocess/seur_clean_*.rds .
ln -s /path/to/exp2/preprocess/seur_clean_*.rds .
module load scrnaseq/1.0 && integrate
```

Alternatively, you can specify the input Seurat objects explicitly:

```bash
module load scrnaseq/1.0 && integrate \
  --inputRds "/path/to/exp1/preprocess/seur_clean_sampleA.rds,/path/to/exp2/preprocess/seur_clean_sampleB.rds" \
  --resDir "/path/to/integration_results"
```

See full usage:

```bash
module load scrnaseq/1.0 && integrate --help
```

Full usage example:

```bash
module load scrnaseq/1.0 && integrate \
  --inputRds "/path/to/exp1/preprocess/seur_clean_sampleA.rds,/path/to/exp2/preprocess/seur_clean_sampleB.rds" \
  --resDir "/path/to/integration_results" \
  --integration_method "RPCA" \
  --nFeatures 3000 \
  --species "Arabidopsis thaliana" \
  --memory 16 \
  --ref_yaml "/path/to/scQC.yaml"
```

### Outputs

- `seuratObjs_<method>_integrated.rds`
- `UMAP_bySample.png`, `UMAP_byCluster.png`
- Violin plots of nCount_RNA, nFeature_RNA for each cluster
- cluster/sample composition CSV files and percentage barplots
- optional annotation plots when marker/reference entries are available for the selected species in `--ref_yaml`

> Note: Plot files named with `Celltype` are anchor-based annotation, while those named with `MarkerBasedCelltype` are marker-based annotation.

---

## Q&A

### What if I want to use different cutoff values for different samples?

You can run each sample in a separate nextflow run with different cutoff values. Then you can generate the [combined summary report](#how-can-i-build-one-summary-report-across-samples-from-different-experiments) to compare the results across samples.


### How can I build one summary report across samples from different experiments?

Render `bin/summary.Rmd` manually and pass combined paths from multiple runs.

Provide:
- `cellranger_dir`: absolute paths to per-sample CellRanger output directories
- `cellranger_sel`: sample IDs to include (must match names in those paths)
- `seur_dir`: absolute paths to preprocess/CELL_FILTERING work directories containing `seur_clean_*.rds` or `seur_diem_*.rds`
- `seur_sel`: sample IDs to include in Seurat/read-count sections
- `species`: full species name used in this combined report

Example:

```bash
cd /path/to/combined_summary_folder

module load scrnaseq/1.0

Rscript -e "rmarkdown::render(
  input = 'bin/summary.Rmd',
  params = list(
    title = 'Combined scRNA-seq summary',
    author = 'Qiaoshan Lin',
    cellranger_dir = c(
      '/absolute/path/to/exp1/cellranger/output/folder/sampleA',
      '/absolute/path/to/exp1/cellranger/output/folder/sampleB',
      '/absolute/path/to/exp2/cellranger/output/folder/sampleC',
      '/absolute/path/to/exp2/cellranger/output/folder/sampleD'
    ),
    cellranger_sel = c('sampleA', 'sampleC'),
    seur_dir = c(
      '/absolute/path/to/exp1/preprocess/output/folder',
      '/absolute/path/to/exp2/preprocess/output/folder'
    ),
    seur_sel = c('sampleA', 'sampleC'),
    species = 'Arabidopsis thaliana'
  ),
  output_dir = getwd(),
  output_file = 'summary_combined.html'
)"
```

Tip: use `cellranger_sel` / `seur_sel` to include only selected samples without moving any files.

---