# snRNAseq Processing Pipeline — User Guide

A Nextflow pipeline for processing single-nucleus RNA sequencing data from plants. It runs end-to-end: reference preparation → CellRanger count → read counting → filtering/demultiplexing → summary report.

---

## Requirements

- Nextflow ≥ 23.04.0
- SLURM cluster (Midway3, account `pi-gzy8899`)
- Modules: `cellranger/9.0.1`, `agat/1.6.1`, `scrnaseq/1.0`, `velocyto/0.17.17`, `samtools/1.22.1`, `deeptools/3.5.6`

---

## Quick Start

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
nextflow run jow30/ROMalleyLab-scRNAseq-nf --input samplesheet.csv --out exp1-nf-results {$otherParams} -resume
```

---

## Samplesheet Format

Provide a CSV with the following columns:

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

---

## Species

The following species have pre-configured references (genome, GTF, CellRanger index, Seurat reference atlas):

| Species                  | Markers available |
|--------------------------|-------------------|
| `Arabidopsis thaliana`   | Yes               |
| `Arabidopsis lyrata`     | No                |
| `Capsella rubella`       | No                |
| `Brassica oleracea`      | No                |

Pre-built multi-species cellranger references exist for:
- `Arabidopsis thaliana,Capsella rubella`
- `Arabidopsis thaliana,Arabidopsis lyrata,Capsella rubella,Brassica oleracea`

For any other species, supply `--genome` and `--gtf` manually (see below).

---

## Parameters

### Mandatory

| Parameter  | Description                  |
|------------|------------------------------|
| `--input`  | Path to samplesheet CSV      |
| `--out`    | Output directory path        |
| `--species`| Species name(s), comma-separated (default: `Arabidopsis thaliana`) |

### Custom Reference (if species not in the list above)

| Parameter  | Description                          |
|------------|--------------------------------------|
| `--genome` | Path to genome FASTA file            |
| `--gtf`    | Path to annotation file (`.gtf`, `.gff3`, or `.gff`) |
| `--ref_yaml` | Path to scQC YAML which contains information about organelles (default: `refs/scQC.yaml`) |

### Library QC

| Parameter                | Default | Description                                                  |
|--------------------------|---------|--------------------------------------------------------------|
| `--bin_size`             | `500`   | Window size to count reads for mapping efficiency analysis   |

### Cell Filtering

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

### Multi-species Demultiplexing and Filtering

These parameters only apply when running with multiple species (`--clean chi` is auto-selected for multi-species).

| Parameter                        | Default | Description                                    |
|----------------------------------|---------|------------------------------------------------|
| `--clean`                        | auto    | Cleaning method: `diem` (single-species) or `chi` (multi-species) |
| `--min_UMI_per_cell_barcode`     | `400`   | Minimum UMIs per barcode for chi demultiplexing |
| `--chisq_pvalues_max`            | `0.01`  | Maximum chi-squared p-value                    |
| `--ambient_rate_max`             | `0.5`   | Maximum ambient RNA rate                       |
| `--multiple_species_per_droplet` | `true`  | Allow multi-species droplets                   |

### Cleanup

| Parameter   | Default | Description                          |
|-------------|---------|--------------------------------------|
| `--cleanup` | `false` | Cleanup large intermediate files     |

---

## Output Structure

```
<out>/
├── cellranger/          # CellRanger count outputs per sample
├── preprocess/          # Seurat objects, QC plots, and tables
│   └── <species>/       # (multi-species only) per-species subdirectories
└── summary/             # HTML summary report
```

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

## Configuring `scQC.yaml`

The file `refs/scQC.yaml` tells the pipeline how to identify organellar and ribosomal genes for each species during QC filtering. You can also point it to a Seurat reference atlas for automated cell-type annotation. Edit this file when you add a new species or need to change organelle definitions.

Each top-level key is a species name (must match `--species` exactly). Under it, provide **one** of two strategies for mitochondrial and chloroplast genes:

### Strategy 1 — Regex patterns (when gene IDs encode organelle origin)

Use this when organellar genes can be identified by a naming prefix (e.g. *A. thaliana* Araport11 uses `ATMG*` / `ATCG*`).

```yaml
Arabidopsis thaliana:
    mitochondrial_pattern: "^ATMG"
    chloroplast_pattern: "^ATCG"
    ribosomal_genes: "/path/to/ribosomal_genes.txt"
    annotation_ref_seurat_obj: "/path/to/reference_atlas.RDS"   # optional
```

### Strategy 2 — Gene list files (when no naming convention exists)

Provide plain-text files with one gene ID per line.

```yaml
Brassica oleracea:
    mitochondrial_genes: "/path/to/mitochondrial_genes.txt"
    chloroplast_genes: "/path/to/chloroplast_genes.txt"
    ribosomal_genes: "/path/to/ribosomal_genes.txt"
```

### Field reference

| Field | Type | Description |
|-------|------|-------------|
| `mitochondrial_pattern` | regex | Grep pattern to match mitochondrial gene IDs (e.g. `"^ATMG"`) |
| `chloroplast_pattern` | regex | Grep pattern to match chloroplast gene IDs (e.g. `"^ATCG"`) |
| `mitochondrial_genes` | file path | One mitochondrial gene ID per line |
| `chloroplast_genes` | file path | One chloroplast gene ID per line |
| `ribosomal_genes` | file path | One ribosomal gene ID per line |
| `annotation_ref_seurat_obj` | file path | Seurat reference RDS with a `celltype` column for anchor-based label transfer (optional) |

> For each species, use **either** the pattern fields **or** the gene-list fields for mitochondria/chloroplast — not both. Ribosomal genes always require a file. The `annotation_ref_seurat_obj` field is optional; if omitted, cell-type annotation is skipped.

---

## Adding a New Species

Edit `nextflow.config` and add an entry under `params.species_map`:

```groovy
'My new species': [
    genome:      '/path/to/genome.fa',
    gtf:         '/path/to/annotation.gtf',
    cellranger:  "${projectDir}/refs/My_new_species/cellranger",  // or null to build on-the-fly
    seurat_ref:  '/path/to/reference_atlas.RDS',                  // or null
    markers:     '/path/to/markers.csv'                            // or null
]
```

If `cellranger` points to an existing directory it will be used directly; otherwise the pipeline builds it from `genome` and `gtf`.

---

## Help

```bash
nextflow run jow30/ROMalleyLab-scRNAseq-nf --help
```
