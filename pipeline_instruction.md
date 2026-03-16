# snRNAseq data processing pipeline instruction

## Purpose

This pipeline aims to process single-nucleus RNA sequencing (snRNAseq) data from plants using Nextflow. It orchestrates the entire workflow from reference preparation to final analysis and visualization.

## Environment

R and python scripts are available in /project/gzy8899/qiaoshan/scRNAseq/nextflow/bin

The following modules are required:
module load agat/1.5.0
module load cellranger/9.0.1
module load scrnaseq/1.0
module load velocyto/0.17.17
module load samtools/1.22.1
module load deeptools/3.5.6

## Input parameters

### Mandatory parameters
- `--input` : path to the input samplesheet
- `--out` : path to the output directory
- `--species` : species name(s), default is 'Arabidopsis thaliana'

### Optional parameters
#### reference parameters
- `--genome` : path(s) to the directory containing cellranger reference(s)
- `--gtf` : path(s) to the gtf file(s)
- `--ref_yaml` : path to the scQC.yaml file, default is 'refs/scQC.yaml'. Pass to `preprocess.R` commands and `demultiplex.R` commands.
#### filtering parameters (universal)
- `--min_diem_debris_score` : minimum debris score to remove debris, default is 1. This parameter is only applicable when `--clean` is set to 'diem'. Pass to `preprocess.R` commands.
- `--min_unsplice_ratio` : minimum Velocyto-derived unsplice ratio, default is 0.1. Pass to `preprocess.R` commands.
- `--min_nCount_RNA` : minimum UMI per cell, default is 400. Pass to `preprocess.R` commands.
- `--min_nFeature_RNA` : minimum genes per cell, default is 300. Pass to `preprocess.R` commands.
- `--max_mt` : maximum percent mitochondrial gene expression, default is 10. Pass to `preprocess.R` commands.
- `--max_cp` : maximum percent chloroplast gene expression, default is 15. Pass to `preprocess.R` commands.
- `--nHVG` : number of highly variable genes (HVGs), default is 3000. Pass to `preprocess.R` commands.
- `--min_ncell_expr` : minimum number of cells expressing a gene, default is 5. Pass to `preprocess.R` commands.
- `--remove_doublet` : whether to remove doublets, default is TRUE. Pass to `preprocess.R` commands.
- `--max_doublet_score` : maximum doubletFinder pANN, default is 0.4. Pass to `preprocess.R` commands.
- `--min_nClusterMarker` : minimum number of markers detected for a valid cluster, default is 5. Pass to `preprocess.R` commands.
#### filtering parameters (demultiplexing)
- `--clean` : method to filter for clean cells or nuclei, options include 'diem' and 'chi'. Default value is 'diem' for one species, 'chi' for more than one species. 
- `--min_UMI_per_cell_barcode` : minimum UMI per cell barcode, default is 400. This parameter is only applicable when `--clean` is set to 'chi'. Pass to `demultiplex.R` commands.
- `--chisq_pvalues_max` : maximum chi-squared p-values, default is 0.01. This parameter is only applicable when `--clean` is set to 'chi'. Pass to `demultiplex.R` commands.
- `--ambient_rate_max` : maximum ambient rate, default is 0.5. This parameter is only applicable when `--clean` is set to 'chi'. Pass to `demultiplex.R` commands.
- `--multiple_species_per_droplet` : whether there are multiple species per droplet, default is TRUE. This parameter is only applicable when `--clean` is set to 'chi'. Pass to `demultiplex.R` commands.

Default values are set in `nextflow.config`.

Parameters will be passed to the corresponding commands in the following steps.

## Processing Steps

### 1) Verify parameters

Verify the parameters to make sure they are correct.

- `--input` : The samplesheet should be in the format of CSV. Make sure the samplesheet path is available. The samplesheet should have the following columns: sample, fastq_dir, fastq_R1, fastq_R2. The sample names should only contain letters, numbers, dashes and underscores. The paths of fastq files in the samplesheet should be correct.
- `--out` : Make sure the output directory contains only letters, numbers, dashes and underscores; if not, throw out an error message to ask the user to rename the output directory. If the output directory already exists, throw out a warning and overwrite the existing files
- `--species` : species name(s), default is 'Arabidopsis thaliana'. Other available species are 'Arabidopsis lyrata', 'Capsella rubella'. If user provide a species name that is not available, throw out an error message to ask the user to either provide an available species name or add a new species name and the corresponding `--genome` and `--gtf` parameters in the nextflow.config file. If user provides more than one species name, split by commas and check if the species names are available. If any species name is not available, throw out an error message to ask the user to either provide available species names or add new species names and the corresponding `--genome` and `--gtf` parameters in the nextflow.config file.
- `--genome` : Make sure the provided path(s) to the file(s) exists
- `--gtf` : Make sure the provided path(s) to the file(s) exists
- `--ref_yaml` : Make sure the provided path to the file exists
- `--min_diem_debris_score` : if user set this parameter while `--clean` is set to 'chi', throw out a warning message to inform the user that this parameter is ignored when `--clean` is set to 'chi'.
- `--min_unsplice_ratio` : Make sure it ranges from 0 to 1
- `--min_nCount_RNA` : Make sure it's a positive integer
- `--min_nFeature_RNA` : Make sure it's a positive integer
- `--max_mt` : Make sure it ranges from 0 to 100
- `--max_cp` : Make sure it ranges from 0 to 100
- `--nHVG` : Make sure it's a positive integer
- `--min_ncell_expr` : Make sure it's a positive integer
- `--remove_doublet` : Make sure it's a boolean value
- `--max_doublet_score` : Make sure it ranges from 0 to 1
- `--min_nClusterMarker` : Make sure it's a positive integer
- `--clean` : Throw out an error message if the user provides a value other than 'diem' or 'chi'. If the user provides one species name and use 'chi' for `--clean`, throw out a warning message to inform the user that only 'diem' is available for processing a single species. Then use 'diem' for `--clean`.
- `--min_UMI_per_cell_barcode` : Make sure it's a positive integer. Throw out a warning message to the user that this parameter will be ignored when `--clean` is set to 'diem' if the user sets this parameter while `--clean` is set to 'diem'. 
- `--chisq_pvalues_max` : Make sure it ranges from 0 to 1. Throw out a warning message to the user that this parameter will be ignored when `--clean` is set to 'diem' if the user sets this parameter while `--clean` is set to 'diem'.
- `--ambient_rate_max` : Make sure it ranges from 0 to 1. Throw out a warning message to the user that this parameter will be ignored when `--clean` is set to 'diem' if the user sets this parameter while `--clean` is set to 'diem'.
- `--multiple_species_per_droplet` : Make sure it's a boolean value. Throw out a warning message to the user that this parameter will be ignored when `--clean` is set to 'diem' if the user sets this parameter while `--clean` is set to 'diem'.

### 2) Make the output directory for the project. 

Output directory should be specified by the user in the `--out` parameter.

### 3) Prepare the reference files

If user-provided species are available, the pipeline will use the corresponding reference files in `refs` directory. If not, the pipeline will use the provided `--genome` and `--gtf` parameters to create a cellranger reference folder named as the user-provided species name in `refs` directory.

Example 1: Create a cellranger reference for a single species
```
genome=/project/gzy8899/references/Arabidopsis_thaliana/Phytozome/PhytozomeV12/early_release/Athaliana_447_Araport11/assembly/Athaliana_447_TAIR10.fa
gff=/project/gzy8899/references/Arabidopsis_thaliana/Phytozome/PhytozomeV13/Athaliana/Araport11/annotation/Athaliana_447_Araport11.gene.gff3

out_dir=/project/gzy8899/qiaoshan/scRNAseq/nextflow/refs/Arabidopsis_thaliana
genome_name=Arabidopsis_thaliana

module load agat/1.5.0

basename=`basename $gff .gff3`

agat_convert_sp_gff2gtf.pl --gff $gff -o $out_dir/$basename.gtf

module load cellranger/9.0.1

cellranger mkgtf \
    $out_dir/$basename.gtf \
    $out_dir/$basename.flt.gtf \
    --attribute=gene_biotype:protein_coding \
    --attribute=gene_biotype:lncRNA \
    --attribute=gene_biotype:miRNA \
    --attribute=gene_biotype:snoRNA \
    --attribute=gene_biotype:snRNA \
    --attribute=transcript_biotype:protein_coding \
    --attribute=transcript_biotype:lncRNA \
    --attribute=transcript_biotype:miRNA \
    --attribute=transcript_biotype:snoRNA \
    --attribute=transcript_biotype:snRNA

wc -l $out_dir/$basename.flt.gtf $out_dir/$basename.gtf
mv $out_dir/$basename.flt.gtf $out_dir/$basename.gtf

export MRO_DISK_SPACE_CHECK=disable

cellranger mkref \
    --genome=$genome_name \
    --fasta=$genome \
    --genes=$out_dir/$basename.gtf \
    --nthreads 16 \
    --memgb 24 \
    --localcores 16 \
    --localmem 24 \
    --output-dir $out_dir/cellranger
```

Example 2: Create a cellranger reference for multiple species
```
out_dir=/project/gzy8899/qiaoshan/scRNAseq/nextflow/refs/Athaliana_Crubella

module load agat/1.5.0
module load cellranger/9.0.1

gff=(/project/gzy8899/references/Arabidopsis_thaliana/Phytozome/PhytozomeV13/Athaliana/Araport11/annotation/Athaliana_447_Araport11.gene.gff3 /project/gzy8899/references/Capsella_rubella/Phytozome/Crubella_474_v1.1.gene.gff3)

for gff in ${gff[@]}; do
    basename=`basename $gff .gff3`
    agat_convert_sp_gff2gtf.pl --gff $gff -o $out_dir/$basename.gtf
    cellranger mkgtf \
        $out_dir/$basename.gtf \
        $out_dir/$basename.flt.gtf \
        --attribute=gene_biotype:protein_coding \
        --attribute=gene_biotype:lncRNA \
        --attribute=gene_biotype:miRNA \
        --attribute=gene_biotype:snoRNA \
        --attribute=gene_biotype:snRNA \
        --attribute=transcript_biotype:protein_coding \
        --attribute=transcript_biotype:lncRNA \
        --attribute=transcript_biotype:miRNA \
        --attribute=transcript_biotype:snoRNA \
        --attribute=transcript_biotype:snRNA

    wc -l $out_dir/$basename.flt.gtf $out_dir/$basename.gtf
    mv $out_dir/$basename.flt.gtf $out_dir/$basename.gtf
done

export MRO_DISK_SPACE_CHECK=disable

cellranger mkref \
    --genome=Arabidopsis_thaliana --fasta=/project/gzy8899/references/Arabidopsis_thaliana/Phytozome/PhytozomeV12/early_release/Athaliana_447_Araport11/assembly/Athaliana_447_TAIR10.fa --genes=$out_dir/Athaliana_447_Araport11.gene.gtf \
    --genome=Capsella_rubella --fasta=/project/gzy8899/references/Capsella_rubella/Phytozome/Crubella_474_v1.fa --genes=$out_dir/Crubella_474_v1.1.gene.gtf \
    --nthreads 16 \
    --memgb 64 \
    --localcores 16 \
    --localmem 64 \
    --output-dir $out_dir/cellranger
```

### 4) Run `cellranger_count`

Create a folder called `cellranger` in the output directory. Allocate adequate computational resources for the job. In principle, each sample should be processed in a separate job. 

**Computational resources allocation criteria:**

- core: 16
- memory: 8G (increase if the job failed with out-of-memory error)
- time: 12h (increase if the job failed with time limit error)

**How to set the parameters:**

- `--id` : the sample name specified by user in the input samplesheet to name the cellranger output directory
- `--fastqs` : the path to the directory containing fastq files
- `--sample` : the prefix that identify the sample in the fastq files, which is the string before the first underscore in the file name
- `--transcriptome` : the path to the directory containing cellranger reference, when multiple species are present, use the combined reference according to the instructions in the nextflow.config file
- `--localcores` : the number of cores to use
- `--localmem` : the amount of memory to use

Example 1: one library/lane per sample
```
module load cellranger/9.0.1

TRANSCRIPTOME=/project/gzy8899/qiaoshan/scRNAseq/references/Arabidopsis_thaliana/cellranger
raw_data_dir=/project/gzy8899/data/202509_Uchicago_DAP_10X_RNA

id=(Athaliana_5d_seedling_XC Athaliana_5d_seedling_LG)
prefix=(RO-LG-1s-200 RO-LG-1s-201)

export MRO_DISK_SPACE_CHECK=disable

cellranger count --create-bam true --id=${id[$SLURM_ARRAY_TASK_ID]} --fastqs=$raw_data_dir --sample=${prefix[$SLURM_ARRAY_TASK_ID]} --transcriptome=$TRANSCRIPTOME --localcores=16 --localmem=24
```
Example 2: multiple libraries/lanes per sample
```
module load cellranger/9.0.1

TRANSCRIPTOME=/project/gzy8899/qiaoshan/scRNAseq/references/Arabidopsis_thaliana/cellranger
raw_data_dir=/project/gzy8899/qiaoshan/scRNAseq/experiments/At_5d_seedling_2025/raw_data

id=(At_5d_seedling_LG At_5d_seedling_XC)
sample_lib1=(RO-LG-1s-01 RO-LG-1s-02)
sample_lib2=(XY-186 XY-187)
sample_lib3=(RO-GY-2s-pl1-10x-At-5d-seedling-LG RO-GY-2s-pl1-10x-At-5d-seedling-xc)

export MRO_DISK_SPACE_CHECK=disable

cellranger count --create-bam true --id=${id[$SLURM_ARRAY_TASK_ID]} --fastqs=$raw_data_dir --sample=${sample_lib1[$SLURM_ARRAY_TASK_ID]},${sample_lib2[$SLURM_ARRAY_TASK_ID]},${sample_lib3[$SLURM_ARRAY_TASK_ID]} --transcriptome=$TRANSCRIPTOME --localcores=16 --localmem=20
``` 
Example 3: multi-species sample
```
module load cellranger/9.0.1

TRANSCRIPTOME=/project/gzy8899/qiaoshan/scRNAseq/references/brassica-mix/cellranger
raw_data_dir=/project/gzy8899/qiaoshan/scRNAseq/experiments/control_exp_brassica-mix_shoot

id=(brassica_mix_rep1 brassica_mix_rep2)
prefix=(rep1 rep2)

export MRO_DISK_SPACE_CHECK=disable

cellranger count --create-bam true --id=${id[$SLURM_ARRAY_TASK_ID]} --fastqs=$raw_data_dir --sample=${prefix[$SLURM_ARRAY_TASK_ID]} --transcriptome=$TRANSCRIPTOME --localcores=16 --localmem=24
``` 

Run jobs within the `cellranger` folder. When jobs are completed, each sample will have a folder named as the sample name specified by user that contains the cellranger output files.

### 5) Filtering for clean cells/nuclei (and optional demultiplexing)

Create a folder called `preprocess` in the output directory. 

#### Condition 1: only one species is present

Run `preprocess.R` to filter for clean cells/nuclei using cellranger output. 

Example:
```
module load scrnaseq/1.0

cellranger_outDir=/project/gzy8899/qiaoshan/scRNAseq/control_exp/cellranger
id=(SRR33560988 SRR33560965 SRR33560964) # correspond to sample names specified by user in input samplesheet

Rscript bin/preprocess.R \
        -i $cellranger_outDir/${id[$SLURM_ARRAY_TASK_ID]}/outs \
        -s ${id[$SLURM_ARRAY_TASK_ID]} \
        -S 'Arabidopsis thaliana' \
        -o 'results/preprocess' \
        --min_diem_debris_score 1 # keep only cells/nuclei with debris score <= 1
```

#### Condition 2: multiple species are present

Run `demultiplex.R` to demultiplex and filter for clean cells/nuclei.

Example: 
```
module load scrnaseq/1.0

input=(/project/gzy8899/qiaoshan/scRNAseq/experiments/At_Cr_5d_seedling_2025/At_Cr_5d_seedling_XC/outs /project/gzy8899/qiaoshan/scRNAseq/experiments/At_Cr_5d_seedling_2025/At_Cr_5d_seedling_LG/outs)
sample=(At_Cr_5d_seedling_XC At_Cr_5d_seedling_LG)
outdir=preprocess
species=("Arabidopsis thaliana" "Capsella rubella")

species_csv=$(IFS=,; echo "${species[*]}")
Rscript bin/demultiplex.R -i ${input[$SLURM_ARRAY_TASK_ID]} -o $outdir -s ${sample[$SLURM_ARRAY_TASK_ID]} --species "$species_csv" \
        --ref_yaml refs/scQC.yaml \
        --min_UMI_per_cell_barcode 150 \
        --chisq_pvalues_max 0.01 \
        --ambient_rate_max 0.5 \
        --multiple_species_per_droplet TRUE \
        --demux_matrix TRUE
```

The output of `demultiplex.R` are subfolders named as the species name inside the preprocess folder. In each species folder, there are outputs for each sample as indicated by the sample name. The subfolders named as the sample name contains demultiplexed raw_feature_bc_matrix, including barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz, which will be used in the downstream analysis.

Example:
```
preprocess/
├── Arabidopsis_thaliana/
│   ├── At_Cr_5d_seedling_LG/
│   │── seurat_diem_barcodes_At_Cr_5d_seedling_LG.tsv.gz
│   │── seurat_objs_At_Cr_5d_seedling_LG.rds
│   │── summary_dims_At_Cr_5d_seedling_LG.rds
│   │── summary_opts_At_Cr_5d_seedling_LG.rds
│   │── summary_plts_At_Cr_5d_seedling_LG.rds
│   │── summary_tbs_At_Cr_5d_seedling_LG.rds
│   │── At_Cr_5d_seedling_XC/
│   ├── seurat_diem_barcodes_At_Cr_5d_seedling_XC.tsv.gz
│   ├── seurat_objs_At_Cr_5d_seedling_XC.rds
│   ├── summary_dims_At_Cr_5d_seedling_XC.rds
│   ├── summary_opts_At_Cr_5d_seedling_XC.rds
│   ├── summary_plts_At_Cr_5d_seedling_XC.rds
│   └── summary_tbs_At_Cr_5d_seedling_XC.rds
└── Capsella_rubella/
    ├── At_Cr_5d_seedling_LG/
    │── seurat_diem_barcodes_At_Cr_5d_seedling_LG.tsv.gz
    │── seurat_objs_At_Cr_5d_seedling_LG.rds
    │── summary_dims_At_Cr_5d_seedling_LG.rds
    │── summary_opts_At_Cr_5d_seedling_LG.rds
    │── summary_plts_At_Cr_5d_seedling_LG.rds
    │── summary_tbs_At_Cr_5d_seedling_LG.rds
    │── At_Cr_5d_seedling_XC/
    │── seurat_diem_barcodes_At_Cr_5d_seedling_XC.tsv.gz
    │── seurat_objs_At_Cr_5d_seedling_XC.rds
    │── summary_dims_At_Cr_5d_seedling_XC.rds
    │── summary_opts_At_Cr_5d_seedling_XC.rds
    │── summary_plts_At_Cr_5d_seedling_XC.rds
    └── summary_tbs_At_Cr_5d_seedling_XC.rds
```

Then demultiplex bam files for each species.

Example:
```

input=(/project/gzy8899/qiaoshan/scRNAseq/experiments/At_Cr_5d_seedling_2025/At_Cr_5d_seedling_XC/outs /project/gzy8899/qiaoshan/scRNAseq/experiments/At_Cr_5d_seedling_2025/At_Cr_5d_seedling_LG/outs)
species=("Arabidopsis thaliana" "Capsella rubella")
sample=("At_Cr_5d_seedling_XC" "At_Cr_5d_seedling_LG")
outdir=preprocess

set -euo pipefail
module load samtools

INPUT_BAM=${input[$SLURM_ARRAY_TASK_ID]}/possorted_genome_bam.bam
THREADS="${SLURM_CPUS_PER_TASK:-1}"

for S in "${species[@]}"; do
  PREFIX_REGEX="${S// /_}"
  OUTPUT_BAM=$outdir/${PREFIX_REGEX}/${sample[$SLURM_ARRAY_TASK_ID]}/possorted_genome_bam.bam

  echo "Start demultiplexing bam file for ${S}"
  samtools view -h -@ "${THREADS}" "${INPUT_BAM}" |
  awk -v re="(${PREFIX_REGEX})_+" '
    BEGIN { OFS="\t" }
    /^@/ {
      gsub(re, "", $0)
      print
      next
    }
    {
        if($3 ~ "^" re) {
                for (i = 1; i <= NF; i++) { gsub(re, "", $i) }
                print
        }
    }
  ' | samtools view -@ "${THREADS}" -b -o "${OUTPUT_BAM}" -
  samtools index -@ "${THREADS}" "${OUTPUT_BAM}"
  echo "Finish demultiplexing for ${S}."

done
```

If `--clean` is set to 'chi', go to next step.
If `--clean` is set to 'diem', run `preprocess.R` using the demultiplexed raw_feature_bc_matrix to filter for clean cells/nuclei. Example:
```
module load scrnaseq/1.0

sp=(Arabidopsis_thaliana Arabidopsis_thaliana Capsella_rubella Capsella_rubella)
id=(At_Cr_5d_seedling_LG At_Cr_5d_seedling_XC At_Cr_5d_seedling_LG At_Cr_5d_seedling_XC)
cellranger_outDir=/project/gzy8899/qiaoshan/scRNAseq/experiments/At_Cr_5d_seedling_2025/demultiplex/${sp[$SLURM_ARRAY_TASK_ID]}/${id[$SLURM_ARRAY_TASK_ID]} # note that this must be absolute path
sp_name=("Arabidopsis thaliana" "Arabidopsis thaliana" "Capsella rubella" "Capsella rubella")

Rscript /project/gzy8899/qiaoshan/scRNAseq/scripts/preprocess.R \
        -i $cellranger_outDir \
        -s ${id[$SLURM_ARRAY_TASK_ID]} \
        -S "${sp_name[$SLURM_ARRAY_TASK_ID]}" \
        -o /project/gzy8899/qiaoshan/scRNAseq/experiments/At_Cr_5d_seedling_2025/demultiplex/${sp[$SLURM_ARRAY_TASK_ID]} \
        --min_diem_debris_score 1 # keep only cells/nuclei with debris score <= 1
```

### 6) Count reads in windows of 2kb across the genome

Once demultiplexed bam files are ready, separate mapped, multi-mapped, and uniquely mapped reads, and count the number of reads in windows of 2kb across the genome for each species and each sample.

Example for separating mapped, multi-mapped, and uniquely mapped reads:
```
module load samtools

id=(At_Cr_5d_seedling_LG At_Cr_5d_seedling_XC)

samtools view -b -F 4 -@ 16 ${id[$SLURM_ARRAY_TASK_ID]}/outs/possorted_genome_bam.bam > ${id[$SLURM_ARRAY_TASK_ID]}.mapped.bam
samtools view -b -q 255 -@ 16 ${id[$SLURM_ARRAY_TASK_ID]}.mapped.bam -o ${id[$SLURM_ARRAY_TASK_ID]}.uniq.bam -U ${id[$SLURM_ARRAY_TASK_ID]}.multi.bam
```

Example for counting reads in windows of 2kb across the genome:
```
module load samtools

# ---- Input/Output ----
BAM=(At_Cr_5d_seedling_LG.mapped.bam At_Cr_5d_seedling_XC.mapped.bam At_Cr_5d_seedling_LG.multi.bam At_Cr_5d_seedling_XC.multi.bam At_Cr_5d_seedling_LG.uniq.bam At_Cr_5d_seedling_XC.uniq.bam)
OUT=(At_Cr_5d_seedling_LG.mapped.tsv At_Cr_5d_seedling_XC.mapped.tsv At_Cr_5d_seedling_LG.multi.tsv At_Cr_5d_seedling_XC.multi.tsv At_Cr_5d_seedling_LG.uniq.tsv At_Cr_5d_seedling_XC.uniq.tsv)

samtools index ${BAM[$SLURM_ARRAY_TASK_ID]}

module load deeptools/3.5.6

# ---- Parameters ----
BIN_SIZE=2000           # Window size in bases (2 kb)
THREADS=8               # Match #SBATCH -c above

# Step 1: Generate bedgraph (chr, start, end, count)
BEDGRAPH="${OUT[$SLURM_ARRAY_TASK_ID]%.tsv}.bedgraph"

bamCoverage \
    -b "${BAM[$SLURM_ARRAY_TASK_ID]}" \
    -o "$BEDGRAPH" \
    --binSize "$BIN_SIZE" \
    --normalizeUsing None \
    --outFileFormat bedgraph \
    -p "$THREADS"

# Step 2: Add header and ratio column (count / total_counts)
awk 'BEGIN {OFS="\t"}
     NR==FNR {total += $4; next}
     FNR==1  {print "chr", "start", "end", "count", "ratio"}
             {print $1, $2, $3, $4, $4/total}
' "$BEDGRAPH" "$BEDGRAPH" > "${OUT[$SLURM_ARRAY_TASK_ID]}"

echo "Done. Total bins: $(tail -n +2 "${OUT[$SLURM_ARRAY_TASK_ID]}" | wc -l)"

# Clean up intermediate bedgraph
rm -f "$BEDGRAPH"
```

If only one species, the path to save the output files should be just `preprocess/`; if multiple species, the path should be like `preprocess/Arabidopsis_thaliana/At_Cr_5d_seedling_LG/`.

### 7) Run velocyto

Example:
```
module load velocyto/0.17.17
module load samtools/1.22.1

cellranger_resDir=./cellranger
seurat_resDir=./preprocess

id=(SRR33560988 SRR33560965 SRR33560964)

gtf=/project/gzy8899/qiaoshan/scRNAseq/references/Arabidopsis_thaliana/Athaliana_447_Araport11.gene.gtf

# velocyto run10x -@ $SLURM_CPUS_PER_TASK $cellranger_resDir/${id[$SLURM_ARRAY_TASK_ID]} $gtf

velocyto run -@ $SLURM_CPUS_PER_TASK \
 -b $seurat_resDir/seur_diem_barcodes_${id[$SLURM_ARRAY_TASK_ID]}.tsv.gz \
 -o $cellranger_resDir/${id[$SLURM_ARRAY_TASK_ID]} \
 $cellranger_resDir/${id[$SLURM_ARRAY_TASK_ID]}/outs/possorted_genome_bam.bam \
 $gtf

mv $cellranger_resDir/${id[$SLURM_ARRAY_TASK_ID]}/possorted_genome_bam_*.loom $cellranger_resDir/${id[$SLURM_ARRAY_TASK_ID]}/possorted_genome_bam.loom
python bin/unsplice_ratio.py $cellranger_resDir/${id[$SLURM_ARRAY_TASK_ID]}/possorted_genome_bam.loom $seurat_resDir/${id[$SLURM_ARRAY_TASK_ID]}.txt # this txt file will be used in the -v option of preprocess.R
```

Note that for multi-species samples, the `cellranger_resDir` should be the demultiplexed cellranger output directory, which is like `cellranger/Arabidopsis_thaliana/At_Cr_5d_seedling_LG/`; the `seurat_resDir` should be the demultiplexed seurat output directory, which is like `preprocess/Arabidopsis_thaliana/`.

### 8) Run `preprocess.slurm` with the `-v` option to include velocyto results

```
module load scrnaseq/1.0

cellranger_outDir=/project/gzy8899/qiaoshan/scRNAseq/control_exp/cellranger
id=(SRR33560988 SRR33560965 SRR33560964)

Rscript /project/gzy8899/qiaoshan/scRNAseq/scripts/preprocess.R \
        -i $cellranger_outDir/${id[$SLURM_ARRAY_TASK_ID]}/outs \
        -s ${id[$SLURM_ARRAY_TASK_ID]} \
        -S 'Arabidopsis thaliana' \
        -o '/project/gzy8899/qiaoshan/scRNAseq/control_exp/seurat' \
        -v /project/gzy8899/qiaoshan/scRNAseq/control_exp/${id[$SLURM_ARRAY_TASK_ID]}.txt \
        --ref_yaml /project/gzy8899/qiaoshan/scRNAseq/references/scQC.yaml \
        --min_diem_debris_score 1 \
        --min_unsplice_ratio 0.1 \
        --min_nCount_RNA 400 \
        --min_nFeature_RNA 300 \
        --max_mt 10 \
        --max_cp 15 \
        --nHVG 3000 \
        --min_ncell_expr 5 \
        --remove_doublet TRUE \
        --max_doublet_score 0.4 \
        --min_nClusterMarker 5 \
        --threads 8
```

### 9) Run `summary.Rmd`

Example for single-species samples:
```
module load scrnaseq/1.0

sp=Arabidopsis_thaliana

Rscript -e "rmarkdown::render(
  input = '/project/gzy8899/qiaoshan/scRNAseq/scripts/summary.Rmd',
  params = list(
    title = 'scRNA-seq summary report',
    author = 'Qiaoshan Lin',
    cellranger_dir = c(
      '/project/gzy8899/qiaoshan/scRNAseq/control_exp/cellranger/SRR33560988/outs',
      '/project/gzy8899/qiaoshan/scRNAseq/control_exp/cellranger/SRR33560965/outs',
      '/project/gzy8899/qiaoshan/scRNAseq/control_exp/cellranger/SRR33560964/outs'
    ),
    cellranger_sel = c(
      'SRR33560988',
      'SRR33560965',
      'SRR33560964'
    ),
    seur_dir = c(
      '/project/gzy8899/qiaoshan/scRNAseq/control_exp/seurat/$sp'
    ),
    seur_sel = c(
      'SRR33560988',
      'SRR33560965',
      'SRR33560964'
    ),
    output_path = '/project/gzy8899/qiaoshan/scRNAseq/control_exp/summary.html'
  ),
  output_file = '/project/gzy8899/qiaoshan/scRNAseq/control_exp/summary.html'
)"
```

Example for multi-species samples:
```
module load scrnaseq/1.0

sp=(Arabidopsis_thaliana Capsella_rubella)

Rscript -e "rmarkdown::render(
  input = '/project/gzy8899/qiaoshan/scRNAseq/scripts/summary.Rmd',
  params = list(
    title = 'scRNA-seq summary report',
    author = 'Qiaoshan Lin',
    cellranger_dir = c(
      '/project/gzy8899/qiaoshan/scRNAseq/experiments/At_Cr_5d_seedling_2025/At_Cr_5d_seedling_XC',
      '/project/gzy8899/qiaoshan/scRNAseq/experiments/At_Cr_5d_seedling_2025/At_Cr_5d_seedling_LG'
    ),
    cellranger_sel = c(
      'At_Cr_5d_seedling_XC',
      'At_Cr_5d_seedling_LG'
    ),
    seur_dir = c(
      '/project/gzy8899/qiaoshan/scRNAseq/experiments/At_Cr_5d_seedling_2025/demultiplex/${sp[$SLURM_ARRAY_TASK_ID]}'
    ),
    seur_sel = c(
      'At_Cr_5d_seedling_XC',
      'At_Cr_5d_seedling_LG'
    )
  ),
  output_file = '/project/gzy8899/qiaoshan/scRNAseq/experiments/At_Cr_5d_seedling_2025/summary_${sp[$SLURM_ARRAY_TASK_ID]}.html'
)"
```

### 10) Integration

Example:
```
module load scrnaseq/1.0

inputRds="output_dir/preprocess/sample1.rds,output_dir/preprocess/sample2.rds,output_dir/preprocess/sample3.rds"
resDir="output_dir/integration"

Rscript /project/gzy8899/qiaoshan/scRNAseq/nextflow/bin/integration.R \
        -i $inputRds \
        -o $resDir \
        --nFeatures 3000
```

## Output folder structure

For single-species samples:
```
output_dir/
├── cellranger/
│   ├── sample_1/
│   │   ├── outs/
│   │   │   ├── filtered_feature_bc_matrix/
│   │   │   ├── raw_feature_bc_matrix/
│   │   │   ├── possorted_genome_bam.bam
│   │   │   ├── possorted_genome_bam.loom
│   │   │   └── ...
│   │   └── ...
│   └── sample_2/
│       └── ...
├── preprocess/
│   ├── summary.html
│   ├── sample_1/
│   │   ├── Seurat_object.rds
│   │   └── ...
│   └── sample_2/
│       └── ...
├── integration/
│   ├── Seurat_object.rds
│   └── ...
└── ...
```

For multi-species samples:
```
output_dir/
├── cellranger/
│   ├── sample_1/
│   │   ├── outs/
│   │   │   ├── filtered_feature_bc_matrix/
│   │   │   ├── raw_feature_bc_matrix/
│   │   │   ├── possorted_genome_bam.bam
│   │   │   └── ...
│   │   └── ...
│   └── sample_2/
│       └── ...
├── preprocess/
│   ├── species_1/
│   │   ├── summary.html
│   │   ├── sample_1/
│   │   │   ├── raw_feature_bc_matrix/
│   │   │   ├── possorted_genome_bam.bam
│   │   │   ├── possorted_genome_bam.loom
│   │   │   └── ...
│   │   ├── sce_sample_1.rds
│   │   ├── seur_objs_sample_1.rds
│   │   ├── summary_opts_sample_1.rds
│   │   ├── summary_dims_sample_1.rds
│   │   ├── summary_plts_sample_1.rds
│   │   ├── summary_tbs_sample_1.rds
│   │   ├── seur_diem_barcodes_sample_1.tsv.gz
│   │   └── ...
│   │   └── sample_2/
│   │   │   ├── raw_feature_bc_matrix/
│   │   │   ├── possorted_genome_bam.bam
│   │   │   ├── possorted_genome_bam.loom
│   │   │   └── ...
│   │   ├── sce_sample_2.rds
│   │   ├── seur_objs_sample_2.rds
│   │   ├── summary_opts_sample_2.rds
│   │   ├── summary_dims_sample_2.rds
│   │   ├── summary_plts_sample_2.rds
│   │   ├── summary_tbs_sample_2.rds
│   │   ├── seur_diem_barcodes_sample_2.tsv.gz
│   │   └── ...
│   └── species_2/
│       ├── summary.html
│       ├── sample_1/
│       │   ├── raw_feature_bc_matrix/
│       │   └── ...
│       ├── ...
│       ├── sample_2/
│       │   └── ...
│       └── ...
├── integration/
│   ├── Seurat_object.rds
│   └── ...
└── ...
```

## Other notes

In R scripts or the example commands given above, sometimes absolute paths are used, but in the nextflow pipeline, relative paths should be used. Please modify the R scripts to use relative paths. However, some paths need to be absolute, such as the preprocess.R -i option, which is the absolute path to the cellranger output directory. Therefore, for the data-related paths, please use absolute paths.

The example commands given above come from many different experiments, so the paths are not consistent. In the nextflow pipeline, the paths should be consistent. Please aim to build an output directory structure as shown in the example output folder structure above.

Reference files like genome sequences, GTF files, cellranger reference files, etc. should be placed in the refs directory. The nextflow pipeline should be able to find these files using relative paths. When user creates a new reference species by providing the genome sequence and GTF file, the pipeline should be able to create the cellranger reference files automatically and place all the reference files in the refs directory, so that it can be used by the pipeline in the future.

