/*
 * Preprocessing module: DIEM filtering and demultiplexing
 */

process PREPROCESS_INITIAL {
    label 'process_medium'
    label 'scrnaseq'
    publishDir "${publish_dir}", mode: 'copy'

    input:
    tuple val(sample_id), val(cellranger_outs_path)
    val species_name
    val publish_dir

    output:
    tuple val(sample_id), path("seur_diem_barcodes_${sample_id}.tsv.gz"), emit: barcodes
    tuple val(sample_id), path("seur_objs_${sample_id}.rds"),            emit: seur_objs
    tuple val(sample_id), path("summary_dims_${sample_id}.rds"),         emit: summary_dims
    tuple val(sample_id), path("summary_tbs_${sample_id}.rds"),          emit: summary_tbs
    tuple val(sample_id), path("summary_plts_${sample_id}.rds"),         emit: summary_plts
    tuple val(sample_id), path("summary_opts_${sample_id}.rds"),         emit: summary_opts
    tuple val(sample_id), path("sce_${sample_id}.rds"),                  emit: sce

    script:
    """
    Rscript ${projectDir}/bin/preprocess_initial.R \\
        -i ${cellranger_outs_path} \\
        -s ${sample_id} \\
        -S '${species_name}' \\
        -o . \\
        --ref_yaml ${params.ref_yaml} \\
        --min_diem_debris_score ${params.min_diem_debris_score} \\
        --threads ${task.cpus}
    """
}

process DEMULTIPLEX {
    label 'process_medium'
    label 'scrnaseq'
    publishDir "${params.out}/preprocess", mode: 'copy'

    input:
    tuple val(sample_id), val(cellranger_outs_path)
    val species_csv
    val species_list

    output:
    tuple val(sample_id), path("**/*.rds"),                                        emit: rds_files
    tuple val(sample_id), path("**/seur_diem_barcodes_${sample_id}.tsv.gz"),       emit: barcodes
    tuple val(sample_id), path("**/raw_feature_bc_matrix"),                        emit: demux_matrices

    script:
    """
    Rscript ${projectDir}/bin/demultiplex.R \\
        -i ${cellranger_outs_path} \\
        -o . \\
        -s ${sample_id} \\
        --species "${species_csv}" \\
        --ref_yaml ${params.ref_yaml} \\
        --min_UMI_per_cell_barcode ${params.min_UMI_per_cell_barcode} \\
        --chisq_pvalues_max ${params.chisq_pvalues_max} \\
        --ambient_rate_max ${params.ambient_rate_max} \\
        --multiple_species_per_droplet ${params.multiple_species_per_droplet ? 'TRUE' : 'FALSE'} \\
        --demux_matrix TRUE
    """
}

process DEMUX_BAM {
    label 'process_medium'
    label 'samtools'
    publishDir "${params.out}/preprocess", mode: 'copy'

    input:
    tuple val(sample_id), val(cellranger_outs_path)
    val species_list

    output:
    tuple val(sample_id), path("**/possorted_genome_bam.bam"),     emit: demux_bams
    tuple val(sample_id), path("**/possorted_genome_bam.bam.bai"), emit: demux_bam_indices

    script:
    def species_cmds = species_list.collect { sp ->
        def prefix_regex = sp.replaceAll(' ', '_')
        """echo "Start demultiplexing bam file for ${sp}"
mkdir -p ${prefix_regex}/${sample_id}
INPUT_BAM=${cellranger_outs_path}/possorted_genome_bam.bam
OUTPUT_BAM=${prefix_regex}/${sample_id}/possorted_genome_bam.bam
samtools view -h -@ ${task.cpus} "\${INPUT_BAM}" | \\
awk -v re="(${prefix_regex})_+" '
    BEGIN { OFS="\\t" }
    /^@/ {
        gsub(re, "", \$0)
        print
        next
    }
    {
        if(\$3 ~ "^" re) {
            for (i = 1; i <= NF; i++) { gsub(re, "", \$i) }
            print
        }
    }
' | samtools view -@ ${task.cpus} -b -o "\${OUTPUT_BAM}" -
samtools index -@ ${task.cpus} "\${OUTPUT_BAM}"
echo "Finish demultiplexing for ${sp}."
"""
    }.join('\n')
    """
    ${species_cmds}
    """
}

process PREPROCESS_DEMUX {
    label 'process_medium'
    label 'scrnaseq'
    publishDir "${publish_dir}", mode: 'copy'

    input:
    tuple val(sample_id), val(species_name), val(demux_matrix_path)
    val publish_dir

    output:
    tuple val(sample_id), val(species_name), path("seur_diem_barcodes_${sample_id}.tsv.gz"), emit: barcodes
    tuple val(sample_id), val(species_name), path("seur_objs_${sample_id}.rds"),             emit: seur_objs
    tuple val(sample_id), val(species_name), path("summary_dims_${sample_id}.rds"),          emit: summary_dims
    tuple val(sample_id), val(species_name), path("summary_tbs_${sample_id}.rds"),           emit: summary_tbs
    tuple val(sample_id), val(species_name), path("summary_plts_${sample_id}.rds"),          emit: summary_plts
    tuple val(sample_id), val(species_name), path("summary_opts_${sample_id}.rds"),          emit: summary_opts
    tuple val(sample_id), val(species_name), path("sce_${sample_id}.rds"),                   emit: sce

    script:
    def display_species = species_name.replaceAll('_', ' ')
    """
    Rscript ${projectDir}/bin/preprocess_initial.R \\
        -i ${demux_matrix_path} \\
        -s ${sample_id} \\
        -S '${display_species}' \\
        -o . \\
        --ref_yaml ${params.ref_yaml} \\
        --min_diem_debris_score ${params.min_diem_debris_score} \\
        --threads ${task.cpus}
    """
}
