/*
 * CellRanger count module
 */

process CELLRANGER_COUNT {
    label 'process_cellranger'
    label 'cellranger'
    publishDir "${params.out}/cellranger", mode: 'copy'

    input:
    tuple val(sample_id), val(fastq_dir), val(fastq_prefix)
    path transcriptome

    output:
    tuple val(sample_id), path("${sample_id}"), emit: cellranger_out

    script:
    """
    export MRO_DISK_SPACE_CHECK=disable
    cellranger count \\
        --create-bam true \\
        --id=${sample_id} \\
        --fastqs=${fastq_dir} \\
        --sample=${fastq_prefix} \\
        --transcriptome=${transcriptome} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()}
    """
}
