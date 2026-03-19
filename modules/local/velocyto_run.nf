process VELOCYTO_RUN {
    tag "${sample_id}"
    label 'process_high'
    label 'velocyto'
    publishDir "${publish_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(cellranger_dir), path(barcodes_tsv), path(gtf), val(publish_dir)

    output:
    tuple val(sample_id), path("${sample_id}/possorted_genome_bam.loom"), emit: loom

    script:
    """
    command -v velocyto || { echo "ERROR: velocyto not found after module load"; exit 1; }

    mkdir -p ${sample_id}

    velocyto run -@ ${task.cpus} \\
        -b ${barcodes_tsv} \\
        -o ${sample_id} \\
        ${cellranger_dir}/outs/possorted_genome_bam.bam \\
        ${gtf}

    # Rename the loom file to a consistent name
    for f in ${sample_id}/possorted_genome_bam_*.loom; do
        mv "\$f" ${sample_id}/possorted_genome_bam.loom
        break
    done
    """
}
