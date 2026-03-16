/*
 * Velocyto module: RNA velocity analysis
 */

process VELOCYTO_RUN {
    label 'process_medium'
    label 'velocyto'
    publishDir "${publish_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(cellranger_dir), path(barcodes_tsv)
    path gtf
    val publish_dir

    output:
    tuple val(sample_id), path("${sample_id}/possorted_genome_bam.loom"), emit: loom

    script:
    """
    module load velocyto/0.17.17
    module load samtools/1.22.1
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

process UNSPLICE_RATIO {
    label 'process_low'
    label 'velocyto'
    publishDir "${publish_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(loom_file)
    val publish_dir

    output:
    tuple val(sample_id), path("${sample_id}.txt"), emit: unsplice_txt

    script:
    """
    module load velocyto/0.17.17
    python ${projectDir}/bin/unsplice_ratio.py ${loom_file} ${sample_id}.txt
    """
}

/*
 * Multi-species velocyto: runs on demultiplexed BAM files directly
 * (no outs/ subdirectory structure). Publishes to
 * ${publish_dir_base}/${sp_dir}/${sample_id}/
 */
process VELOCYTO_RUN_MULTI {
    label 'process_medium'
    label 'velocyto'
    publishDir "${publish_dir_base}/${sp_dir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sp_dir), val(sample_id), path(bam_file), path(bam_index), path(barcodes_tsv), path(gtf_file)
    val publish_dir_base

    output:
    tuple val(sp_dir), val(sample_id), path("possorted_genome_bam.loom"), emit: loom

    script:
    """
    module load velocyto/0.17.17
    module load samtools/1.22.1
    command -v velocyto || { echo "ERROR: velocyto not found after module load"; exit 1; }

    velocyto run -@ ${task.cpus} \\
        -b ${barcodes_tsv} \\
        -o . \\
        ${bam_file} \\
        ${gtf_file}

    # Rename the loom file to a consistent name
    for f in *.loom; do
        [ "\$f" = "possorted_genome_bam.loom" ] && break
        mv "\$f" possorted_genome_bam.loom
        break
    done
    """
}

process UNSPLICE_RATIO_MULTI {
    label 'process_low'
    label 'velocyto'
    publishDir "${publish_dir_base}/${sp_dir}", mode: 'copy'

    input:
    tuple val(sp_dir), val(sample_id), path(loom_file)
    val publish_dir_base

    output:
    tuple val(sp_dir), val(sample_id), path("${sample_id}.txt"), emit: unsplice_txt

    script:
    """
    module load velocyto/0.17.17
    python ${projectDir}/bin/unsplice_ratio.py ${loom_file} ${sample_id}.txt
    """
}
