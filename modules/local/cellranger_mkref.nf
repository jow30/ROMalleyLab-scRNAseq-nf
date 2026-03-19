process CELLRANGER_MKREF {
    tag "${species_name}"
    label 'process_high'
    label 'cellranger'
    publishDir "${projectDir}/refs/${species_dir}", mode: 'copy'

    input:
    tuple val(species_name), path(genome_fasta), path(filtered_gtf)
    val species_dir

    output:
    tuple val(species_name), path("cellranger"), emit: cellranger_ref

    script:
    genome_name = species_name.replaceAll(' ', '_')
    """
    export MRO_DISK_SPACE_CHECK=disable
    cellranger mkref \\
        --genome=${genome_name} \\
        --fasta=${genome_fasta} \\
        --genes=${filtered_gtf} \\
        --nthreads ${task.cpus} \\
        --memgb ${task.memory.toGiga()} \\
        --localcores ${task.cpus} \\
        --localmem ${task.memory.toGiga()} \\
        --output-dir cellranger
    """
}
