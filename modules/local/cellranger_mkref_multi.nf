process CELLRANGER_MKREF_MULTI {
    tag "${species_dir}"
    label 'process_high'
    label 'cellranger'
    publishDir "${projectDir}/refs/${species_dir}", mode: 'copy'

    input:
    val genome_args   // list of '--genome=X --fasta=Y --genes=Z' strings
    val species_dir

    output:
    path "cellranger", emit: cellranger_ref

    script:
    genome_args_str = genome_args.join(' \\\n        ')
    """
    export MRO_DISK_SPACE_CHECK=disable
    cellranger mkref \\
        ${genome_args_str} \\
        --nthreads ${task.cpus} \\
        --memgb ${task.memory.toGiga()} \\
        --localcores ${task.cpus} \\
        --localmem ${task.memory.toGiga()} \\
        --output-dir cellranger
    """
}
