/*
 * Reference preparation module
 */

process GFF_TO_GTF {
    label 'process_medium'
    label 'agat'

    input:
    tuple val(species_name), path(gff)

    output:
    tuple val(species_name), path("${basename}.gtf"), emit: gtf

    script:
    basename = gff.baseName.replaceAll(/\.gff3?$/, '').replaceAll(/\.gene$/, '')
    """
    agat_convert_sp_gff2gtf.pl --gff ${gff} -o ${basename}.gtf
    """
}

process CELLRANGER_MKGTF {
    label 'process_low'
    label 'cellranger'

    input:
    tuple val(species_name), path(gtf)

    output:
    tuple val(species_name), path("${gtf.baseName}.flt.gtf"), emit: filtered_gtf

    script:
    """
    cellranger mkgtf \\
        ${gtf} \\
        ${gtf.baseName}.flt.gtf \\
        --attribute=gene_biotype:protein_coding \\
        --attribute=gene_biotype:lncRNA \\
        --attribute=gene_biotype:miRNA \\
        --attribute=gene_biotype:snoRNA \\
        --attribute=gene_biotype:snRNA \\
        --attribute=transcript_biotype:protein_coding \\
        --attribute=transcript_biotype:lncRNA \\
        --attribute=transcript_biotype:miRNA \\
        --attribute=transcript_biotype:snoRNA \\
        --attribute=transcript_biotype:snRNA
    """
}

process CELLRANGER_MKREF {
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

process CELLRANGER_MKREF_MULTI {
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
