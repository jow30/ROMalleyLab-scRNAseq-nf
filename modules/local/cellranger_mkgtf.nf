process CELLRANGER_MKGTF {
    tag "${species_name}"
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
