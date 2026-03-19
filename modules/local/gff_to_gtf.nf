process GFF_TO_GTF {
    tag "${species_name}"
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
