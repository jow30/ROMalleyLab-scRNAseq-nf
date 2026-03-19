process DEMULTIPLEX {
    tag "${sample_id}"
    label 'process_medium'
    label 'scrnaseq'
    // Publishes to demultiplex/<sp_dir>/<sample_id>/...
    publishDir "${params.out}/demultiplex", mode: 'copy'

    input:
    tuple val(sample_id), val(cellranger_outs_path)
    val species_csv

    output:
    // chi mode: demultiplex.R writes per-species RDS files
    tuple val(sample_id), path("**/*.rds"),                                       emit: rds_files,  optional: true
    // barcodes for velocyto (chi mode); may not exist in diem mode
    tuple val(sample_id), path("**/seur_diem_barcodes_${sample_id}.tsv.gz"),     emit: barcodes,   optional: true
    // demuxed count matrices: <sp_dir>/<sample_id>/raw_feature_bc_matrix
    tuple val(sample_id), path("*/*/raw_feature_bc_matrix"),                     emit: demux_matrices

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
