/*
 * Integration module
 */

process INTEGRATION {
    label 'process_medium'
    label 'scrnaseq'
    publishDir "${out_dir}", mode: 'copy'

    input:
    val input_rds_paths  // list of absolute paths to seur_clean_*.rds files
    val markers_file     // path to markers CSV (can be null/empty)
    val seurat_ref_file  // path to Seurat reference RDS (can be null/empty)
    val out_dir          // output directory for publishDir

    output:
    path "seuratObjs_integrated.rds",       emit: integrated_rds
    path "*.png",                            emit: plots, optional: true
    path "*.svg",                            emit: svg_plots, optional: true
    path "*.csv",                            emit: csv_files, optional: true

    script:
    rds_csv = input_rds_paths.join(',')
    markers_arg = markers_file ? "--markers ${markers_file}" : ''
    ref_arg     = seurat_ref_file ? "--seurat_reference ${seurat_ref_file}" : ''
    """
    Rscript ${projectDir}/bin/integration.R \\
        -i "${rds_csv}" \\
        -o . \\
        --nFeatures ${params.nFeatures} \\
        ${markers_arg} \\
        ${ref_arg}
    """
}
