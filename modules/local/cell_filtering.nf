process CELL_FILTERING {
    tag "${sample_id}"
    label 'process_high'
    label 'scrnaseq'
    publishDir "${publish_dir}", mode: 'copy'

    input:
    tuple val(key), val(sample_id), val(species_name), path(velocyto_txt),
          path(seur_objs_rds), path(summary_dims_rds), path(summary_tbs_rds), path(summary_plts_rds),
          val(publish_dir)

    output:
    tuple val(key), path("seur_clean_${sample_id}.rds"),       emit: seur_clean,    optional: true
    tuple val(key), path("seur_objs_${sample_id}.rds"),        emit: seur_objs,     optional: true
    tuple val(key), path("summary_dims_${sample_id}.rds"),     emit: summary_dims,  optional: true
    tuple val(key), path("summary_tbs_${sample_id}.rds"),      emit: summary_tbs,   optional: true
    tuple val(key), path("summary_plts_${sample_id}.rds"),     emit: summary_plts,  optional: true
    tuple val(key), path("summary_opts_${sample_id}.rds"),     emit: summary_opts,  optional: true

    script:
    """
    Rscript ${projectDir}/bin/preprocess_velocyto.R \\
        -s ${sample_id} \\
        -S '${species_name}' \\
        -o . \\
        -v ${velocyto_txt} \\
        --ref_yaml ${params.ref_yaml} \\
        --min_unsplice_ratio ${params.min_unsplice_ratio} \\
        --min_nCount_RNA ${params.min_nCount_RNA} \\
        --min_nFeature_RNA ${params.min_nFeature_RNA} \\
        --max_mt ${params.max_mt} \\
        --max_cp ${params.max_cp} \\
        --nHVG ${params.nHVG} \\
        --min_ncell_expr ${params.min_ncell_expr} \\
        --remove_doublet ${params.remove_doublet ? 'TRUE' : 'FALSE'} \\
        --max_doublet_score ${params.max_doublet_score} \\
        --min_nClusterMarker ${params.min_nClusterMarker} \\
        --min_cells ${params.min_cells} \\
        --threads ${task.cpus}
    """
}
