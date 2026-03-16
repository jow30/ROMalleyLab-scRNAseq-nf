/*
 * Full preprocessing with velocyto results module
 *
 * PREPROCESS_WITH_VELOCYTO       : single-species; publish_dir is the preprocess dir
 * PREPROCESS_WITH_VELOCYTO_FOR_SPECIES : multi-species; keyed by [sp_dir, sample_id];
 *                                   publishes to ${publish_dir_base}/${sp_dir}
 */

process PREPROCESS_WITH_VELOCYTO {
    label 'process_medium'
    label 'scrnaseq'
    publishDir "${publish_dir}", mode: 'copy'

    input:
    tuple val(sample_id), val(velocyto_txt_path),
          path(seur_objs_rds), path(summary_dims_rds), path(summary_tbs_rds), path(summary_plts_rds)
    val species_name
    val publish_dir

    output:
    tuple val(sample_id), path("seur_clean_${sample_id}.rds"),       emit: seur_clean
    tuple val(sample_id), path("seur_objs_${sample_id}.rds"),        emit: seur_objs
    tuple val(sample_id), path("summary_dims_${sample_id}.rds"),     emit: summary_dims
    tuple val(sample_id), path("summary_tbs_${sample_id}.rds"),      emit: summary_tbs
    tuple val(sample_id), path("summary_plts_${sample_id}.rds"),     emit: summary_plts
    tuple val(sample_id), path("summary_opts_${sample_id}.rds"),     emit: summary_opts

    script:
    """
    Rscript ${projectDir}/bin/preprocess_velocyto.R \\
        -s ${sample_id} \\
        -S '${species_name}' \\
        -o . \\
        -v ${velocyto_txt_path} \\
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
        --threads ${task.cpus}
    """
}

process PREPROCESS_WITH_VELOCYTO_FOR_SPECIES {
    label 'process_medium'
    label 'scrnaseq'
    publishDir "${publish_dir_base}/${sp_dir}", mode: 'copy'

    input:
    tuple val(sp_dir), val(sample_id), val(species_name), val(velocyto_txt_path),
          path(seur_objs_rds), path(summary_dims_rds), path(summary_tbs_rds), path(summary_plts_rds)
    val publish_dir_base

    output:
    tuple val(sp_dir), val(sample_id), path("seur_clean_${sample_id}.rds"),   emit: seur_clean
    tuple val(sp_dir), val(sample_id), path("seur_objs_${sample_id}.rds"),    emit: seur_objs
    tuple val(sp_dir), val(sample_id), path("summary_dims_${sample_id}.rds"), emit: summary_dims
    tuple val(sp_dir), val(sample_id), path("summary_tbs_${sample_id}.rds"),  emit: summary_tbs
    tuple val(sp_dir), val(sample_id), path("summary_plts_${sample_id}.rds"), emit: summary_plts
    tuple val(sp_dir), val(sample_id), path("summary_opts_${sample_id}.rds"), emit: summary_opts

    script:
    """
    Rscript ${projectDir}/bin/preprocess_velocyto.R \\
        -s ${sample_id} \\
        -S '${species_name}' \\
        -o . \\
        -v ${velocyto_txt_path} \\
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
        --threads ${task.cpus}
    """
}
