/*
 * DIEM_DEBRIS_REMOVAL — runs preprocess_initial.R to perform DIEM-based
 * ambient RNA / debris filtering on a raw feature-barcode matrix.
 *
 * Used in two contexts:
 *   Single-species : sp_dir = '' (empty string)
 *   Multi-species  : sp_dir = species directory name (e.g. "Arabidopsis_thaliana")
 *
 * sp_dir is carried through the tuple so Nextflow channel keys stay unique
 * when the same sample_id is processed for multiple species simultaneously.
 * The R script always receives the original sample_id for output file naming.
 */

process DIEM_DEBRIS_REMOVAL {
    tag { sp_dir ? "${sp_dir}/${sample_id}" : sample_id }
    label 'process_medium'
    label 'scrnaseq'
    publishDir "${publish_dir}", mode: 'copy', saveAs: { fn -> fn.startsWith('summary_') || fn.startsWith('seur_diem_barcodes_') || fn.startsWith('seur_objs_') ? null : fn }

    input:
    tuple val(sample_id), val(sp_dir), val(input_path), val(species_name), val(publish_dir)

    output:
    tuple val(sample_id), val(sp_dir), path("seur_diem_barcodes_${sample_id}.tsv.gz"), emit: barcodes
    tuple val(sample_id), val(sp_dir), path("seur_objs_${sample_id}.rds"),             emit: seur_objs
    tuple val(sample_id), val(sp_dir), path("summary_dims_${sample_id}.rds"),          emit: summary_dims
    tuple val(sample_id), val(sp_dir), path("summary_tbs_${sample_id}.rds"),           emit: summary_tbs
    tuple val(sample_id), val(sp_dir), path("summary_plts_${sample_id}.rds"),          emit: summary_plts
    tuple val(sample_id), val(sp_dir), path("summary_opts_${sample_id}.rds"),          emit: summary_opts

    script:
    def script_hash = file("${projectDir}/bin/preprocess_initial.R").text.md5()
    """
    # script_hash: ${script_hash}
    Rscript ${projectDir}/bin/preprocess_initial.R \\
        -i ${input_path} \\
        -s ${sample_id} \\
        -S '${species_name}' \\
        -o . \\
        --ref_yaml ${params.ref_yaml} \\
        --min_diem_debris_score ${params.min_diem_debris_score} \\
        --threads ${task.cpus}
    """
}
