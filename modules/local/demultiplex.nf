process DEMULTIPLEX {
    tag "${sample_id}"
    label 'process_medium'
    label 'scrnaseq'
    publishDir "${params.out}/demultiplex", mode: 'copy', saveAs: { fn -> fn.tokenize('/').last().startsWith('summary_') ? null : fn }

    input:
    tuple val(sample_id), val(cellranger_outs_path)
    val species_csv

    output:
    // chi mode: demultiplex.R writes per-species RDS files → moved to <sp_dir>/preprocess/
    tuple val(sample_id), path("*/preprocess/*.rds"),                                          emit: rds_files,      optional: true
    // barcodes for velocyto (chi mode) → moved to <sp_dir>/preprocess/
    tuple val(sample_id), path("*/preprocess/seur_diem_barcodes_${sample_id}.tsv.gz"),        emit: barcodes,       optional: true
    // demuxed count matrices → moved to <sp_dir>/cellranger/<sample_id>/outs/raw_feature_bc_matrix
    tuple val(sample_id), path("*/cellranger/*/outs/raw_feature_bc_matrix"),                  emit: demux_matrices

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

    # Reorganise demultiplex.R output to match the target directory layout:
    #   <sp_dir>/<sample_id>/raw_feature_bc_matrix  →  <sp_dir>/cellranger/<sample_id>/raw_feature_bc_matrix
    #   <sp_dir>/<sample_id>/<chi-mode files>        →  <sp_dir>/preprocess/<chi-mode files>
    for sp_dir in */; do
        sp_dir="\${sp_dir%/}"
        [ -d "\$sp_dir" ] || continue
        mkdir -p "\${sp_dir}/cellranger" "\${sp_dir}/preprocess"
        for sample_dir in "\${sp_dir}"/*/; do
            [ -d "\$sample_dir" ] || continue
            sid="\$(basename "\$sample_dir")"
            [ "\$sid" = "cellranger" ] && continue
            [ "\$sid" = "preprocess" ] && continue
            # count matrix → cellranger/<sample_id>/outs/
            if [ -d "\${sample_dir}raw_feature_bc_matrix" ]; then
                mkdir -p "\${sp_dir}/cellranger/\${sid}/outs"
                mv "\${sample_dir}raw_feature_bc_matrix" "\${sp_dir}/cellranger/\${sid}/outs/"
            fi
            # remaining files (chi mode: barcodes, RDS) → preprocess/
            for f in "\${sample_dir}"*; do
                [ -e "\$f" ] && mv "\$f" "\${sp_dir}/preprocess/"
            done
            rmdir "\${sample_dir}" 2>/dev/null || true
        done
        # chi mode: RDS + barcodes written directly to sp_dir by demultiplex.R → move to preprocess/
        for f in "\${sp_dir}"/*.rds "\${sp_dir}"/*.tsv.gz; do
            [ -e "\$f" ] && mv "\$f" "\${sp_dir}/preprocess/"
        done
    done
    """
}
