/*
 * DEMULTIPLEX_WF — runs species demultiplexing and BAM splitting for
 * multi-species mode.  Always executes:
 *   DEMULTIPLEX_PROCESS  (demultiplex.R — writes per-species count matrices;
 *                         in chi mode also writes barcodes + RDS files)
 *   DEMUX_BAM            (splits CellRanger BAM by species)
 *
 * Included in multi_species.nf as:
 *   include { DEMULTIPLEX_WF as DEMULTIPLEX } from '...'
 * so the calling workflow still calls it as DEMULTIPLEX(...).
 *
 * Downstream channel transformations (flatMap, key construction) are
 * intentionally left to multi_species.nf.  Only direct process outputs
 * are emitted here — Nextflow requires emit: entries to reference process
 * outputs directly (def channel variables are not reliably resolved).
 */

nextflow.enable.dsl = 2

include { DEMULTIPLEX as DEMULTIPLEX_PROCESS } from '../../modules/local/demultiplex'
include { DEMUX_BAM }                          from '../../modules/local/demux_bam'

workflow DEMULTIPLEX_WF {

    take:
    cr_outs_ch    // [sample_id, cellranger_outs_path]
    species_csv   // val: comma-separated species string
    species_list  // val: list of species names

    main:
    DEMULTIPLEX_PROCESS(cr_outs_ch, species_csv)
    DEMUX_BAM(cr_outs_ch, species_list)

    emit:
    // [sample_id, path(*/cellranger/<sample_id>/outs/possorted_genome_bam.bam)] — list when >1 species
    demux_sample_dirs = DEMUX_BAM.out.demux_sample_dirs
    // [sample_id, path(*/cellranger/*/raw_feature_bc_matrix)]  — used as dependency signal in diem mode
    demux_matrices    = DEMULTIPLEX_PROCESS.out.demux_matrices
    // chi-mode outputs from demultiplex.R (optional: absent in diem mode)
    barcodes          = DEMULTIPLEX_PROCESS.out.barcodes
    rds_files         = DEMULTIPLEX_PROCESS.out.rds_files
}
