/*
 * PREPROCESS — shared downstream workflow used by both single- and
 * multi-species modes.  Runs everything from SEPARATE_READS through
 * CELL_FILTERING and emits seur_clean for SUMMARY_REPORT / INTEGRATION,
 * which are handled at the calling-workflow level.
 *
 * Channel key convention
 * ──────────────────────
 *   Single-species : key = sample_id                    (plain string)
 *   Multi-species  : key = "${sp_dir}:::${sample_id}"  (composite string)
 *
 * ":::" cannot appear in validated sample/species names, making it an
 * unambiguous separator.  All input channels share this key so .join() works
 * uniformly regardless of mode.
 */

nextflow.enable.dsl = 2

include { SEPARATE_READS }  from '../../modules/local/separate_reads'
include { COUNT_READS }     from '../../modules/local/count_reads'
include { VELOCYTO_RUN }    from '../../modules/local/velocyto_run'
include { UNSPLICE_RATIO }  from '../../modules/local/unsplice_ratio'
include { CELL_FILTERING }  from '../../modules/local/cell_filtering'

workflow PREPROCESS {

    take:
    // [key, cellranger_dir]   — dir must have outs/ with bam + matrix
    cr_dirs_ch

    // [key, barcodes_tsv, seur_objs_rds, summary_dims_rds, summary_tbs_rds, summary_plts_rds]
    // — precomputed by DIEM_DEBRIS_REMOVAL (single) or by the DEMULTIPLEX subworkflow (multi)
    preprocess_rds_ch

    // [key, gtf_path]
    gtf_ch

    // [key, species_name]   — display string, e.g. "Arabidopsis thaliana"
    species_ch

    // [key, preprocess_publish_dir]   — where preprocess outputs land
    preprocess_pub_ch

    // [key, cellranger_publish_dir]   — where velocyto loom lands
    cellranger_pub_ch

    main:

    // ── SEPARATE_READS ────────────────────────────────────────────────────────
    SEPARATE_READS(cr_dirs_ch.join(preprocess_pub_ch))

    def all_bams_ch = SEPARATE_READS.out.mapped_bam
        .join(preprocess_pub_ch)
        .map { key, bam, pub -> [key, 'mapped', bam, pub] }
        .mix(
            SEPARATE_READS.out.uniq_bam
                .join(preprocess_pub_ch)
                .map { key, bam, pub -> [key, 'uniq', bam, pub] },
            SEPARATE_READS.out.multi_bam
                .join(preprocess_pub_ch)
                .map { key, bam, pub -> [key, 'multi', bam, pub] }
        )

    // ── COUNT_READS ───────────────────────────────────────────────────────────
    COUNT_READS(all_bams_ch)

    // ── VELOCYTO_RUN ──────────────────────────────────────────────────────────
    def barcodes_ch = preprocess_rds_ch.map { key, bc, s, d, t, p -> [key, bc] }

    def velo_ch = cr_dirs_ch
        .join(barcodes_ch)
        .join(gtf_ch)
        .join(cellranger_pub_ch)
        .map { key, cr_dir, bc, gtf, pub -> [key, cr_dir, bc, gtf, pub] }

    VELOCYTO_RUN(velo_ch)

    // ── UNSPLICE_RATIO ────────────────────────────────────────────────────────
    def unsplice_ch = VELOCYTO_RUN.out.loom
        .join(preprocess_pub_ch)
        .map { key, loom, pub -> [key, loom, pub] }

    UNSPLICE_RATIO(unsplice_ch)

    // ── CELL_FILTERING ────────────────────────────────────────────────────────
    def rds_ch = preprocess_rds_ch
        .map { key, bc, seur, dims, tbs, plts -> [key, seur, dims, tbs, plts] }

    def cell_filter_ch = UNSPLICE_RATIO.out.unsplice_txt
        .join(rds_ch)
        .join(species_ch)
        .join(preprocess_pub_ch)
        .map { key, txt, seur, dims, tbs, plts, sp, pub ->
            [key, sp, txt, seur, dims, tbs, plts, pub]
        }

    CELL_FILTERING(cell_filter_ch)

    emit:
    seur_clean   = CELL_FILTERING.out.seur_clean   // [key, rds]
    seur_objs    = CELL_FILTERING.out.seur_objs
    summary_dims = CELL_FILTERING.out.summary_dims
    summary_tbs  = CELL_FILTERING.out.summary_tbs
    summary_plts = CELL_FILTERING.out.summary_plts
    summary_opts = CELL_FILTERING.out.summary_opts
}
