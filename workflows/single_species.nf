/*
 * Single-species snRNAseq workflow
 */

nextflow.enable.dsl = 2

include { CELLRANGER_COUNT }   from '../modules/local/cellranger_count'
include { DIEM_DEBRIS_REMOVAL } from '../modules/local/diem_debris_removal'
include { SUMMARY_REPORT }     from '../modules/local/summary_report'
include { INTEGRATION }        from '../modules/local/integration'
include { PREPROCESS }         from '../subworkflows/local/preprocess'

workflow SINGLE_SPECIES_WF {
    take:
    sample_ch          // [ sample_id, fastq_dir, fastq_prefix ]
    transcriptome_ch   // path to cellranger reference
    species_name       // string
    gtf_file           // path to GTF file for velocyto
    markers            // path to markers CSV, or '' if none
    seurat_ref         // path to Seurat reference RDS, or '' if none

    main:
    def preprocess_dir = "${params.out}/preprocess"
    def cellranger_dir = "${params.out}/cellranger"

    // ── CellRanger count ──────────────────────────────────────────────────────
    CELLRANGER_COUNT(sample_ch, transcriptome_ch)

    // Barrier: re-emit each per-sample tuple only after ALL CellRanger jobs finish.
    def cr_barrier_signal = CELLRANGER_COUNT.out.cellranger_out
        .map { sid, dir -> sid }
        .collect()
        .map { sids -> [sids] }

    def cr_barrier_ch = CELLRANGER_COUNT.out.cellranger_out
        .combine(cr_barrier_signal)
        .map { sid, dir, signal -> [sid, dir] }

    // ── DIEM debris removal ───────────────────────────────────────────────────
    // sp_dir = '' for single-species (keeps tuple shape consistent with multi)
    def diem_input_ch = cr_barrier_ch.map { sid, dir ->
        [sid, '', "${dir.toAbsolutePath()}/outs", species_name, preprocess_dir]
    }

    DIEM_DEBRIS_REMOVAL(diem_input_ch)

    // Build the unified preprocess_rds channel required by PREPROCESS:
    // [key, barcodes_tsv, seur_objs_rds, summary_dims_rds, summary_tbs_rds, summary_plts_rds]
    def preprocess_rds_ch = DIEM_DEBRIS_REMOVAL.out.barcodes
        .join(DIEM_DEBRIS_REMOVAL.out.seur_objs,    by: [0, 1])
        .join(DIEM_DEBRIS_REMOVAL.out.summary_dims,  by: [0, 1])
        .join(DIEM_DEBRIS_REMOVAL.out.summary_tbs,   by: [0, 1])
        .join(DIEM_DEBRIS_REMOVAL.out.summary_plts,  by: [0, 1])
        .map { sid, sp_dir, bc, seur, dims, tbs, plts ->
            [sid, bc, seur, dims, tbs, plts]
        }

    // Per-sample publish-dir channels (same dir for all samples in single-species)
    def preprocess_pub_ch  = cr_barrier_ch.map { sid, dir -> [sid, preprocess_dir] }
    def cellranger_pub_ch  = cr_barrier_ch.map { sid, dir -> [sid, cellranger_dir] }
    def gtf_ch             = cr_barrier_ch.map { sid, dir -> [sid, gtf_file] }
    def species_ch         = cr_barrier_ch.map { sid, dir -> [sid, species_name] }

    // ── Shared downstream subworkflow ─────────────────────────────────────────
    PREPROCESS(
        cr_barrier_ch,
        preprocess_rds_ch,
        gtf_ch,
        species_ch,
        preprocess_pub_ch,
        cellranger_pub_ch
    )

    // ── Summary report ────────────────────────────────────────────────────────
    def cr_dirs_ch = cr_barrier_ch
        .map { sid, dir -> dir.toAbsolutePath().toString() }
        .collect()
    def cr_sels_ch = cr_barrier_ch
        .map { sid, dir -> sid }
        .collect()
    def seur_sels_ch = PREPROCESS.out.seur_clean
        .map { sid, rds -> sid }
        .collect()

    def summary_ch = cr_dirs_ch.map { dirs -> [dirs] }
        .combine(cr_sels_ch.map { sels -> [sels] })
        .combine(seur_sels_ch.map { sels -> [sels] })
        .map { cr_dirs, cr_sels, seur_sels ->
            [cr_dirs, cr_sels, [file(preprocess_dir).toAbsolutePath().toString()], seur_sels, species_name, preprocess_dir]
        }
    SUMMARY_REPORT(summary_ch)

    // ── Integration ───────────────────────────────────────────────────────────
    def rds_paths_ch = PREPROCESS.out.seur_clean
        .map { sid, rds -> rds.toAbsolutePath().toString() }
        .collect()

    def integration_ch = rds_paths_ch.map { rds_paths ->
        [rds_paths, markers, seurat_ref, "${params.out}/integration"]
    }
    INTEGRATION(integration_ch)
}
