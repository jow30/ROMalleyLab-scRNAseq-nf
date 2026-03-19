/*
 * Multi-species snRNAseq workflow
 *
 * Output structure
 * ────────────────
 *   output/
 *   ├── cellranger/<sample>/
 *   └── demultiplex/
 *       └── <sp_dir>/
 *           ├── <sample>/
 *           │   ├── raw_feature_bc_matrix/    (DEMULTIPLEX process)
 *           │   └── outs/
 *           │       ├── possorted_genome_bam.bam    (DEMUX_BAM)
 *           │       └── possorted_genome_bam.bam.bai
 *           ├── preprocess/
 *           └── integration/
 *
 * Channel key convention: "${sp_dir}:::${sample_id}"
 * ":::" cannot appear in validated species/sample names.
 */

nextflow.enable.dsl = 2

include { CELLRANGER_COUNT }              from '../modules/local/cellranger_count'
include { DIEM_DEBRIS_REMOVAL }           from '../modules/local/diem_debris_removal'
include { SUMMARY_REPORT }               from '../modules/local/summary_report'
include { INTEGRATION }                  from '../modules/local/integration'
include { DEMULTIPLEX_WF as DEMULTIPLEX } from '../subworkflows/local/demultiplex'
include { PREPROCESS }                   from '../subworkflows/local/preprocess'

workflow MULTI_SPECIES_WF {
    take:
    sample_ch          // [ sample_id, fastq_dir, fastq_prefix ]
    transcriptome_ch   // path to combined cellranger reference
    species_list       // list of species names
    clean_method       // 'diem' or 'chi'  — plain Groovy string
    anno_map           // Map[ species_name -> [gtf, markers, seurat_ref] ]

    main:
    def species_csv = species_list.join(',')
    def demux_base  = "${params.out}/demultiplex"

    // ── CellRanger count ──────────────────────────────────────────────────────
    CELLRANGER_COUNT(sample_ch, transcriptome_ch)

    def cr_barrier_signal = CELLRANGER_COUNT.out.cellranger_out
        .map { sid, dir -> sid }
        .collect()
        .map { sids -> [sids] }

    def cr_barrier_ch = CELLRANGER_COUNT.out.cellranger_out
        .combine(cr_barrier_signal)
        .map { sid, dir, signal -> [sid, dir] }

    def cr_outs_ch = cr_barrier_ch.map { sid, dir ->
        [sid, "${dir.toAbsolutePath()}/outs"]
    }

    // ── Demultiplex subworkflow (DEMULTIPLEX_PROCESS + DEMUX_BAM) ─────────────
    DEMULTIPLEX(cr_outs_ch, species_csv, species_list)

    // Flatten DEMUX_BAM output: one [key, demux_sample_dir] per (species, sample).
    // key = "${sp_dir}:::${sample_id}" — unique across all species/sample combos.
    // Note: flatMap is done here (not in the subworkflow emit block) because
    // Nextflow requires emit: entries to reference process outputs directly.
    def demux_dirs_flat = DEMULTIPLEX.out.demux_sample_dirs
        .flatMap { sid, dirs ->
            (dirs instanceof List ? dirs : [dirs]).collect { d ->
                ["${d.parent.name}:::${sid}", d]
            }
        }

    // Lookup channel: [key, sample_id] — recovers the plain sample_id from
    // the composite key without string-splitting in downstream closures.
    def key_to_sid_ch = demux_dirs_flat.map { key, d -> [key, d.name] }

    // ── Build unified preprocess_rds_ch (diem vs chi) ────────────────────────
    // Both branches produce: [key, barcodes_tsv, seur_objs, dims, tbs, plts]
    def preprocess_rds_ch

    if (clean_method == 'diem') {
        // Run DIEM_DEBRIS_REMOVAL on each demuxed species/sample matrix.
        // DEMULTIPLEX_PROCESS has already published raw_feature_bc_matrix there;
        // preprocess_initial.R reads it from the published path.
        def diem_input_ch = demux_dirs_flat
            .map { key, d ->
                def sp_dir_name = d.parent.name
                def sid         = d.name
                def sp          = species_list.find { it.replaceAll(' ', '_') == sp_dir_name }
                def display_sp  = sp ?: sp_dir_name.replaceAll('_', ' ')
                def matrix_dir  = file("${demux_base}/${sp_dir_name}/${sid}").toAbsolutePath().toString()
                def pub_dir     = "${demux_base}/${sp_dir_name}/preprocess"
                [sid, sp_dir_name, matrix_dir, display_sp, pub_dir]
            }

        DIEM_DEBRIS_REMOVAL(diem_input_ch)

        preprocess_rds_ch = DIEM_DEBRIS_REMOVAL.out.barcodes
            .join(DIEM_DEBRIS_REMOVAL.out.seur_objs,    by: [0, 1])
            .join(DIEM_DEBRIS_REMOVAL.out.summary_dims,  by: [0, 1])
            .join(DIEM_DEBRIS_REMOVAL.out.summary_tbs,   by: [0, 1])
            .join(DIEM_DEBRIS_REMOVAL.out.summary_plts,  by: [0, 1])
            .map { sid, sp_dir, bc, seur, dims, tbs, plts ->
                ["${sp_dir}:::${sid}", bc, seur, dims, tbs, plts]
            }

    } else {
        // chi: demultiplex.R already ran preprocessing — barcodes + RDS come
        // directly from DEMULTIPLEX output
        def bc_flat = DEMULTIPLEX.out.barcodes
            .flatMap { sid, bcs ->
                (bcs instanceof List ? bcs : [bcs]).collect { bc ->
                    def sp_dir_name = bc.parent.parent.name
                    ["${sp_dir_name}:::${sid}", bc]
                }
            }

        def rds_flat = DEMULTIPLEX.out.rds_files
            .flatMap { sid, rds_list ->
                def entries = []
                for (sp in species_list) {
                    def sp_dir_name  = sp.replaceAll(' ', '_')
                    def sp_rds       = (rds_list instanceof List ? rds_list : [rds_list])
                        .findAll { it.parent.parent.name == sp_dir_name }
                    def seur_obj     = sp_rds.find { it.name == "seur_objs_${sid}.rds" }
                    def summary_dims = sp_rds.find { it.name == "summary_dims_${sid}.rds" }
                    def summary_tbs  = sp_rds.find { it.name == "summary_tbs_${sid}.rds" }
                    def summary_plts = sp_rds.find { it.name == "summary_plts_${sid}.rds" }
                    if (seur_obj && summary_dims && summary_tbs && summary_plts) {
                        entries << ["${sp_dir_name}:::${sid}", seur_obj, summary_dims, summary_tbs, summary_plts]
                    }
                }
                return entries
            }

        preprocess_rds_ch = bc_flat
            .join(rds_flat)
            .map { key, bc, seur, dims, tbs, plts -> [key, bc, seur, dims, tbs, plts] }
    }

    // ── Per-item channels for PREPROCESS ─────────────────────────────────────
    def gtf_ch = demux_dirs_flat.map { key, d ->
        def sp_dir_name = d.parent.name
        def sp = species_list.find { it.replaceAll(' ', '_') == sp_dir_name }
        [key, file(anno_map[sp].gtf)]
    }

    def species_ch = demux_dirs_flat.map { key, d ->
        def sp_dir_name = d.parent.name
        def sp = species_list.find { it.replaceAll(' ', '_') == sp_dir_name }
        [key, sp]
    }

    def preprocess_pub_ch = demux_dirs_flat.map { key, d ->
        def sp_dir_name = d.parent.name
        [key, "${demux_base}/${sp_dir_name}/preprocess"]
    }

    def cellranger_pub_ch = demux_dirs_flat.map { key, d ->
        def sp_dir_name = d.parent.name
        [key, "${demux_base}/${sp_dir_name}/cellranger"]
    }

    // ── Shared downstream subworkflow ─────────────────────────────────────────
    PREPROCESS(
        demux_dirs_flat,
        preprocess_rds_ch,
        gtf_ch,
        species_ch,
        preprocess_pub_ch,
        cellranger_pub_ch
    )

    // ── Per-species summary + integration ─────────────────────────────────────
    def cr_dirs_collected = cr_barrier_ch
        .map { sid, dir -> dir.toAbsolutePath().toString() }
        .collect()
    def cr_sels_collected = cr_barrier_ch
        .map { sid, dir -> sid }
        .collect()

    def summary_ch     = Channel.empty()
    def integration_ch = Channel.empty()

    for (sp in species_list) {
        def sp_dir        = sp.replaceAll(' ', '_')
        def sp_preprocess = "${demux_base}/${sp_dir}/preprocess"
        def anno          = anno_map[sp]

        // Filter to this species, recover plain sample_id via key_to_sid_ch
        def sp_seur_sels = PREPROCESS.out.seur_clean
            .filter { key, rds -> key.startsWith("${sp_dir}:::") }
            .join(key_to_sid_ch)
            .map    { key, rds, sid -> sid }
            .collect()

        def sp_summary = cr_dirs_collected.map { dirs -> [dirs] }
            .combine(cr_sels_collected.map { sels -> [sels] })
            .combine(sp_seur_sels.map { sels -> [sels] })
            .map { cr_dirs, cr_sels, sp_sels ->
                [cr_dirs, cr_sels, [file(sp_preprocess).toAbsolutePath().toString()], sp_sels, sp, sp_preprocess]
            }
        summary_ch = summary_ch.mix(sp_summary)

        def sp_rds = PREPROCESS.out.seur_clean
            .filter { key, rds -> key.startsWith("${sp_dir}:::") }
            .map    { key, rds -> rds.toAbsolutePath().toString() }
            .collect()
            .map    { rds_list ->
                [rds_list, anno.markers, anno.seurat_ref, "${demux_base}/${sp_dir}/integration"]
            }
        integration_ch = integration_ch.mix(sp_rds)
    }

    SUMMARY_REPORT(summary_ch)
    INTEGRATION(integration_ch)
}
