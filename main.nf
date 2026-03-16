#!/usr/bin/env nextflow

/*
 * ============================================================================
 *  snRNAseq Processing Pipeline for Plants
 * ============================================================================
 *  Pipeline for processing single-nucleus RNA sequencing data from plants.
 *  Orchestrates: reference preparation → cellranger count → filtering/
 *  demultiplexing → read counting → velocyto → full preprocessing →
 *  summary report → integration.
 * ============================================================================
 */

nextflow.enable.dsl = 2

// ── Module includes ─────────────────────────────────────────────────────────

include { GFF_TO_GTF; CELLRANGER_MKGTF; CELLRANGER_MKREF; CELLRANGER_MKREF_MULTI } from './modules/ref_prep'
include { CELLRANGER_COUNT }          from './modules/cellranger'
include { PREPROCESS_INITIAL; DEMULTIPLEX; DEMUX_BAM; PREPROCESS_DEMUX } from './modules/preprocess'
include { SEPARATE_READS; COUNT_READS_2KB } from './modules/read_counts'
include {
    SEPARATE_READS as SEPARATE_READS_MULTI;
    COUNT_READS_2KB as COUNT_READS_2KB_MULTI
} from './modules/read_counts'
include { VELOCYTO_RUN; UNSPLICE_RATIO; VELOCYTO_RUN_MULTI; UNSPLICE_RATIO_MULTI } from './modules/velocyto'
include { PREPROCESS_WITH_VELOCYTO; PREPROCESS_WITH_VELOCYTO_FOR_SPECIES } from './modules/preprocess_velocyto'
include { SUMMARY_REPORT }                 from './modules/summary'
include { SUMMARY_REPORT as SUMMARY_REPORT_MULTI } from './modules/summary'
include { INTEGRATION }                    from './modules/integration'

// ── Help message ────────────────────────────────────────────────────────────

def helpMessage() {
    log.info """
    =========================================================
     snRNAseq Processing Pipeline  v${workflow.manifest.version}
    =========================================================

    Usage:
      nextflow run main.nf --input <samplesheet.csv> --out <output_dir> [options]

    Mandatory:
      --input           Path to samplesheet CSV (columns: sample, fastq_dir, fastq_R1, fastq_R2)
      --out             Output directory path
      --species         Species name(s), comma-separated (default: 'Arabidopsis thaliana')

    Reference:
      --genome          Path(s) to genome FASTA (auto-resolved from species if available)
      --gtf             Path(s) to GTF file(s) (auto-resolved from species if available)
      --ref_yaml        Path to scQC.yaml (default: refs/scQC.yaml)

    Filtering (preprocess.R):
      --min_diem_debris_score   Min debris score for DIEM (default: 1)
      --min_unsplice_ratio      Min unsplice ratio (default: 0.1)
      --min_nCount_RNA          Min UMI per cell (default: 400)
      --min_nFeature_RNA        Min genes per cell (default: 300)
      --max_mt                  Max percent MT (default: 10)
      --max_cp                  Max percent CP (default: 15)
      --nHVG                    Number of HVGs (default: 3000)
      --min_ncell_expr          Min cells expressing gene (default: 5)
      --remove_doublet          Remove doublets (default: true)
      --max_doublet_score       Max doublet pANN (default: 0.4)
      --min_nClusterMarker      Min cluster markers (default: 5)

    Demultiplexing (multi-species):
      --clean                         Method: 'diem' or 'chi' (auto-selected)
      --min_UMI_per_cell_barcode      Min UMI per barcode for chi (default: 400)
      --chisq_pvalues_max             Max chi-squared p-value (default: 0.01)
      --ambient_rate_max              Max ambient rate (default: 0.5)
      --multiple_species_per_droplet  Allow multi-species droplets (default: true)

    Integration:
      --nFeatures       Number of integration features (default: 3000)

    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// ── Parameter validation ────────────────────────────────────────────────────

def validateParams() {
    def errors = []
    def warnings = []

    // Mandatory parameters
    if (!params.input)   errors << "ERROR: --input is required."
    if (!params.out)     errors << "ERROR: --out is required."

    // Validate input samplesheet
    if (params.input) {
        def input_file = file(params.input)
        if (!input_file.exists()) errors << "ERROR: Input samplesheet not found: ${params.input}"
    }

    // Validate output directory name
    if (params.out) {
        def out_basename = file(params.out).name
        if (!(out_basename ==~ /^[a-zA-Z0-9_\-]+$/)) {
            errors << "ERROR: Output directory name '${out_basename}' contains invalid characters. Use only letters, numbers, dashes and underscores."
        }
        if (file(params.out).exists()) {
            warnings << "WARNING: Output directory '${params.out}' already exists. Existing files may be overwritten."
        }
    }

    // Parse species list
    def species_list = params.species.toString().split(',').collect { it.trim() }
    def n_species = species_list.size()
    def available_species = params.species_map.keySet()

    for (sp in species_list) {
        if (!(sp in available_species)) {
            if (!params.genome || !params.gtf) {
                errors << "ERROR: Species '${sp}' is not available. Either provide an available species or supply --genome and --gtf parameters. Available species: ${available_species.join(', ')}"
            }
        }
    }

    // Determine clean method
    def clean_method = params.clean
    if (clean_method == null) {
        clean_method = n_species > 1 ? 'chi' : 'diem'
    }
    if (!(clean_method in ['diem', 'chi'])) {
        errors << "ERROR: --clean must be 'diem' or 'chi', got: '${clean_method}'"
    }
    if (n_species == 1 && clean_method == 'chi') {
        warnings << "WARNING: --clean 'chi' is not available for single-species. Using 'diem' instead."
        clean_method = 'diem'
    }

    // Validate parameter ranges
    if (params.min_unsplice_ratio < 0 || params.min_unsplice_ratio > 1) {
        errors << "ERROR: --min_unsplice_ratio must be between 0 and 1."
    }
    if (params.min_nCount_RNA < 1) errors << "ERROR: --min_nCount_RNA must be a positive integer."
    if (params.min_nFeature_RNA < 1) errors << "ERROR: --min_nFeature_RNA must be a positive integer."
    if (params.max_mt < 0 || params.max_mt > 100) errors << "ERROR: --max_mt must be between 0 and 100."
    if (params.max_cp < 0 || params.max_cp > 100) errors << "ERROR: --max_cp must be between 0 and 100."
    if (params.nHVG < 1) errors << "ERROR: --nHVG must be a positive integer."
    if (params.min_ncell_expr < 1) errors << "ERROR: --min_ncell_expr must be a positive integer."
    if (params.max_doublet_score < 0 || params.max_doublet_score > 1) {
        errors << "ERROR: --max_doublet_score must be between 0 and 1."
    }
    if (params.min_nClusterMarker < 1) errors << "ERROR: --min_nClusterMarker must be a positive integer."
    if (params.chisq_pvalues_max < 0 || params.chisq_pvalues_max > 1) {
        errors << "ERROR: --chisq_pvalues_max must be between 0 and 1."
    }
    if (params.ambient_rate_max < 0 || params.ambient_rate_max > 1) {
        errors << "ERROR: --ambient_rate_max must be between 0 and 1."
    }

    // Boolean parameter validation
    def rb = params.remove_doublet
    if (!(rb instanceof Boolean) && !(rb.toString().toLowerCase() in ['true', 'false'])) {
        errors << "ERROR: --remove_doublet must be a boolean value (true or false)."
    }
    def ms = params.multiple_species_per_droplet
    if (!(ms instanceof Boolean) && !(ms.toString().toLowerCase() in ['true', 'false'])) {
        errors << "ERROR: --multiple_species_per_droplet must be a boolean value (true or false)."
    }

    // Validate genome/gtf paths if provided
    if (params.genome) {
        params.genome.toString().split(',').each { g ->
            if (!file(g.trim()).exists()) errors << "ERROR: --genome path not found: ${g.trim()}"
        }
    }
    if (params.gtf) {
        params.gtf.toString().split(',').each { g ->
            if (!file(g.trim()).exists()) errors << "ERROR: --gtf path not found: ${g.trim()}"
        }
    }

    // Cross-parameter warnings
    if (clean_method == 'chi' && params.min_diem_debris_score != 1) {
        warnings << "WARNING: --min_diem_debris_score is ignored when --clean is 'chi'."
    }
    if (clean_method == 'diem') {
        if (params.min_UMI_per_cell_barcode != 400) {
            warnings << "WARNING: --min_UMI_per_cell_barcode is ignored when --clean is 'diem'."
        }
        if (params.chisq_pvalues_max != 0.01) {
            warnings << "WARNING: --chisq_pvalues_max is ignored when --clean is 'diem'."
        }
        if (params.ambient_rate_max != 0.5) {
            warnings << "WARNING: --ambient_rate_max is ignored when --clean is 'diem'."
        }
        if (params.multiple_species_per_droplet != true) {
            warnings << "WARNING: --multiple_species_per_droplet is ignored when --clean is 'diem'."
        }
    }

    // Validate ref_yaml
    if (params.ref_yaml && !file(params.ref_yaml).exists()) {
        errors << "ERROR: --ref_yaml file not found: ${params.ref_yaml}"
    }

    // Print warnings
    for (w in warnings) log.warn(w)

    // Fail on errors
    if (errors) {
        for (e in errors) log.error(e)
        exit 1
    }

    return [
        species_list: species_list,
        n_species:    n_species,
        clean_method: clean_method
    ]
}

// ── Parse samplesheet ───────────────────────────────────────────────────────

def parseSamplesheet(samplesheet_path) {
    def samples = []
    def lines = file(samplesheet_path).readLines()
    def header = lines[0].split(',').collect { it.trim() }

    // Validate header
    def required_cols = ['sample', 'fastq_dir', 'fastq_R1', 'fastq_R2']
    for (col in required_cols) {
        if (!(col in header)) {
            log.error "ERROR: Samplesheet missing required column: ${col}"
            exit 1
        }
    }

    def col_idx = [:]
    for (col in required_cols) {
        col_idx[col] = header.indexOf(col)
    }

    for (int i = 1; i < lines.size(); i++) {
        def line = lines[i].trim()
        if (line == '') continue
        def fields = line.split(',').collect { it.trim() }

        def sample_name = fields[col_idx['sample']]
        def fastq_dir   = fields[col_idx['fastq_dir']]
        def fastq_r1    = fields[col_idx['fastq_R1']]

        // Validate sample name
        if (!(sample_name ==~ /^[a-zA-Z0-9_\-]+$/)) {
            log.error "ERROR: Sample name '${sample_name}' (row ${i+1}) contains invalid characters. Use only letters, numbers, dashes and underscores."
            exit 1
        }

        // Validate fastq paths
        if (!file(fastq_dir).isDirectory()) {
            log.error "ERROR: FASTQ directory not found for sample '${sample_name}': ${fastq_dir}"
            exit 1
        }

        // Extract sample prefix from R1 filename (string before _S\d+_L\d+_R\d+_\d+.fastq.gz)
        // Use only the filename in case fastq_R1 is a full path
        def fastq_r1_name = file(fastq_r1).name
        def fastq_prefix = fastq_r1_name.replaceAll(/_S\d+_L\d+_R\d+_\d+\.fastq\.gz$/, '')
        if (fastq_prefix == fastq_r1_name) {
            log.warn "WARNING: Sample '${sample_name}' fastq_R1 filename '${fastq_r1_name}' does not match standard Illumina naming pattern (_S*_L*_R*_*.fastq.gz). Using full filename as prefix."
        }

        samples << [sample_name, fastq_dir, fastq_prefix]
    }

    return samples
}

// ── Resolve cellranger reference ────────────────────────────────────────────

def resolveReference(species_list) {
    def n_species = species_list.size()

    if (n_species == 1) {
        def sp = species_list[0]
        def ref_info = params.species_map[sp]

        if (ref_info && file(ref_info.cellranger).isDirectory()) {
            log.info "Using existing cellranger reference for ${sp}: ${ref_info.cellranger}"
            return [type: 'existing', ref_path: ref_info.cellranger]
        } else if (params.genome && params.gtf) {
            log.info "Will build cellranger reference for ${sp}"
            return [type: 'build_single', species: sp, genome: params.genome, gtf: params.gtf]
        } else if (ref_info) {
            log.info "Will build cellranger reference for ${sp} from species_map"
            return [type: 'build_single', species: sp, genome: ref_info.genome, gtf: ref_info.gtf]
        } else {
            log.error "Cannot resolve reference for species '${sp}'. Provide --genome and --gtf."
            exit 1
        }
    } else {
        def sorted_species = species_list.sort().join(',')
        if (params.combined_ref_map.containsKey(sorted_species) && file(params.combined_ref_map[sorted_species]).isDirectory()) {
            log.info "Using existing combined cellranger reference: ${params.combined_ref_map[sorted_species]}"
            return [type: 'existing', ref_path: params.combined_ref_map[sorted_species]]
        } else {
            log.info "Will build combined cellranger reference for: ${species_list.join(', ')}"
            return [type: 'build_multi', species_list: species_list]
        }
    }
}

// ── Resolve annotation references for a species ────────────────────────────

def resolveAnnotationRefs(species_name) {
    def sp_key = species_name.replaceAll('_', ' ')
    def ref_info = params.species_map.containsKey(sp_key) ? params.species_map[sp_key] : null
    return [
        markers:    ref_info?.markers ?: '',
        seurat_ref: ref_info?.seurat_ref ?: ''
    ]
}

// ═══════════════════════════════════════════════════════════════════════════
// SINGLE-SPECIES WORKFLOW
// ═══════════════════════════════════════════════════════════════════════════

workflow SINGLE_SPECIES_WF {
    take:
        sample_ch          // [ sample_id, fastq_dir, fastq_prefix ]
        transcriptome_ch   // path to cellranger reference
        species_name       // string

    main:
        def preprocess_dir = "${params.out}/preprocess"
        def species_dir    = species_name.replaceAll(' ', '_')
        def sp_info        = params.species_map[species_name]
        def gtf_path       = sp_info ? sp_info.gtf : params.gtf
        if (!gtf_path) {
            log.error "ERROR: Cannot determine GTF file for velocyto. Provide --gtf or add species '${species_name}' to species_map."
            exit 1
        }
        def gtf_file = file(gtf_path)

        // Step 4: CellRanger count
        CELLRANGER_COUNT(sample_ch, transcriptome_ch)

        // Step 5: Initial preprocess (DIEM only, no velocyto)
        // Pass absolute path strings (not file objects) since PREPROCESS_INITIAL uses val inputs
        def cr_outs_ch = CELLRANGER_COUNT.out.cellranger_out.map { sid, dir ->
            [sid, "${dir.toAbsolutePath()}/outs"]
        }
        PREPROCESS_INITIAL(cr_outs_ch, species_name, preprocess_dir)

        // Step 6: Separate reads and count in 2kb windows
        def bam_ch = CELLRANGER_COUNT.out.cellranger_out.map { sid, dir ->
            [sid, file("${dir}/outs/possorted_genome_bam.bam")]
        }
        SEPARATE_READS(bam_ch, preprocess_dir)

        def all_bams_ch = SEPARATE_READS.out.mapped_bam.map { sid, bam -> [sid, 'mapped', bam] }
            .mix(
                SEPARATE_READS.out.uniq_bam.map   { sid, bam -> [sid, 'uniq', bam] },
                SEPARATE_READS.out.multi_bam.map  { sid, bam -> [sid, 'multi', bam] }
            )
        COUNT_READS_2KB(all_bams_ch, preprocess_dir)

        // Step 7: Velocyto
        def velocyto_input_ch = CELLRANGER_COUNT.out.cellranger_out
            .join(PREPROCESS_INITIAL.out.barcodes)
            .map { sid, cr_dir, barcode_file -> [sid, cr_dir, barcode_file] }

        VELOCYTO_RUN(velocyto_input_ch, gtf_file, "${params.out}/cellranger")
        UNSPLICE_RATIO(VELOCYTO_RUN.out.loom, preprocess_dir)

        // Step 8: Full preprocess with velocyto
        // Pass absolute path strings since PREPROCESS_WITH_VELOCYTO uses val inputs
        def preprocess_velo_ch = CELLRANGER_COUNT.out.cellranger_out
            .join(UNSPLICE_RATIO.out.unsplice_txt)
            .map { sid, cr_dir, velo_txt -> [sid, "${cr_dir.toAbsolutePath()}/outs", velo_txt.toAbsolutePath().toString()] }

        PREPROCESS_WITH_VELOCYTO(preprocess_velo_ch, species_name, preprocess_dir)

        // Step 9: Summary report
        // Must run AFTER step 8 (PREPROCESS_WITH_VELOCYTO). Derive seur_sels_ch from its
        // output so Nextflow waits for all preprocessing before scheduling the report.
        def cr_dirs_ch = CELLRANGER_COUNT.out.cellranger_out
            .map { sid, dir -> dir.toAbsolutePath().toString() }
            .collect()
        def cr_sels_ch = CELLRANGER_COUNT.out.cellranger_out
            .map { sid, dir -> sid }
            .collect()
        def seur_sels_ch = PREPROCESS_WITH_VELOCYTO.out.seur_clean
            .map { sid, rds -> sid }
            .collect()

        SUMMARY_REPORT(
            cr_dirs_ch,
            cr_sels_ch,
            [file(preprocess_dir).toAbsolutePath().toString()],
            seur_sels_ch,
            species_name,
            preprocess_dir
        )

        // Step 10: Integration
        def rds_paths_ch = PREPROCESS_WITH_VELOCYTO.out.seur_clean
            .map { sid, rds -> rds.toAbsolutePath().toString() }
            .collect()

        def anno = resolveAnnotationRefs(species_name)
        INTEGRATION(rds_paths_ch, anno.markers, anno.seurat_ref, "${params.out}/integration")
}

// ═══════════════════════════════════════════════════════════════════════════
// MULTI-SPECIES WORKFLOW
// ═══════════════════════════════════════════════════════════════════════════

workflow MULTI_SPECIES_WF {
    take:
        sample_ch          // [ sample_id, fastq_dir, fastq_prefix ]
        transcriptome_ch   // path to cellranger reference
        species_list       // list of species names
        clean_method       // 'diem' or 'chi'

    main:
        def species_csv   = species_list.join(',')
        def preprocess_base = "${params.out}/preprocess"

        // Step 4: CellRanger count
        CELLRANGER_COUNT(sample_ch, transcriptome_ch)

        // Step 5: Demultiplex
        // Pass absolute path strings since DEMULTIPLEX and DEMUX_BAM use val inputs
        def cr_outs_ch = CELLRANGER_COUNT.out.cellranger_out.map { sid, dir ->
            [sid, "${dir.toAbsolutePath()}/outs"]
        }
        DEMULTIPLEX(cr_outs_ch, species_csv, species_list)
        DEMUX_BAM(cr_outs_ch, species_list)

        // If diem, run preprocess on demultiplexed matrices
        if (clean_method == 'diem') {
            // Create [sample, species_dir, matrix_dir] channel from demux outputs
            // preprocess.R expects the directory CONTAINING raw_feature_bc_matrix
            // (it appends "/raw_feature_bc_matrix" itself). Use the published absolute path.
            def demux_preprocess_input = DEMULTIPLEX.out.demux_matrices
                .flatMap { sid, dirs ->
                    def entries = []
                    for (sp in species_list) {
                        def sp_dir = sp.replaceAll(' ', '_')
                        def any_match = (dirs instanceof List ? dirs : [dirs]).any { it.toString().contains(sp_dir) }
                        if (any_match) {
                            // Pass the parent directory that contains raw_feature_bc_matrix
                            def matrix_parent = file("${params.out}/preprocess/${sp_dir}/${sid}").toAbsolutePath().toString()
                            entries << [sid, sp_dir, matrix_parent]
                        }
                    }
                    return entries
                }

            def publish_dirs = demux_preprocess_input.map { sid, sp, dir ->
                "${preprocess_base}/${sp}"
            }
            PREPROCESS_DEMUX(demux_preprocess_input, publish_dirs)
        }

        // Step 6: Separate and count reads per species
        // Flatten demuxed BAMs into [sample_id, species_dir, bam]
        def demux_bams_flat = DEMUX_BAM.out.demux_bams
            .flatMap { sid, bams ->
                def entries = []
                for (sp in species_list) {
                    def sp_dir = sp.replaceAll(' ', '_')
                    for (b in (bams instanceof List ? bams : [bams])) {
                        if (b.toString().contains(sp_dir)) {
                            entries << ["${sp_dir}_${sid}", b]
                        }
                    }
                }
                return entries
            }

        SEPARATE_READS_MULTI(demux_bams_flat, preprocess_base)

        def all_bams_multi = SEPARATE_READS_MULTI.out.mapped_bam.map { sid, bam -> [sid, 'mapped', bam] }
            .mix(
                SEPARATE_READS_MULTI.out.uniq_bam.map   { sid, bam -> [sid, 'uniq', bam] },
                SEPARATE_READS_MULTI.out.multi_bam.map  { sid, bam -> [sid, 'multi', bam] }
            )
        COUNT_READS_2KB_MULTI(all_bams_multi, preprocess_base)

        // Step 7: Velocyto per species
        // Build a single channel: [sp_dir, sample_id, bam, bam_index, barcodes, gtf]
        // VELOCYTO_RUN_MULTI uses direct BAM files (no outs/ subdirectory).
        def multi_velo_ch = Channel.empty()
        for (sp in species_list) {
            def sp_dir  = sp.replaceAll(' ', '_')
            def sp_info = params.species_map[sp]
            def sp_gtf  = file(sp_info.gtf)

            def sp_bam_ch = DEMUX_BAM.out.demux_bams
                .flatMap { sid, bams ->
                    (bams instanceof List ? bams : [bams])
                        .findAll { it.toString().contains("${sp_dir}/") }
                        .collect { [sid, it] }
                }
            def sp_bai_ch = DEMUX_BAM.out.demux_bam_indices
                .flatMap { sid, bais ->
                    (bais instanceof List ? bais : [bais])
                        .findAll { it.toString().contains("${sp_dir}/") }
                        .collect { [sid, it] }
                }

            def barcodes_ch
            if (clean_method == 'diem') {
                barcodes_ch = PREPROCESS_DEMUX.out.barcodes
                    .filter { sid, sp_name, bc -> sp_name == sp_dir }
                    .map    { sid, sp_name, bc -> [sid, bc] }
            } else {
                barcodes_ch = DEMULTIPLEX.out.barcodes
                    .flatMap { sid, bcs ->
                        (bcs instanceof List ? bcs : [bcs])
                            .findAll { it.toString().contains("${sp_dir}/") }
                            .collect { [sid, it] }
                    }
            }

            def sp_velo_ch = sp_bam_ch
                .join(sp_bai_ch)
                .join(barcodes_ch)
                .map { sid, bam, bai, bc -> [sp_dir, sid, bam, bai, bc, sp_gtf] }

            multi_velo_ch = multi_velo_ch.mix(sp_velo_ch)
        }

        VELOCYTO_RUN_MULTI(multi_velo_ch, preprocess_base)
        UNSPLICE_RATIO_MULTI(VELOCYTO_RUN_MULTI.out.loom, preprocess_base)

        // Step 8: Full preprocess with velocyto for each species/sample
        // matrix_path: published dir containing raw_feature_bc_matrix (preprocess.R appends it)
        // velo_txt:    published unsplice ratio txt (preprocess/sp_dir/sample_id.txt)
        def preprocess_multi_ch = UNSPLICE_RATIO_MULTI.out.unsplice_txt
            .map { sp_dir, sid, txt ->
                def sp = species_list.find { it.replaceAll(' ', '_') == sp_dir }
                def matrix_path = file("${preprocess_base}/${sp_dir}/${sid}").toAbsolutePath().toString()
                def velo_path   = file("${preprocess_base}/${sp_dir}/${sid}.txt").toAbsolutePath().toString()
                [sp_dir, sid, sp, matrix_path, velo_path]
            }
        PREPROCESS_WITH_VELOCYTO_FOR_SPECIES(preprocess_multi_ch, preprocess_base)

        // Steps 9-10: Summary and integration per species
        def cr_dirs_ch = CELLRANGER_COUNT.out.cellranger_out
            .map { sid, dir -> dir.toAbsolutePath().toString() }
            .collect()
        def cr_sels_ch = CELLRANGER_COUNT.out.cellranger_out
            .map { sid, dir -> sid }
            .collect()

        for (sp in species_list) {
            def sp_dir         = sp.replaceAll(' ', '_')
            def sp_preprocess  = "${preprocess_base}/${sp_dir}"

            // Step 9: Summary report for this species.
            // seur_sels_ch is derived from PREPROCESS_WITH_VELOCYTO_FOR_SPECIES outputs
            // so Nextflow waits for all preprocessing (steps 5-8) before running the report.
            def sp_seur_sels_ch = PREPROCESS_WITH_VELOCYTO_FOR_SPECIES.out.seur_clean
                .filter { sp_d, sid, rds -> sp_d == sp_dir }
                .map    { sp_d, sid, rds -> sid }
                .collect()

            SUMMARY_REPORT_MULTI(
                cr_dirs_ch,
                cr_sels_ch,
                [file(sp_preprocess).toAbsolutePath().toString()],
                sp_seur_sels_ch,
                sp,
                sp_preprocess
            )

            // Step 10: Integration for this species
            def sp_rds_ch = PREPROCESS_WITH_VELOCYTO_FOR_SPECIES.out.seur_clean
                .filter { sp_d, sid, rds -> sp_d == sp_dir }
                .map    { sp_d, sid, rds -> rds.toAbsolutePath().toString() }
                .collect()

            def anno = resolveAnnotationRefs(sp)
            INTEGRATION(sp_rds_ch, anno.markers, anno.seurat_ref, "${params.out}/integration/${sp_dir}")
        }
}

// ── Main workflow ───────────────────────────────────────────────────────────

workflow {

    // Step 1: Validate parameters
    def validated   = validateParams()
    def species_list = validated.species_list
    def n_species    = validated.n_species
    def clean_method = validated.clean_method
    def is_multi     = n_species > 1

    log.info """
    =========================================================
     snRNAseq Processing Pipeline  v${workflow.manifest.version}
    =========================================================
     Input        : ${params.input}
     Output       : ${params.out}
     Species      : ${species_list.join(', ')}
     Mode         : ${is_multi ? 'Multi-species' : 'Single-species'}
     Clean method : ${clean_method}
    =========================================================
    """.stripIndent()

    // Step 2: Parse samplesheet and create channel
    def samples = parseSamplesheet(params.input)
    def sample_ch = Channel.from(samples)

    // Step 3: Resolve reference
    def ref_info = resolveReference(species_list)
    def transcriptome_ch

    if (ref_info.type == 'existing') {
        transcriptome_ch = Channel.value(file(ref_info.ref_path))
    } else if (ref_info.type == 'build_single') {
        def gtf_path = ref_info.gtf
        if (gtf_path.endsWith('.gff3') || gtf_path.endsWith('.gff')) {
            GFF_TO_GTF(Channel.value([ref_info.species, file(gtf_path)]))
            CELLRANGER_MKGTF(GFF_TO_GTF.out.gtf)
        } else {
            CELLRANGER_MKGTF(Channel.value([ref_info.species, file(gtf_path)]))
        }
        def species_dir_name = ref_info.species.replaceAll(' ', '_')
        CELLRANGER_MKREF(
            CELLRANGER_MKGTF.out.filtered_gtf.map { sp, gtf -> [sp, file(ref_info.genome), gtf] },
            species_dir_name
        )
        transcriptome_ch = CELLRANGER_MKREF.out.cellranger_ref.map { sp, ref -> ref }
    } else {
        // build_multi
        def genome_args = species_list.collect { sp ->
            def sp_info = params.species_map[sp]
            def genome_name = sp.replaceAll(' ', '_')
            "--genome=${genome_name} --fasta=${sp_info.genome} --genes=${sp_info.gtf}"
        }
        def species_dir_name = species_list.collect { it.replaceAll(' ', '_') }.join('_')
        CELLRANGER_MKREF_MULTI(genome_args, species_dir_name)
        transcriptome_ch = CELLRANGER_MKREF_MULTI.out.cellranger_ref
    }

    // Route to single-species or multi-species workflow
    if (!is_multi) {
        SINGLE_SPECIES_WF(sample_ch, transcriptome_ch, species_list[0])
    } else {
        MULTI_SPECIES_WF(sample_ch, transcriptome_ch, species_list, clean_method)
    }
}
