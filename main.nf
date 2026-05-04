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

// ── Module includes (reference preparation only) ─────────────────────────────

include { CELLRANGER_MKREF }      from './modules/local/cellranger_mkref'
include { CELLRANGER_MKREF_MULTI } from './modules/local/cellranger_mkref_multi'

// ── Workflow includes ────────────────────────────────────────────────────────

include { SINGLE_SPECIES_WF } from './workflows/single_species'
include { MULTI_SPECIES_WF }  from './workflows/multi_species'

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

    Read counting:
      --bin_size        Genomic window size (bp) for bamCoverage read counting (default: 2000)

    Cleanup:
      --cleanup         Remove large intermediate files (*.bam, *.bai, *.loom) from
                        the output directory after successful pipeline completion (default: false)

    """.stripIndent()
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
    def species_list = params.species.toString().split(',').collect { token -> token.trim() }
    def n_species = species_list.size()
    def available_species = params.species_map.keySet()

    species_list.each { sp ->
        if (!(sp in available_species)) {
            errors << "ERROR: Species '${sp}' is not available. Please provide available species names in --species. Available species: ${available_species.join(', ')}"
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
    if (!(params.remove_doublet.toString().toLowerCase() in ['true', 'false'])) {
        errors << "ERROR: --remove_doublet must be a boolean value (true or false)."
    }
    if (!(params.multiple_species_per_droplet.toString().toLowerCase() in ['true', 'false'])) {
        errors << "ERROR: --multiple_species_per_droplet must be a boolean value (true or false)."
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

    // Validate ref_yaml. NB: do NOT try to mutate params.ref_yaml here —
    // params is read-only at runtime in modern Nextflow and the assignment
    // would silently no-op. Modules absolutize via `${file(params.ref_yaml)}`
    // at script-render time instead, which resolves relative paths against
    // the launch dir.
    if (params.ref_yaml && !file(params.ref_yaml).exists()) {
        errors << "ERROR: --ref_yaml file not found: ${params.ref_yaml}"
    }

    // Print warnings
    warnings.each { warning_msg -> log.warn(warning_msg) }

    // Fail on errors
    if (errors) {
        errors.each { error_msg -> log.error(error_msg) }
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
    def header = lines[0].split(',').collect { token -> token.trim() }

    // Validate header
    def required_cols = ['sample', 'fastq_dir', 'fastq_R1', 'fastq_R2']
    required_cols.each { col ->
        if (!(col in header)) {
            log.error "ERROR: Samplesheet missing required column: ${col}"
            exit 1
        }
    }

    def col_idx = [:]
    required_cols.each { col ->
        col_idx[col] = header.indexOf(col)
    }

    (1..<lines.size()).each { i ->
        def line = lines[i].trim()
        if (line == '') {
            return
        }
        def fields = line.split(',').collect { token -> token.trim() }

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

        if (!ref_info) {
            log.error "Cannot resolve reference for species '${sp}'. Please add species '${sp}' to params.species_map in nextflow.config."
            exit 1
        }

        // Prefer a pre-built cellranger reference if one is configured AND exists on disk
        if (ref_info.cellranger && file(ref_info.cellranger).isDirectory()) {
            log.info "Using existing cellranger reference for ${sp}: ${ref_info.cellranger}"
            return [type: 'existing', ref_path: ref_info.cellranger, gtf: ref_info.gtf]
        }

        // Otherwise build one — require genome + gtf to be present and to exist
        if (!ref_info.genome || !ref_info.gtf) {
            log.error "species_map['${sp}'] is missing required keys. Need 'genome' and 'gtf' (got: ${ref_info})."
            exit 1
        }
        if (!file(ref_info.genome).exists()) {
            log.error "Genome FASTA for species '${sp}' not found: ${ref_info.genome}"
            exit 1
        }
        if (!file(ref_info.gtf).exists()) {
            log.error "GTF for species '${sp}' not found: ${ref_info.gtf}"
            exit 1
        }
        log.info "Will build cellranger reference for ${sp} from species_map"
        return [type: 'build_single', species: sp, genome: ref_info.genome, gtf: ref_info.gtf]

    } else {
        def species_list_str = species_list.sort().join(',')
        def combined_ref = params.combined_ref_map[species_list_str]

        // Prefer a pre-built combined cellranger reference if one is configured AND exists on disk
        if (combined_ref && file(combined_ref).isDirectory()) {
            log.info "Using existing combined cellranger reference: ${combined_ref}"
            return [type: 'existing', ref_path: combined_ref]
        }

        // Otherwise build one — every species needs a complete species_map entry
        def missing = species_list.findAll { sp -> !params.species_map.containsKey(sp) }
        if (missing) {
            log.error "Cannot resolve reference for species '${species_list_str}'. Missing species in species_map: ${missing.join(', ')}"
            log.error "Available species: ${params.species_map.keySet().join(', ')}"
            exit 1
        }
        def incomplete = species_list.findAll { sp ->
            def info = params.species_map[sp]
            !info?.genome || !info?.gtf
        }
        if (incomplete) {
            log.error "species_map entries for [${incomplete.join(', ')}] must have both 'genome' and 'gtf' keys to build a combined reference."
            exit 1
        }
        log.info "Will build combined cellranger reference for: ${species_list_str}"
        return [type: 'build_multi', species_list: species_list]
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// MAIN WORKFLOW
// ═══════════════════════════════════════════════════════════════════════════

workflow {

    // Handle help request inside workflow scope (DSL2-safe).
    if (params.help) {
        helpMessage()
        exit 0
    }

    // Step 1: Validate parameters
    def validated    = validateParams()
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
    def samples   = parseSamplesheet(params.input)
    def sample_ch = channel.from(samples)

    // Step 3: Resolve reference + per-species GTF map
    def ref_info = resolveReference(species_list) 
    // [type: 'existing', ref_path: '/project/gzy8899/qiaoshan/scRNAseq/refs/Athaliana_Crubella/cellranger'] or 
    // [type: 'build_single', species: 'Arabidopsis thaliana', genome: '/project/gzy8899/.../Athaliana_447_TAIR10.fa', gtf: '/project/gzy8899/.../Athaliana_447_TAIR10.gene.gtf'] or 
    // [type: 'build_multi', species_list: ['Arabidopsis thaliana', 'Capsella rubella']]
    def cellranger_index_ch
    def gtf_ch

    if (ref_info.type == 'existing') {
        cellranger_index_ch = channel.value(file(ref_info.ref_path))
        if (!is_multi) {
            gtf_ch = channel.value(file(ref_info.gtf))
        }
    } else if (ref_info.type == 'build_single') {
        def species_dir_name = ref_info.species.replaceAll(' ', '_')
        CELLRANGER_MKREF(
            channel.value([ref_info.species, file(ref_info.genome), file(ref_info.gtf)]),
            species_dir_name
        )
        cellranger_index_ch = CELLRANGER_MKREF.out.cellranger_ref.map { _sp, ref -> ref }
        gtf_ch = channel.value(file(ref_info.gtf))
    } else if (ref_info.type == 'build_multi') {
        def cellranger_args = species_list.collect { sp ->
            def sp_info     = params.species_map[sp]
            def genome_path = file(sp_info.genome).toAbsolutePath()
            def gtf_path    = file(sp_info.gtf).toAbsolutePath()
            def genome_name = sp.replaceAll(' ', '_')
            "--genome=${genome_name} --fasta=${genome_path} --genes=${gtf_path}"
        }
        def cellranger_args_ch = channel.value(cellranger_args)

        def species_dir_name = species_list.collect { sp -> sp.replaceAll(' ', '_') }.join('_')
        CELLRANGER_MKREF_MULTI(cellranger_args_ch, species_dir_name)
        cellranger_index_ch = CELLRANGER_MKREF_MULTI.out.cellranger_ref
    }

    // Step 4+: Route to single-species or multi-species workflow
    if (!is_multi) {
        def sp   = species_list[0]
        SINGLE_SPECIES_WF(
            sample_ch,
            cellranger_index_ch,
            sp,
            gtf_ch
        )
    } else {
        // Per-species GTF lookup (species name -> GTF file path) for downstream
        // velocyto, etc. Plain Groovy Map (not a channel) so the workflow can
        // index by species name directly.
        def gtf_map = species_list.collectEntries { sp ->
            [(sp): file(params.species_map[sp].gtf)]
        }
        MULTI_SPECIES_WF(
            sample_ch,
            cellranger_index_ch,
            species_list,
            clean_method,
            gtf_map
        )
    }

    // ── Post-run cleanup ────────────────────────────────────────────────────
    // Registered inside the entry workflow (required by Nextflow's strict
    // syntax / v2 parser, which forbids top-level statements) and via the
    // assignment form `workflow.onComplete = { ... }` rather than the method-
    // call form `workflow.onComplete { ... }`. Assignment was also the
    // official workaround for the historical 'workflow/params null in
    // onComplete' bug (nextflow-io/nextflow#5261, fixed in 24.10), so this
    // spelling is portable across all supported Nextflow versions.
    // Wrapped in try/catch so cleanup can never fail the run reporting.
    workflow.onComplete = {
        try {
            if (params.cleanup && workflow.success) {
                log.info "Removing large intermediate files (*.bam, *.bai, *.loom) under ${params.out} ..."
                ['bash', '-c', "find '${params.out}' \\( -name '*.bam' -o -name '*.bai' -o -name '*.loom' \\) -delete"].execute().waitFor()
                log.info "Cleanup complete."
            }
        } catch (Throwable t) {
            log.warn "onComplete cleanup encountered an error (non-fatal): ${t.message}"
        }
    }
}
