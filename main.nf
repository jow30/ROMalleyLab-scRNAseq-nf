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

include { GFF_TO_GTF }            from './modules/local/gff_to_gtf'
include { CELLRANGER_MKGTF }      from './modules/local/cellranger_mkgtf'
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

    Read counting:
      --bin_size        Genomic window size (bp) for bamCoverage read counting (default: 2000)

    Integration:
      --nFeatures       Number of integration features (default: 3000)

    Cleanup:
      --cleanup         Remove large intermediate files (*.bam, *.bai, *.loom) from
                        the output directory after successful pipeline completion (default: false)

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

// ── Resolve per-species annotation, GTF, and integration references ─────────

def resolveSpeciesInfo(species_name) {
    def sp_key   = species_name.replaceAll('_', ' ')
    def ref_info = params.species_map.containsKey(sp_key) ? params.species_map[sp_key] : null
    return [
        gtf:        ref_info?.gtf        ?: (params.gtf ?: ''),
        markers:    ref_info?.markers    ?: '',
        seurat_ref: ref_info?.seurat_ref ?: ''
    ]
}

// ═══════════════════════════════════════════════════════════════════════════
// MAIN WORKFLOW
// ═══════════════════════════════════════════════════════════════════════════

workflow {

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

    // Step 4+: Route to single-species or multi-species workflow
    if (!is_multi) {
        def sp   = species_list[0]
        def info = resolveSpeciesInfo(sp)
        if (!info.gtf) {
            log.error "ERROR: Cannot determine GTF file for velocyto. Provide --gtf or add species '${sp}' to species_map."
            exit 1
        }
        SINGLE_SPECIES_WF(
            sample_ch,
            transcriptome_ch,
            sp,
            file(info.gtf),
            info.markers,
            info.seurat_ref
        )
    } else {
        def anno_map = species_list.collectEntries { sp -> [sp, resolveSpeciesInfo(sp)] }
        MULTI_SPECIES_WF(
            sample_ch,
            transcriptome_ch,
            species_list,
            clean_method,
            anno_map
        )
    }
}

// ── Post-run cleanup ─────────────────────────────────────────────────────────

workflow.onComplete {
    if (params.cleanup && workflow.success) {
        log.info "Removing large intermediate files (*.bam, *.bai, *.loom) under ${params.out} ..."
        ['*.bam', '*.bai', '*.loom'].each { pattern ->
            ['bash', '-c', "find '${params.out}' -name '${pattern}' -delete"].execute().waitFor()
        }
        log.info "Cleanup complete."
    }
}
