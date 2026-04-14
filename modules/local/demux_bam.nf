process DEMUX_BAM {
    tag "${sample_id}"
    label 'process_medium'
    label 'samtools'
    // Publish BAM + BAI as explicit path outputs (not the whole <sp>/cellranger/<sample>/
    // directory). Nextflow often does not apply publishDir `pattern` / `saveAs` to files
    // nested inside a directory output, so those approaches published nothing. File outputs
    // copy into demultiplex/ without replacing sibling outs/raw_feature_bc_matrix from DEMULTIPLEX.
    publishDir "${params.out}/demultiplex", mode: 'copy'

    input:
    tuple val(sample_id), val(cellranger_outs_path)
    val species_list

    output:
    // One path per species; workflow derives <sp>/cellranger/<sample_id> as bam.parent.parent
    tuple val(sample_id), path("*/cellranger/${sample_id}/outs/possorted_genome_bam.bam"), emit: demux_sample_dirs
    path("*/cellranger/${sample_id}/outs/possorted_genome_bam.bam.bai"), emit: demux_bai

    script:
    def species_cmds = species_list.collect { sp ->
        def sp_dir = sp.replaceAll(' ', '_')
        """echo "Start demultiplexing BAM for ${sp}"
mkdir -p ${sp_dir}/cellranger/${sample_id}/outs
INPUT_BAM=${cellranger_outs_path}/possorted_genome_bam.bam
OUTPUT_BAM=${sp_dir}/cellranger/${sample_id}/outs/possorted_genome_bam.bam
samtools view -h -@ ${task.cpus} "\${INPUT_BAM}" | \\
awk -v re="(${sp_dir})_+" '
    BEGIN { OFS="\\t" }
    /^@/ {
        if (/^@SQ/) {
            if (\$0 ~ re) {
                gsub(re, "", \$0)
                print
            }
        } else {
            gsub(re, "", \$0)
            print
        }
        next
    }
    {
        if(\$3 ~ "^" re) {
            for (i = 1; i <= NF; i++) { gsub(re, "", \$i) }
            print
        }
    }
' | samtools view -@ ${task.cpus} -b -o "\${OUTPUT_BAM}" -
samtools index -@ ${task.cpus} "\${OUTPUT_BAM}"
echo "Finish demultiplexing BAM for ${sp}."
"""
    }.join('\n')
    """
    ${species_cmds}
    """
}
