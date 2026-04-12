process DEMUX_BAM {
    tag "${sample_id}"
    label 'process_medium'
    label 'samtools'
    // Publish only BAM + index. Do not publish the whole cellranger/<sample>/ tree:
    // DEMULTIPLEX already publishes outs/raw_feature_bc_matrix into the same path;
    // copying the entire directory here would replace outs/ and drop the matrices.
    publishDir "${params.out}/demultiplex", mode: 'copy', saveAs: { fn ->
        (fn.endsWith('.bam') || fn.endsWith('.bai')) ? fn : null
    }

    input:
    tuple val(sample_id), val(cellranger_outs_path)
    val species_list

    output:
    // Per-species sample directory: <sp_dir>/cellranger/<sample_id>/
    // Contains outs/possorted_genome_bam.bam — mirrors CellRanger layout so
    // VELOCYTO_RUN can read ${cellranger_dir}/outs/possorted_genome_bam.bam
    tuple val(sample_id), path("*/cellranger/${sample_id}"), emit: demux_sample_dirs

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
