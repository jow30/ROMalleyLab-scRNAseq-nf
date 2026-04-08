process COUNT_READS {
    tag "${sample_id} (${read_type})"
    label 'process_low'
    label 'deeptools'
    publishDir "${publish_dir}/read_counts_in_genomic_bins", mode: 'copy'

    input:
    tuple val(key), val(read_type), path(bam_file), val(sample_id), val(publish_dir)

    output:
    tuple val(key), val(read_type), path("${sample_id}.${read_type}.tsv"), emit: count_table

    script:
    """
    samtools index -@ ${task.cpus} ${bam_file}

    bamCoverage \\
        -b ${bam_file} \\
        -o ${sample_id}.${read_type}.bedgraph \\
        --binSize ${params.bin_size} \\
        --normalizeUsing None \\
        --outFileFormat bedgraph \\
        -p ${task.cpus}

    awk 'BEGIN {OFS="\\t"}
         NR==FNR {total += \$4; next}
         FNR==1  {print "chr", "start", "end", "count", "ratio"}
                 {print \$1, \$2, \$3, \$4, \$4/total}
    ' ${sample_id}.${read_type}.bedgraph ${sample_id}.${read_type}.bedgraph > ${sample_id}.${read_type}.tsv

    echo "Done. Total bins: \$(tail -n +2 ${sample_id}.${read_type}.tsv | wc -l)"

    rm -f ${sample_id}.${read_type}.bedgraph
    """
}
