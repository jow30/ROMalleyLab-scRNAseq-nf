/*
 * Read counting module: separate read types and count in 2kb windows
 */

process SEPARATE_READS {
    label 'process_medium'
    label 'samtools'
    publishDir "${publish_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file)
    val publish_dir

    output:
    tuple val(sample_id), path("${sample_id}.mapped.bam"), emit: mapped_bam
    tuple val(sample_id), path("${sample_id}.uniq.bam"),   emit: uniq_bam
    tuple val(sample_id), path("${sample_id}.multi.bam"),  emit: multi_bam

    script:
    """
    samtools view -b -F 4 -@ ${task.cpus} ${bam_file} > ${sample_id}.mapped.bam
    samtools view -b -q 255 -@ ${task.cpus} ${sample_id}.mapped.bam \
        -o ${sample_id}.uniq.bam -U ${sample_id}.multi.bam
    """
}

process COUNT_READS_2KB {
    label 'process_low'
    label 'deeptools'
    publishDir "${publish_dir}", mode: 'copy'

    input:
    tuple val(sample_id), val(read_type), path(bam_file)
    val publish_dir

    output:
    tuple val(sample_id), val(read_type), path("${sample_id}.${read_type}.tsv"), emit: count_table

    script:
    """
    samtools index -@ ${task.cpus} ${bam_file}

    bamCoverage \\
        -b ${bam_file} \\
        -o ${sample_id}.${read_type}.bedgraph \\
        --binSize 2000 \\
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
