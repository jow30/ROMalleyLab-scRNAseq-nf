process SEPARATE_READS {
    tag "${sample_id}"
    label 'process_medium'
    label 'samtools'

    input:
    tuple val(key), path(cr_dir), val(sample_id), val(publish_dir)

    output:
    tuple val(key), path("${sample_id}.mapped.bam"), emit: mapped_bam
    tuple val(key), path("${sample_id}.uniq.bam"),   emit: uniq_bam
    tuple val(key), path("${sample_id}.multi.bam"),  emit: multi_bam

    script:
    """
    samtools view -b -F 4 -@ ${task.cpus} ${cr_dir}/outs/possorted_genome_bam.bam > ${sample_id}.mapped.bam
    samtools view -b -q 255 -@ ${task.cpus} ${sample_id}.mapped.bam \
        -o ${sample_id}.uniq.bam -U ${sample_id}.multi.bam
    """
}
