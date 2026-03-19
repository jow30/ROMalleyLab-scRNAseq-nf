process SEPARATE_READS {
    tag "${sample_id}"
    label 'process_medium'
    label 'samtools'
    publishDir "${publish_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(cr_dir), val(publish_dir)

    output:
    tuple val(sample_id), path("${sample_id}.mapped.bam"), emit: mapped_bam
    tuple val(sample_id), path("${sample_id}.uniq.bam"),   emit: uniq_bam
    tuple val(sample_id), path("${sample_id}.multi.bam"),  emit: multi_bam

    script:
    """
    samtools view -b -F 4 -@ ${task.cpus} ${cr_dir}/outs/possorted_genome_bam.bam > ${sample_id}.mapped.bam
    samtools view -b -q 255 -@ ${task.cpus} ${sample_id}.mapped.bam \
        -o ${sample_id}.uniq.bam -U ${sample_id}.multi.bam
    """
}
