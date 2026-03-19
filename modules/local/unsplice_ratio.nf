process UNSPLICE_RATIO {
    tag "${sample_id}"
    label 'process_low'
    label 'velocyto'
    publishDir "${publish_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(loom_file), val(publish_dir)

    output:
    tuple val(sample_id), path("${sample_id}.txt"), emit: unsplice_txt

    script:
    """
    python ${projectDir}/bin/unsplice_ratio.py ${loom_file} ${sample_id}.txt
    """
}
