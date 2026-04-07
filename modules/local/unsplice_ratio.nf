process UNSPLICE_RATIO {
    tag "${sample_id}"
    label 'process_low'
    label 'velocyto'
    publishDir "${publish_dir}", mode: 'copy', enabled: false

    input:
    tuple val(key), path(loom_file), val(sample_id), val(publish_dir)

    output:
    tuple val(key), path("${sample_id}.txt"), emit: unsplice_txt

    script:
    """
    python ${projectDir}/bin/unsplice_ratio.py ${loom_file} ${sample_id}.txt
    """
}
