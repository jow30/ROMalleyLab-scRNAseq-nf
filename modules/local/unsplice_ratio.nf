process UNSPLICE_RATIO {
    tag "${sample_id}"
    label 'process_low'
    label 'velocyto'

    input:
    tuple val(key), path(loom_file), val(sample_id), val(publish_dir)

    output:
    tuple val(key), path("${sample_id}.txt"), emit: unsplice_txt

    script:
    def script_hash = file("${projectDir}/bin/unsplice_ratio.py").text.md5()
    """
    # script_hash: ${script_hash}
    python ${projectDir}/bin/unsplice_ratio.py ${loom_file} ${sample_id}.txt
    """
}
