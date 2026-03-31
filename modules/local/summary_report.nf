process SUMMARY_REPORT {
    label 'process_high'
    label 'scrnaseq'
    publishDir "${publish_dir}", mode: 'copy'

    input:
    tuple val(cellranger_dirs), val(cellranger_sel), val(seur_dirs), val(seur_sel), val(species_name), val(publish_dir)

    output:
    path "summary.html", emit: summary_html

    script:
    cr_dirs_r = cellranger_dirs.collect { "'" + it + "'" }.join(', ')
    cr_sel_r  = cellranger_sel.collect  { "'" + it + "'" }.join(', ')
    sr_dirs_r = seur_dirs.collect       { "'" + it + "'" }.join(', ')
    sr_sel_r  = seur_sel.collect        { "'" + it + "'" }.join(', ')
    """
    mkdir -p '${publish_dir}'
    Rscript -e "rmarkdown::render(
      input = '${projectDir}/bin/summary.Rmd',
      params = list(
        title = 'scRNA-seq summary report',
        author = 'Qiaoshan Lin',
        cellranger_dir = c(${cr_dirs_r}),
        cellranger_sel = c(${cr_sel_r}),
        seur_dir = c(${sr_dirs_r}),
        seur_sel = c(${sr_sel_r}),
        species = '${species_name}'
      ),
      output_file = '${publish_dir}/summary.html'
    )"
    cp '${publish_dir}/summary.html' summary.html
    """
}
