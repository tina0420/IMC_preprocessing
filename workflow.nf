#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input  = "data/biof501_data.h5ad"
params.outdir = "figures"

process BATCH_EFFECT_CHECK {
    tag "batch_effect"

    input:
    path input_file

    output:
    path "batch_umap.png"

    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    python ${baseDir}/batch_effect_check_plot.py \
        --input ${input_file}/ --figdir .
    """
}


process B_GMM_CLUSTER {
    tag "clustering"

    input:
    path input_file

    output:
    path "new_adata.h5ad"

    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    python ${baseDir}/run_leiden_clustering.py \
        --input ${input_file}/ --output .
    """
}

process PLOT_HEATMAP {
    tag "heatmap"

    input:
    path input_file

    output:
    path "cluster_result.png"

    publishDir "${params.outdir}", mode: 'copy'

    script:
    """
    python ${baseDir}/plot_cluster_heatmap.py \
        --input ${input_file} --figdir .
    """
}

workflow {
    input_ch = Channel.fromPath(params.input)

    umap      = BATCH_EFFECT_CHECK(input_ch)
    clustered = B_GMM_CLUSTER(input_ch)
    PLOT_HEATMAP(clustered)
}
