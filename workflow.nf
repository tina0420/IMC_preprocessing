#!/usr/bin/env nextflow

// parameter setting
params.input_adata = "$projectDir/data/biof501_data.h5ad"
params.outputDir = "$projectDir/output"


process batch_effect_check_plot {
    tag "Batch effect check plot"

    input:
        path input_adata
        val figdir

    output:
        path "${figdir}/batch_umap.png"
    script:
    """
    python batch_effect_check_plot.py \
    --input ${input_adata} \
    --figdir ${figdir}
    """
}


process run_leiden_clustering {
    tag "Run leiden clustering"

    input:
        path input_adata
        val outputDir

    output:
        path "${outputDir}/new_adata.h5ad", emit: processed_h5ad

    script:
    """
    run_leiden_clustering.py \
    --input ${input_adata} \
    --output ${outputDir}
    """
}


process plot_cluster_heatmap {
    tag "Plot cluster heatmap"

    input:
        path input_adata
        val figdir

    output:
        path "${figdir}/cluster_result.png",  emit: heatmap_png

    script:
    """
    plot_cluster_heatmap.py \
    --input ${input_adata} \
    --output ${outputDir}
    """
}


workflow {
    batch_png = batch_effect_check_plot(
        params.input_adata,
        params.outputDir
    )
    // batch_effect_check_plot(params.input_adata, params.outputDir)

    // new_adata = run_leiden_clustering(params.input_adata, params.outputDir)
    processed_h5ad = run_leiden_clustering(
        params.input_adata,
        params.outputDir
    )

    // plot_cluster_heatmap(new_adata, params.outputDir)
    
    heatmap_png = plot_cluster_heatmap(
        processed_h5ad,
        params.outputDir
    )

}


workflow.onComplete {
    log.info ( workflow.success ? "\nDone! The ouput is stored in --> $projectDir/$params.outputDir\n" : "Error in the workflow!" )
}

