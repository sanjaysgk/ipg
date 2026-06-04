process PEPQUERY_ANNOTATE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-bd0ce31f08f73d42bb0edebd8f5c7dcd9b1eefb0:1a671ad7f1b9f7bae925b3fe7baf49e06e59e70e-0' :
        'biocontainers/mulled-v2-bd0ce31f08f73d42bb0edebd8f5c7dcd9b1eefb0:1a671ad7f1b9f7bae925b3fe7baf49e06e59e70e-0' }"

    input:
    tuple val(meta), path(peptides_tsv), path(psm_rank)

    output:
    tuple val(meta), path("integrated_peptides_validated.tsv"), emit: peptides
    path "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'pepquery_annotate.py'

    stub:
    """
    touch integrated_peptides_validated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
