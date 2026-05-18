process IMMUNOINFORMATICS_REPORT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-bd0ce31f08f73d42bb0edebd8f5c7dcd9b1eefb0:1a671ad7f1b9f7bae925b3fe7baf49e06e59e70e-0' :
        'biocontainers/mulled-v2-bd0ce31f08f73d42bb0edebd8f5c7dcd9b1eefb0:1a671ad7f1b9f7bae925b3fe7baf49e06e59e70e-0' }"

    input:
    tuple val(meta),
        path(peptides),
        path(netmhcpan,   stageAs: 'netmhcpan.tsv'),
        path(netmhciipan, stageAs: 'netmhciipan.tsv'),
        path(gibbs,       stageAs: 'gibbs.tsv'),
        path(flashlfq,    stageAs: 'flashlfq.tsv'),
        path(blastp,      stageAs: 'blastp.tsv')

    output:
    tuple val(meta), path("${meta.id}_immunoinformatics_report.html"), emit: html
    path "versions.yml",                                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Optional inputs arrive as the NO_FILE sentinel when the upstream
    // module did not run. The template handles that check internally.
    template 'generate_report.py'

    stub:
    """
    touch ${meta.id}_immunoinformatics_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
