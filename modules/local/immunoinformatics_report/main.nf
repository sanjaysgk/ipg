process IMMUNOINFORMATICS_REPORT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-77320db00eefbbf8c599692102c3d387a37ef02a:08144a66f00dc7684fad061f1466033c0176e7ad-0' :
        'biocontainers/mulled-v2-77320db00eefbbf8c599692102c3d387a37ef02a:08144a66f00dc7684fad061f1466033c0176e7ad-0' }"

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
    pep_len = params.peptide_length ?: '9'
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
