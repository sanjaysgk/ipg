process CRYPTIC_REPORT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-77320db00eefbbf8c599692102c3d387a37ef02a:08144a66f00dc7684fad061f1466033c0176e7ad-0' :
        'biocontainers/mulled-v2-77320db00eefbbf8c599692102c3d387a37ef02a:08144a66f00dc7684fad061f1466033c0176e7ad-0' }"

    input:
    tuple val(meta), path(peptides_tsv), path(rescored_tsvs, stageAs: 'rescored/??/*')

    output:
    tuple val(meta), path("${meta.id}_cryptic_report.html"), emit: html
    tuple val(meta), path("${meta.id}_cryptic_summary.tsv"), emit: summary
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    sample = meta.id
    canonical_prefixes = params.canonical_protein_prefixes ?: 'sp|,tr|'
    template 'cryptic_report.py'

    stub:
    """
    touch ${meta.id}_cryptic_report.html ${meta.id}_cryptic_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
