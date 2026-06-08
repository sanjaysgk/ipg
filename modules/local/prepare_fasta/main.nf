process PREPARE_FASTA {
    tag "prepare_fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    path(fasta)

    output:
    path("*_tda.fasta"), emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = fasta.baseName
    template 'prepare_fasta.py'

    stub:
    def prefix = fasta.baseName
    """
    touch ${prefix}_tda.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
        prepare_fasta: "sanjaysgk/ipg"
    END_VERSIONS
    """
}
