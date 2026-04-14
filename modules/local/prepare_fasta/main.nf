process PREPARE_FASTA {
    tag "prepare_fasta"
    label 'process_single'

    conda "conda-forge::python=3.12"

    input:
    path(fasta)

    output:
    path("*_tda.fasta"), emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = fasta.baseName
    """
    prepare_fasta.py \\
        -i ${fasta} \\
        -o ${prefix}_tda.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        prepare_fasta: "sanjaysgk/ipg"
    END_VERSIONS
    """

    stub:
    def prefix = fasta.baseName
    """
    touch ${prefix}_tda.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        prepare_fasta: "sanjaysgk/ipg"
    END_VERSIONS
    """
}
