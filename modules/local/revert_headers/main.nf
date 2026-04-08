process REVERT_HEADERS {
    tag "${meta.id}"
    label 'process_single'

    container "ghcr.io/sanjaysgk/ipg-tools:0.1.0"

    input:
    tuple val(meta), path(alt_fasta)
    path reference_fasta

    output:
    tuple val(meta), path("${meta.id}_${alt_fasta.baseName}_reverted.fasta"), emit: fasta
    path "versions.yml",                                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # revert_headers writes to tmpc.fasta in cwd
    revert_headers ${reference_fasta} ${alt_fasta}
    mv tmpc.fasta ${meta.id}_${alt_fasta.baseName}_reverted.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        revert_headers: "kescull/immunopeptidogenomics@fef8e68"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${alt_fasta.baseName}_reverted.fasta
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        revert_headers: stub
    END_VERSIONS
    """
}
