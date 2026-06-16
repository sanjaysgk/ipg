process REVERT_HEADERS {
    tag "${meta.id}"
    label 'process_single'

    container "ghcr.io/sanjaysgk/ipg-tools@sha256:e6999bae43e9a2b1b85497bbf02c3e0ba9f40c53a266d00a4270e2aa4fc7e7d5"

    input:
    tuple val(meta), path(alt_fasta)
    path reference_fasta

    output:
    tuple val(meta), path("${meta.id}_${alt_fasta.baseName}_reverted.fasta"), emit: fasta
    path "versions.yml",                                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}_${alt_fasta.baseName}_reverted"
    """
    # The improved revert_headers (sanjaysgk/immunopeptidogenomics@a09a74c)
    # accepts an optional 3rd positional arg for the output file prefix.
    # Output is written directly to <prefix>.fasta — no mv tmpc.fasta dance.
    revert_headers ${reference_fasta} ${alt_fasta} ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        revert_headers: "sanjaysgk/immunopeptidogenomics@a09a74c"
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
