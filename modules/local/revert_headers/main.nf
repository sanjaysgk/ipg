process REVERT_HEADERS {
    tag "${meta.id}"
    label 'process_single'

    container "ghcr.io/sanjaysgk/ipg-tools@sha256:165874f737838f65b62b3b9421a52d13fa59a6724f15fd3ac06bc4a646026755"

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
