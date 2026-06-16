process ALT_LIFTOVER {
    tag "${meta.id}"
    label 'process_single'

    container "ghcr.io/sanjaysgk/ipg-tools@sha256:e6999bae43e9a2b1b85497bbf02c3e0ba9f40c53a266d00a4270e2aa4fc7e7d5"

    input:
    tuple val(meta), path(vcf), path(gtf)
    val suffix

    output:
    tuple val(meta), path("*${suffix}.gtf"), emit: gtf
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    alt_liftover -s ${suffix} -v ${vcf} -g ${gtf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alt_liftover: "kescull/immunopeptidogenomics@fef8e68"
    END_VERSIONS
    """

    stub:
    def gtf_base = gtf.baseName
    """
    touch ${gtf_base}${suffix}.gtf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alt_liftover: stub
    END_VERSIONS
    """
}
