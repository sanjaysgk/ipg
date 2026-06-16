process ALT_LIFTOVER {
    tag "${meta.id}"
    label 'process_single'

    container "ghcr.io/sanjaysgk/ipg-tools@sha256:5301688e40a8ea3e0ebb108f32e1b8ece95353eef444a54d9597a19ed6a0a8c8"

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
