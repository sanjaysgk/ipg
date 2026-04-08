process ALT_LIFTOVER {
    tag "${meta.id}"
    label 'process_single'

    container "ghcr.io/sanjaysgk/ipg-tools:0.1.0"

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
