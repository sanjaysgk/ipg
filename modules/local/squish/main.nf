process SQUISH {
    tag "${meta.id}"
    label 'process_high'

    container "ghcr.io/sanjaysgk/ipg-tools:0.1.0"

    input:
    tuple val(meta), path(fastas)

    output:
    tuple val(meta), path("${meta.id}_cryptic.fasta"), emit: fasta
    tuple val(meta), path("*.log"),                    emit: log
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // `fastas` may be a single Path (when called with one input,
    // e.g. params.include_variant_peptides == false → reference-only)
    // or a List<Path> (when called with all three branches). Normalise
    // so the -d arg construction is uniform.
    def fasta_list = fastas instanceof List ? fastas : [fastas]
    def db_args    = fasta_list.collect { f -> "-d ${f}" }.join(' ')
    """
    squish ${db_args} \\
        -t ${task.cpus} \\
        -o ${meta.id}_cryptic.fasta \\
        > ${meta.id}_squish.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        squish: "kescull/immunopeptidogenomics@fef8e68"
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_cryptic.fasta
    touch ${meta.id}_squish.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        squish: stub
    END_VERSIONS
    """
}
