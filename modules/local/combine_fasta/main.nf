process COMBINE_FASTA {
    tag "${prefix}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    path(fastas, stageAs: 'input/*')   // 2+ FASTAs concatenated in order; .gz or plain

    output:
    path("${prefix}.fasta"), emit: fasta
    path "versions.yml",     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: 'search_db'
    // Plain concatenation (NOT a dedup-merge): a peptide present in both the
    // canonical and cryptic FASTA must keep both records so INTEGRATE_ENGINES
    // can class it canonical by sp|/tr| prefix. zcat -f transparently passes
    // through uncompressed inputs, so gzipped (URL-fetched) and plain mix.
    """
    zcat -f ${fastas} > ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zcat: \$(zcat --version 2>/dev/null | head -1 | grep -oE '[0-9]+\\.[0-9]+' | head -1 || echo NA)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: 'search_db'
    """
    zcat -f ${fastas} > ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        zcat: "stub"
    END_VERSIONS
    """
}
