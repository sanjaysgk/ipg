process DENOVO_CLASSIFY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(predictions), path(cryptic_fasta)
    path canonical_fasta

    output:
    tuple val(meta), path("${prefix}.denovo_classified.tsv"), emit: classified
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    // Look up each de novo peptide (I/L-normalised substring match) in the cryptic
    // and (optional) canonical FASTA to assign class canonical/cryptic/novel.
    // Confidence filtering + column names come via ext.args.
    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def canon_arg = canonical_fasta ? "--canonical-fasta ${canonical_fasta}" : ''
    """
    denovo_classify.py \\
        --predictions ${predictions} \\
        --cryptic-fasta ${cryptic_fasta} \\
        ${canon_arg} \\
        --out ${prefix}.denovo_classified.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'spectrum_id\\tpeptide\\tlength\\tlog_prob\\tclass\\tengine\\nstub:1\\tPEPTIDEK\\t8\\t-0.1\\tcryptic\\tinstanovo\\n' > ${prefix}.denovo_classified.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
