process SPLIT_FASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(fasta)
    val num_chunks

    output:
    tuple val(meta), path("split_db/*/*.fasta"), emit: fasta_chunks
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    // Split the FASTA into num_chunks parts at protein boundaries for MSFragger
    // split-database search (each chunk is searched with --partial, then merged
    // with full-DB e-value recompute — statistically identical to an unsplit search).
    // bin/split_fasta.py is a standalone port of FragPipe's get_fasta_offsets().
    script:
    """
    split_fasta.py --fasta ${fasta} --num_chunks ${num_chunks} --outdir split_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    for i in \$(seq 0 \$((${num_chunks} - 1))); do
        mkdir -p split_db/\$i
        printf '>chunk%s_protein\\nPEPTIDEK\\n' "\$i" > split_db/\$i/${fasta.baseName}.fasta
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
