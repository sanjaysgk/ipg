process COMET {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::comet-ms=2024010"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/comet-ms:2024010--h43eeafb_0' :
        'biocontainers/comet-ms:2024010--h43eeafb_0' }"

    input:
    tuple val(meta), path(mzml_files)
    path(fasta)
    path(params_file)

    output:
    tuple val(meta), path("*.pin"), emit: pin
    path("comet_search_log.txt"),   emit: log
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Comet reads FASTA path and params; pass params file directly.
    # The FASTA path inside params is ignored when -D is given on the CLI.
    comet \\
        -D${fasta} \\
        -P${params_file} \\
        ${mzml_files} \\
        > comet_search_log.txt 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        comet: \$(comet 2>&1 | grep -oE 'Comet version "[^"]+"' | head -1 | sed 's/Comet version //;s/"//g')
    END_VERSIONS
    """

    stub:
    """
    for f in ${mzml_files}; do
        base=\$(basename "\$f"); base="\${base%.*}"
        touch "\${base}.pin"
    done
    touch comet_search_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        comet: "stub"
    END_VERSIONS
    """
}
