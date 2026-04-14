process IMMUNOINFORMATICS_REPORT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-bd0ce31f08f73d42bb0edebd8f5c7dcd9b1eefb0:1a671ad7f1b9f7bae925b3fe7baf49e06e59e70e-0' :
        'biocontainers/mulled-v2-bd0ce31f08f73d42bb0edebd8f5c7dcd9b1eefb0:1a671ad7f1b9f7bae925b3fe7baf49e06e59e70e-0' }"

    input:
    tuple val(meta),
        path(peptides),
        path(netmhcpan,   stageAs: 'netmhcpan.tsv'),
        path(netmhciipan, stageAs: 'netmhciipan.tsv'),
        path(gibbs,       stageAs: 'gibbs.tsv'),
        path(flashlfq,    stageAs: 'flashlfq.tsv'),
        path(blastp,      stageAs: 'blastp.tsv')

    output:
    tuple val(meta), path("${meta.id}_immunoinformatics_report.html"), emit: html
    path "versions.yml",                                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def netmhcpan_arg   = netmhcpan.name   == 'NO_FILE' ? ''  : "--netmhcpan ${netmhcpan}"
    def netmhciipan_arg = netmhciipan.name == 'NO_FILE' ? ''  : "--netmhciipan ${netmhciipan}"
    def gibbs_arg       = gibbs.name       == 'NO_FILE' ? ''  : "--gibbs ${gibbs}"
    def flashlfq_arg    = flashlfq.name    == 'NO_FILE' ? ''  : "--flashlfq ${flashlfq}"
    def blastp_arg      = blastp.name      == 'NO_FILE' ? ''  : "--blastp ${blastp}"
    """
    generate_report.py \\
        --sample ${meta.id} \\
        --peptides ${peptides} \\
        ${netmhcpan_arg} ${netmhciipan_arg} ${gibbs_arg} ${flashlfq_arg} ${blastp_arg} \\
        --out ${meta.id}_immunoinformatics_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
        logomaker: \$(python3 -c "import logomaker; print(logomaker.__version__)" 2>/dev/null || echo "not installed")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_immunoinformatics_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
