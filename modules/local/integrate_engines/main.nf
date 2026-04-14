process INTEGRATE_ENGINES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-bd0ce31f08f73d42bb0edebd8f5c7dcd9b1eefb0:1a671ad7f1b9f7bae925b3fe7baf49e06e59e70e-0' :
        'biocontainers/mulled-v2-bd0ce31f08f73d42bb0edebd8f5c7dcd9b1eefb0:1a671ad7f1b9f7bae925b3fe7baf49e06e59e70e-0' }"

    input:
    tuple val(meta), path(engine_tsvs, stageAs: '?/*')
    val(engine_names)            // list aligned with engine_tsvs
    path(fasta)

    output:
    tuple val(meta), path("integrated_psms.tsv"),       emit: psms
    tuple val(meta), path("integrated_peptides.tsv"),   emit: peptides
    tuple val(meta), path("chimeric_PSMs.txt"),         emit: chimeric,  optional: true
    tuple val(meta), path("chimera_only_peptides.txt"), emit: chim_only, optional: true
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def len_arg = params.peptide_length ?: '9'
    // Pair staged files with engine names; order is preserved.
    def pairs = [engine_names, engine_tsvs].transpose().collect { name, f -> "${name}=${f}" }.join(' ')
    """
    integrate_engines.py \\
        --engine-tsv ${pairs} \\
        --fasta ${fasta} \\
        --peptide-length ${len_arg} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch integrated_psms.tsv integrated_peptides.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
