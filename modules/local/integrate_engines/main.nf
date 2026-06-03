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
    len_arg = params.peptide_length ?: '9'
    fdr_arg = params.integrate_fdr ?: '0.01'
    canonical_prefixes = params.canonical_protein_prefixes ?: 'sp|,tr|'
    // Pair staged files with engine names; order is preserved. Coerce to lists first:
    // with a single engine Nextflow passes a scalar Path/value (not a list), and
    // transpose() would then iterate the Path's components (pairing name with the
    // '?/*' stage dir '1' instead of the file). Wrapping guarantees n=1 works too.
    def names = engine_names instanceof List ? engine_names : [engine_names]
    def files = engine_tsvs  instanceof List ? engine_tsvs  : [engine_tsvs]
    pairs = [names, files].transpose().collect { name, f -> "${name}=${f}" }.join(' ')
    template 'integrate_engines.py'

    stub:
    """
    touch integrated_psms.tsv integrated_peptides.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
