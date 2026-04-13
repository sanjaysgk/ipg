process DB_COMPARE {
    tag "${meta.id}"
    label 'process_single'

    // R script lives in pipeline bin/ (auto-added to PATH by Nextflow).
    // No container needed for pixi profile; for container profiles,
    // use the rocker/tidyverse image with additional CRAN packages.
    conda "conda-forge::r-base=4.3 conda-forge::r-tidyverse conda-forge::r-optparse conda-forge::r-venndiagram conda-forge::r-upsetr conda-forge::r-hrbrthemes"

    input:
    tuple val(meta), path(cryptic_psm), path(uniprot_psm)
    path(discard_list)
    path(unconventional_list)

    output:
    tuple val(meta), path("*_cryptic_only.txt"),                   emit: cryptic_only,               optional: true
    tuple val(meta), path("*_unambiguous_unconventional.txt"),     emit: unambiguous_unconventional,  optional: true
    tuple val(meta), path("*.png"),                                emit: plots,                       optional: true
    tuple val(meta), path("*.csv"),                                emit: csvs,                        optional: true
    tuple val(meta), path("*.txt"),                                emit: all_txt
    path "versions.yml",                                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def cryptic_score    = meta.cryptic_decoy_score
    def uniprot_score    = meta.uniprot_decoy_score
    def phase2_args      = discard_list.name != 'NO_FILE' ? "-j ${discard_list} -u ${unconventional_list}" : ''
    """
    db_compare_v2.R \\
        -c ${cryptic_psm} \\
        -n ${uniprot_psm} \\
        -d ${cryptic_score} \\
        -m ${uniprot_score} \\
        -p ${prefix} \\
        ${phase2_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
        db_compare: "kescull/immunopeptidogenomics"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_cryptic_only.txt
    touch ${prefix}_unambiguous_unconventional.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
        db_compare: "kescull/immunopeptidogenomics"
    END_VERSIONS
    """
}
