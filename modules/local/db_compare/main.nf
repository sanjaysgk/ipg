process DB_COMPARE {
    tag "${meta.id}"
    label 'process_single'

    // R script lives in pipeline bin/ (auto-added to PATH by Nextflow).
    // No container needed for pixi profile; for container profiles,
    // use the rocker/tidyverse image with additional CRAN packages.
    conda "conda-forge::r-base=4.3 conda-forge::r-tidyverse conda-forge::r-optparse conda-forge::r-venndiagram conda-forge::r-upsetr conda-forge::r-hrbrthemes"

    input:
    tuple val(meta), path(cryptic_psm), path(uniprot_psm)
    path(discard_list,       stageAs: 'discard/*')
    path(unconventional_list, stageAs: 'unconv/*')

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
    // Use file size (NO_FILE is the zero-byte sentinel) rather than name —
    // `stageAs: 'discard/*'` means `.name` returns the staged subpath, not
    // the basename, which broke the previous string-based check.
    def phase2_args      = discard_list.size() > 0 ? "-j ${discard_list} -u ${unconventional_list}" : ''
    """
    # Normalise the score column header so db_compare_v2.R's hard-coded
    # `X.10lgP` (lowercase l) selector matches. Newer PEAKS exports use
    # `-10LgP` (capital L), which R's check.names mangles to `X.10LgP`
    # and breaks the dplyr select. Operate on copies so we don't touch
    # the staged input files in place.
    cp ${cryptic_psm} cryptic_in.csv
    cp ${uniprot_psm} uniprot_in.csv
    sed -i '1s/-10LgP/-10lgP/g' cryptic_in.csv uniprot_in.csv

    db_compare_v2.R \\
        -c cryptic_in.csv \\
        -n uniprot_in.csv \\
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
