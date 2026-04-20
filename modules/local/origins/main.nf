process ORIGINS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/sanjaysgk/ipg-tools:0.2.0"

    input:
    tuple val(meta), path(peptide_list)
    path(tracking)
    path(uniprot_fasta)
    path(transcriptome_fasta)

    output:
    tuple val(meta), path("*_origins_discard.txt"),        emit: discard,        optional: true
    tuple val(meta), path("*_origins_unconventional.txt"), emit: unconventional, optional: true
    tuple val(meta), path("*_origins_*.csv"),              emit: origins_csv,    optional: true
    tuple val(meta), path("*.txt"),                        emit: all_txt
    path "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    origins ${args} \\
        -p ${peptide_list} \\
        -r ${tracking} \\
        -d ${uniprot_fasta} \\
        -n ${transcriptome_fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        origins: "kescull/immunopeptidogenomics"
    END_VERSIONS
    """

    stub:
    def prefix = peptide_list.baseName
    """
    touch ${prefix}_origins_discard.txt
    touch ${prefix}_origins_unconventional.txt
    touch ${prefix}_origins_prot.csv ${prefix}_origins_rna.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        origins: "kescull/immunopeptidogenomics"
    END_VERSIONS
    """
}
