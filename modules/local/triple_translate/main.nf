process TRIPLE_TRANSLATE {
    tag "${meta.id}"
    label 'process_medium'

    container "ghcr.io/sanjaysgk/ipg-tools:0.1.0"

    input:
    tuple val(meta), path(transcriptome_fasta), path(tracking)

    output:
    tuple val(meta), path("*_3translate.fasta"), emit: fasta
    tuple val(meta), path("*.log"),              emit: log
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${transcriptome_fasta.baseName}"
    """
    triple_translate -c ${tracking} ${transcriptome_fasta} > ${prefix}_tt.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        triple_translate: "kescull/immunopeptidogenomics@fef8e68"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${transcriptome_fasta.baseName}"
    """
    touch ${prefix}_3translate.fasta
    touch ${prefix}_tt.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        triple_translate: stub
    END_VERSIONS
    """
}
