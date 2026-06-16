process TRIPLE_TRANSLATE {
    tag "${meta.id}"
    label 'process_medium'

    container "ghcr.io/sanjaysgk/ipg-tools:sha-97c9c7e"

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
    # triple_translate writes <input_basename>_3translate.fasta to the
    # directory of the input file. Nextflow stages inputs as symlinks
    # pointing at upstream work dirs, so we must materialise the input
    # locally to make the output land here (where the *_3translate.fasta
    # glob can find it).
    cp -L ${transcriptome_fasta} ${prefix}.fasta

    # Normalise the gffcompare .tracking column 1: newer gffcompare writes it as
    # TCONS_xxxxxxxx|<isoforms>|<length>, but triple_translate matches the bare
    # TCONS id (from the transcriptome header) against col1 with an exact strcmp,
    # so the suffix silently disables all ENSG|ENST gene/class annotation. Strip
    # it (idempotent: already-bare col1 is left unchanged).
    sed 's/^\\(TCONS_[0-9][0-9]*\\)|[^[:space:]]*/\\1/' ${tracking} > tracking_norm.tracking

    triple_translate -c tracking_norm.tracking ${prefix}.fasta > ${prefix}_tt.log 2>&1
    rm ${prefix}.fasta tracking_norm.tracking

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
