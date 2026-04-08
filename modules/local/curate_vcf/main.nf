process CURATE_VCF {
    tag "${meta.id}"
    label 'process_single'

    container "ghcr.io/sanjaysgk/ipg-tools:0.1.0"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_curated_unmasked.vcf"), emit: unmasked
    tuple val(meta), path("*_curated_indel.vcf"),    emit: indel
    tuple val(meta), path("*.log"),                  emit: log
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when.every { it }

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Unmasked curation (keeps all substitutions + deletions)
    curate_vcf ${vcf} > ${prefix}_curate_unmasked.log 2>&1
    mv \$(basename ${vcf} .vcf)_curated.vcf ${prefix}_curated_unmasked.vcf

    # Indel-only curation (deprioritises deletions that would mask downstream sites)
    curate_vcf -d ${vcf} > ${prefix}_curate_indel.log 2>&1
    mv \$(basename ${vcf} .vcf)_curated.vcf ${prefix}_curated_indel.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curate_vcf: \$(curate_vcf -h 2>&1 | head -1 || echo "kescull/immunopeptidogenomics@fef8e68")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_curated_unmasked.vcf
    touch ${prefix}_curated_indel.vcf
    touch ${prefix}_curate_unmasked.log
    touch ${prefix}_curate_indel.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curate_vcf: stub
    END_VERSIONS
    """
}
