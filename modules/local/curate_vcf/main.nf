process CURATE_VCF {
    tag "${meta.id}"
    label 'process_single'

    container "ghcr.io/sanjaysgk/ipg-tools:0.1.0"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_unmasked.vcf"), emit: unmasked
    tuple val(meta), path("*_indel.vcf"),    emit: indel
    tuple val(meta), path("*.log"),          emit: log
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // curate_vcf writes <stem>_indel.vcf when called without -d and
    // <stem>_unmasked.vcf when called with -d, where <stem> is the input
    // basename with the .vcf extension stripped. See curate_vcf.c:77-79.
    //
    // The C tool reads the VCF as plain text so compressed inputs must be
    // decompressed first. We handle both .vcf and .vcf.gz transparently —
    // the pipeline passes in whatever gatk4/selectvariants emitted (.vcf.gz
    // by default under conf/modules.config).
    """
    set -eo pipefail
    case "${vcf}" in
        *.vcf.gz)
            stem=\$(basename ${vcf} .vcf.gz)
            gunzip -c ${vcf} > "\${stem}.vcf"
            input_vcf="\${stem}.vcf"
            ;;
        *.vcf)
            stem=\$(basename ${vcf} .vcf)
            input_vcf="${vcf}"
            ;;
        *)
            echo "curate_vcf: unsupported input extension for ${vcf}" >&2
            exit 2
            ;;
    esac

    # The kescull curate_vcf binary segfaults on empty VCFs (no variant
    # records). Guard against that by short-circuiting to empty output
    # when there are zero non-header lines — this happens routinely on
    # tiny test data and on real samples that simply have no somatic
    # variants in the analysed region.
    n_variants=\$(grep -vc '^#' "\${input_vcf}" || true)
    echo "curate_vcf: \${n_variants} variant record(s) in input" \\
        | tee ${prefix}_curate_indel.log ${prefix}_curate_unmasked.log

    if [ "\${n_variants}" -eq 0 ]; then
        # Pass through the (empty) input to both downstream branches.
        cp "\${input_vcf}" ${prefix}_indel.vcf
        cp "\${input_vcf}" ${prefix}_unmasked.vcf
    else
        # Indel-preserving curation (no -d): deletions are KEPT, so
        # downstream indel sites that overlap deletions are still callable.
        curate_vcf "\${input_vcf}" >> ${prefix}_curate_indel.log 2>&1
        mv "\${stem}_indel.vcf" ${prefix}_indel.vcf

        # Deletion-deprioritised curation (-d flag): deletions that would
        # mask downstream mutation sites are removed, producing an
        # 'unmasked' VCF.
        curate_vcf -d "\${input_vcf}" >> ${prefix}_curate_unmasked.log 2>&1
        mv "\${stem}_unmasked.vcf" ${prefix}_unmasked.vcf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curate_vcf: "kescull/immunopeptidogenomics@fef8e68"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_unmasked.vcf
    touch ${prefix}_indel.vcf
    touch ${prefix}_curate_unmasked.log
    touch ${prefix}_curate_indel.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curate_vcf: stub
    END_VERSIONS
    """
}
