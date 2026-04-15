process MOKAPOT {
    tag "${meta.id}_${engine}"
    label 'process_medium'

    conda "bioconda::mokapot=0.10.1"

    input:
    tuple val(meta), path(pin_files)
    val(engine)

    output:
    tuple val(meta), path("*.target.psms.txt"),     emit: psms
    tuple val(meta), path("*.target.peptides.txt"), emit: peptides
    tuple val(meta), path("*.decoy.psms.txt"),      emit: decoy_psms,     optional: true
    tuple val(meta), path("*.decoy.peptides.txt"),  emit: decoy_peptides, optional: true
    tuple val(meta), path("*_combined.pin"),        emit: combined_pin
    path("mokapot_log.txt"),                        emit: log
    path "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${engine}_search"
    """
    # Concatenate all PIN files, keeping header from the first only
    awk 'NR == 1 || FNR > 1' ${pin_files} > ${prefix}_combined.pin

    # Engine-specific PIN cleanup
    if [ "${engine}" = "comet" ]; then
        # Comet: merge trailing tab-delimited alternative proteins
        python3 -c "
import sys
with open('${prefix}_combined.pin') as f, open('tmp.pin','w') as o:
    for i, line in enumerate(f):
        s = line.split('\\\\t')
        if i == 0:
            ncol = len(s) - 1
        o.write('\\\\t'.join(s[:ncol]) + '\\\\t' + ';'.join(s[ncol:]))
"
        mv tmp.pin ${prefix}_combined.pin
    elif [ "${engine}" = "msfragger" ]; then
        # MSFragger: remove trailing semicolons
        awk '{sub(/;\$/, ""); print}' ${prefix}_combined.pin > tmp.pin
        mv tmp.pin ${prefix}_combined.pin
    fi

    # Run mokapot with 5% test FDR (needed for low-ID samples)
    mokapot \\
        --keep_decoys \\
        --test_fdr 0.05 \\
        -d . \\
        -r ${prefix} \\
        ${prefix}_combined.pin \\
        > mokapot_log.txt 2>&1

    # Normalise output file names so target/decoy globs are unambiguous.
    [ -f ${prefix}.mokapot.psms.txt ]           && mv ${prefix}.mokapot.psms.txt           ${prefix}.target.psms.txt
    [ -f ${prefix}.mokapot.peptides.txt ]       && mv ${prefix}.mokapot.peptides.txt       ${prefix}.target.peptides.txt
    [ -f ${prefix}.mokapot.decoy.psms.txt ]     && mv ${prefix}.mokapot.decoy.psms.txt     ${prefix}.decoy.psms.txt
    [ -f ${prefix}.mokapot.decoy.peptides.txt ] && mv ${prefix}.mokapot.decoy.peptides.txt ${prefix}.decoy.peptides.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mokapot: \$(mokapot --version 2>&1 | sed 's/mokapot //')
    END_VERSIONS
    """

    stub:
    def prefix = "${engine}_search"
    """
    touch ${prefix}.target.psms.txt ${prefix}.target.peptides.txt
    touch ${prefix}.decoy.psms.txt  ${prefix}.decoy.peptides.txt
    touch ${prefix}_combined.pin
    touch mokapot_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mokapot: "stub"
    END_VERSIONS
    """
}
