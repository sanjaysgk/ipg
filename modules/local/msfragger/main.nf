process MSFRAGGER {
    tag "${meta.id}"
    label 'process_high'

    // MSFragger requires Java 11+ and the user-provided JAR.
    // No conda env needed beyond openjdk — the JAR is self-contained.
    conda "conda-forge::openjdk=17"

    input:
    tuple val(meta), path(ms_files)
    path(fasta)
    path(msfragger_jar)
    path(params_file)

    output:
    tuple val(meta), path("*.pin"),            emit: pin
    tuple val(meta), path("*.mzML"),           emit: mzml
    path("search_log.txt"),                    emit: log
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mem = task.ext.msfragger_mem ?: params.msfragger_mem ?: 8
    def prefix = meta.id
    // Substitute the FASTA path into the params file
    """
    sed 's#your_fasta#${fasta}#g' ${params_file} > used.params

    java -Dfile.encoding=UTF-8 -Xmx${mem}g \\
        -jar ${msfragger_jar} \\
        used.params \\
        ${ms_files} \\
        > search_log.txt 2>&1

    # Rename calibrated mzMLs to clean names (remove _calibrated suffix)
    for f in *_calibrated.mzML; do
        [ -f "\$f" ] && mv "\$f" "\${f/_calibrated/}" || true
    done
    # If calibration failed, use uncalibrated
    for f in *_uncalibrated.mzML; do
        target="\${f/_uncalibrated/}"
        [ ! -f "\$target" ] && mv "\$f" "\$target" || rm -f "\$f"
    done

    # Clean up pepindex files
    rm -f *.pepindex

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -1 | sed 's/.*"\\(.*\\)".*/\\1/')
        msfragger: \$(java -jar ${msfragger_jar} --version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.pin
    touch ${meta.id}.mzML
    touch search_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -1 | sed 's/.*"\\(.*\\)".*/\\1/')
        msfragger: "stub"
    END_VERSIONS
    """
}
