process MSFRAGGER {
    tag "${meta.id}"
    label 'process_high'

    // Two paths supported:
    //   1. Bioconda: `msfragger` wrapper (PATH) — default under -profile pixi
    //      when pixi.toml includes `msfragger`. Needs a one-time license
    //      activation: `pixi run msfragger --key <license-key>`.
    //   2. User-supplied JAR via --msfragger_jar — for containerised runs
    //      or when the bioconda package isn't wanted.
    conda "bioconda::msfragger=4.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/msfragger:4.1--hdfd78af_0' :
        'biocontainers/msfragger:4.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(ms_files)
    path(fasta)
    path(msfragger_jar)    // optional — pass file('NO_FILE') when using bioconda
    path(params_file)

    output:
    tuple val(meta), path("*.pin"),            emit: pin
    // mzml is optional — MSFragger only produces a *_calibrated.mzML when
    // mass calibration runs. On small / DIA / already-calibrated inputs
    // it skips calibration and emits no new mzML. Downstream modules
    // (CONVERT_MZML, COMET, SAGE) fall back to the original mzML via
    // `ch_engine_input = engines.contains('msfragger') ? ch_calibrated_mzml : ch_ms_data`
    // in the subworkflow — left as-is, and COMET/SAGE receive ch_ms_data.
    tuple val(meta), path("*.mzML"),           emit: mzml, optional: true
    path("search_log.txt"),                    emit: log
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mem = task.ext.msfragger_mem ?: params.msfragger_mem ?: 8
    def use_jar = msfragger_jar && msfragger_jar.size() > 0
    def key_arg = (!use_jar && params.msfragger_license) ? "--key ${params.msfragger_license}" : ''
    def run_cmd = use_jar
        ? "java -Dfile.encoding=UTF-8 -Xmx${mem}g -jar ${msfragger_jar}"
        : "msfragger ${key_arg}"
    def version_cmd = use_jar
        ? "java -jar ${msfragger_jar} --version 2>&1 | head -1 || echo unknown"
        : "echo 'MSFragger 4.1 (bioconda)'"
    """
    # bioconda msfragger wrapper hardcodes -Xmx1g which is too small;
    # _JAVA_OPTIONS overrides it without touching the wrapper script.
    export _JAVA_OPTIONS="-Xmx${mem}g -Dfile.encoding=UTF-8"

    sed 's#your_fasta#${fasta}#g' ${params_file} > used.params

    ${run_cmd} \\
        used.params \\
        ${ms_files} \\
        > search_log.txt 2>&1

    # Rename calibrated mzMLs to clean names (remove _calibrated suffix).
    # Use nullglob so unmatched globs produce empty list instead of the
    # literal pattern (bash strict mode -e would kill the script otherwise).
    shopt -s nullglob
    for f in *_calibrated.mzML; do
        mv "\$f" "\${f/_calibrated/}"
    done
    # If calibration failed, use uncalibrated
    for f in *_uncalibrated.mzML; do
        target="\${f/_uncalibrated/}"
        if [ ! -f "\$target" ]; then
            mv "\$f" "\$target"
        else
            rm -f "\$f"
        fi
    done
    shopt -u nullglob


    # Clean up pepindex files
    rm -f *.pepindex

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    msfragger: \$(${version_cmd})
END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.pin
    touch ${meta.id}.mzML
    touch search_log.txt

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    msfragger: "stub"
END_VERSIONS
    """
}
