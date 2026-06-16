process MSFRAGGER {
    tag "${meta.id}"
    label 'process_high'

    // Two paths supported:
    //   1. Bioconda: `msfragger` wrapper (PATH) — default under -profile pixi
    //      when pixi.toml includes `msfragger`. The academic licence key is read
    //      from the MSFRAGGER_LICENSE Nextflow secret (injected as an env var at
    //      runtime, never written to .command.sh). MSFragger does not persist
    //      activation, so the key must be supplied on every run. Set it once with:
    //          nextflow secrets set MSFRAGGER_LICENSE '<key>'
    //   2. User-supplied JAR via --msfragger_jar — for containerised runs
    //      or when the bioconda package isn't wanted.
    secret 'MSFRAGGER_LICENSE'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/msfragger:4.2--py311hdfd78af_0' :
        'biocontainers/msfragger:4.2--py311hdfd78af_0' }"

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
    tuple val(meta), path("search_log.txt"),   emit: log
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mem = task.ext.msfragger_mem ?: params.msfragger_mem ?: 8
    def use_jar = msfragger_jar && msfragger_jar.size() > 0
    // licence key from the MSFRAGGER_LICENSE secret as an env var — single-quoted
    // so Groovy leaves $MSFRAGGER_LICENSE literal; the shell resolves it at runtime,
    // keeping the key out of .command.sh.
    def key_arg = use_jar ? '' : '--key $MSFRAGGER_LICENSE'
    def run_cmd = use_jar
        ? "java -Dfile.encoding=UTF-8 -Xmx${mem}g -jar ${msfragger_jar}"
        : "msfragger ${key_arg}"
    // grep the version token only: with _JAVA_OPTIONS exported the JVM prints a
    // "Picked up _JAVA_OPTIONS: ..." banner whose colons would corrupt versions.yml.
    def version_cmd = use_jar
        ? "java -jar ${msfragger_jar} --version 2>&1 | grep -m1 -oE 'MSFragger-[0-9.]+' || echo unknown"
        : "echo 'MSFragger 4.2 (bioconda)'"
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
