process MSFRAGGER_PARTIAL {
    tag "${meta.id}"
    label 'process_high'

    // Standalone bioconda MSFragger (same calling convention as modules/local/msfragger):
    //   bioconda msfragger with the key from the MSFRAGGER_LICENSE secret, or a user
    //   JAR via --msfragger_jar.
    // NOT FragPipe — split-database search is an MSFragger feature (--partial), so we
    // call MSFragger directly. One invocation searches a sample's spectra against ONE
    // database chunk in partial mode.
    secret 'MSFRAGGER_LICENSE'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/msfragger:4.2--py311hdfd78af_0' :
        'biocontainers/msfragger:4.2--py311hdfd78af_0' }"

    input:
    tuple val(meta), path(ms_files)   // meta.chunk_id = this chunk's index; meta.original_id = sample
    path(fasta)                       // ONE database chunk
    path(msfragger_jar)               // optional — file('NO_FILE') for bioconda
    path(params_file)

    output:
    // Emit ONE directory per chunk, named chunk_<id>, containing MSFragger's
    // outputs named by spectrum-file basename (<run>.pepXML, <run>.pin,
    // <run>_scores_histogram.tsv). merge_split_search.py takes these as
    // --chunk_dirs and finds files by sample name (= the run basename).
    tuple val(meta), path("chunk_${meta.chunk_id}"), emit: chunk_dir
    path("search_log.txt"),                          emit: log
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mem     = task.ext.msfragger_mem ?: params.msfragger_mem ?: 8
    def use_jar = msfragger_jar && msfragger_jar.size() > 0
    // licence key from the MSFRAGGER_LICENSE secret env var (single-quoted so Groovy
    // leaves $MSFRAGGER_LICENSE literal; shell resolves it at runtime, not in .command.sh).
    def key_arg = use_jar ? '' : '--key $MSFRAGGER_LICENSE'
    def run_cmd = use_jar
        ? "java -Dfile.encoding=UTF-8 -Xmx${mem}g -jar ${msfragger_jar}"
        : "msfragger ${key_arg}"
    def version_cmd = use_jar
        ? "java -jar ${msfragger_jar} --version 2>&1 | grep -m1 -oE 'MSFragger-[0-9.]+' || echo unknown"
        : "echo 'MSFragger 4.2 (bioconda)'"
    def chunk = meta.chunk_id
    """
    # bioconda wrapper hardcodes -Xmx1g; override via _JAVA_OPTIONS (don't touch the wrapper).
    export _JAVA_OPTIONS="-Xmx${mem}g -Dfile.encoding=UTF-8"

    sed 's#your_fasta#${fasta}#g' ${params_file} > used.params
    # CRITICAL for statistical validity: disable per-chunk mass calibration so every
    # chunk's score histogram is on the SAME scale. The merge sums histograms across
    # chunks into the full-DB score distribution and recomputes e-values from it —
    # that only holds if the per-chunk distributions are comparable. (Calibrate once
    # up front instead, if/when a calibrate step is added; never per-chunk.)
    sed -i '/^calibrate_mass/d' used.params
    echo "calibrate_mass = 0" >> used.params
    sed -i '/^write_calibrated_mzml/d' used.params
    echo "write_calibrated_mzml = 0" >> used.params

    # --partial <chunk> emits a PARTIAL pepXML (raw hyperscores, NOT per-chunk e-values)
    # + a per-spectrum score histogram, designed for the full-DB-recompute merge.
    ${run_cmd} \\
        used.params \\
        ${ms_files} \\
        --partial ${chunk} \\
        > search_log.txt 2>&1

    # Collect this chunk's outputs into chunk_<id>/ — keeping MSFragger's
    # run-basename filenames, which is what merge_split_search.py expects to find
    # inside each --chunk_dirs entry (sample_name = run basename).
    mkdir -p chunk_${chunk}
    shopt -s nullglob
    for f in *.pepXML *.pin *_scores_histogram.tsv; do
        mv "\$f" "chunk_${chunk}/"
    done
    shopt -u nullglob

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msfragger: \$(${version_cmd})
    END_VERSIONS
    """

    stub:
    def chunk = meta.chunk_id
    """
    mkdir -p chunk_${chunk}
    touch chunk_${chunk}/${meta.id}.pepXML
    touch chunk_${chunk}/${meta.id}.pin
    touch chunk_${chunk}/${meta.id}_scores_histogram.tsv
    touch search_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msfragger: "stub"
    END_VERSIONS
    """
}
