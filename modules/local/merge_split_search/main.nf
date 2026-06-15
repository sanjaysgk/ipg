process MERGE_SPLIT_SEARCH {
    tag "${meta.id}"
    label 'process_medium'

    // Needs BOTH msfragger (for --generate_expect_functions on the combined
    // histogram) AND numpy (histogram math) — so it runs in the msfragger
    // environment with numpy added. NOTE for singularity: confirm the image
    // ships numpy, or supply a mulled msfragger+numpy image via ext.container.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/msfragger:4.2--py311hdfd78af_0' :
        'biocontainers/msfragger:4.2--py311hdfd78af_0' }"

    input:
    tuple val(meta), path(chunk_dirs, stageAs: 'chunks/*')   // N chunk_<id> dirs for ONE run
    val num_chunks
    path(msfragger_jar)                                      // optional — file('NO_FILE') for bioconda
    path(params_file)

    output:
    tuple val(meta), path("merged/*.pin"),    emit: pin
    tuple val(meta), path("merged/*.pepXML"), emit: pepxml, optional: true
    path "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    // bin/merge_split_search.py (faithful port of FragPipe msfragger_pep_split.py):
    //   1. SUM the per-chunk score histograms -> full-DB score distribution
    //   2. RECOMPUTE the e-value function from it via `<msfragger_cmd> --generate_expect_functions`
    //   3. re-rank hits per spectrum and apply the full-DB e-value
    // => e-values mathematically identical to an unsplit full-DB search.
    // sample_name MUST equal the run basename so the script finds <run>.pepXML /
    // <run>.pin / <run>_scores_histogram.tsv inside each chunk dir.
    script:
    def mem         = task.ext.msfragger_mem ?: params.msfragger_mem ?: 8
    def use_jar     = msfragger_jar && msfragger_jar.size() > 0
    def key_arg     = (!use_jar && params.msfragger_license) ? "--key ${params.msfragger_license}" : ''
    def msfragger_cmd = use_jar ? "java -Xmx${mem}g -jar ${msfragger_jar}" : "msfragger ${key_arg}"
    def sample_name = meta.run ?: meta.id
    def topN        = task.ext.output_report_topN ?: 1
    def max_expect  = task.ext.output_max_expect ?: 50.0
    """
    export _JAVA_OPTIONS="-Xmx${mem}g -Dfile.encoding=UTF-8"

    dirs=\$(ls -d chunks/* | paste -sd,)

    merge_split_search.py \\
        --sample_name ${sample_name} \\
        --chunk_dirs "\$dirs" \\
        --num_chunks ${num_chunks} \\
        --msfragger_cmd "${msfragger_cmd}" \\
        --outdir merged \\
        --output_report_topN ${topN} \\
        --output_max_expect ${max_expect}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def sample_name = meta.run ?: meta.id
    """
    mkdir -p merged
    touch merged/${sample_name}.pin
    touch merged/${sample_name}.pepXML

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
