process MS2RESCORE {
    tag "${meta.id}_${engine}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ms2rescore:3.1.5--pyhdfd78af_0' :
        'biocontainers/ms2rescore:3.1.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(mokapot_target), path(mokapot_decoy), path(pin_file)
    val(engine)
    path(scans_pkls, stageAs: 'scans/*')
    path(mgf_files,  stageAs: 'scans/*')
    path(config_json)

    output:
    tuple val(meta), path("${engine}.psms.tsv"),   emit: psms
    tuple val(meta), path("${engine}_rescore_input.tsv"), emit: input_tsv
    tuple val(meta), path("${engine}.report.html"), emit: report, optional: true
    path("${engine}_pipeline_log.txt"),            emit: log
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mod = params.mod_type ?: 'mod'
    def rescore_engine = params.ms2rescore_engine ?: 'mokapot'
    """
    prepare_ms2rescore_input.py \\
        --engine ${engine} \\
        --target ${mokapot_target} \\
        --decoy  ${mokapot_decoy} \\
        --pin    ${pin_file} \\
        --scans  ${scans_pkls.join(' ')} \\
        --mod    ${mod} \\
        --out    ${engine}_rescore_input.tsv

    # Select the final rescoring engine. Default keeps the config's mokapot
    # block; percolator swaps it in (percolator must be on PATH).
    cp ${config_json} run_config.json
    if [ "${rescore_engine}" = "percolator" ]; then
        python3 -c "import json; c=json.load(open('run_config.json')); c['ms2rescore']['rescoring_engine']={'percolator':{}}; json.dump(c, open('run_config.json','w'), indent=4)"
    fi

    ms2rescore \\
        -t tsv \\
        -s scans \\
        -p ${engine}_rescore_input.tsv \\
        -c run_config.json \\
        -o ${engine} \\
        > ${engine}_pipeline_log.txt 2>&1

    # ms2rescore outputs ${engine}.psms.tsv + other files in cwd; normalise name.
    if [ ! -f "${engine}.psms.tsv" ] && [ -f "${engine}.ms2rescore.psms.tsv" ]; then
        mv ${engine}.ms2rescore.psms.tsv ${engine}.psms.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ms2rescore: \$(ms2rescore --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${engine}.psms.tsv ${engine}_rescore_input.tsv ${engine}.report.html ${engine}_pipeline_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ms2rescore: "stub"
    END_VERSIONS
    """
}
