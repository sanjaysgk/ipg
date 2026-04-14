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
    path("${engine}_pipeline_log.txt"),            emit: log
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mod = params.mod_type ?: 'mod'
    """
    prepare_ms2rescore_input.py \\
        --engine ${engine} \\
        --target ${mokapot_target} \\
        --decoy  ${mokapot_decoy} \\
        --pin    ${pin_file} \\
        --scans  ${scans_pkls.join(' ')} \\
        --mod    ${mod} \\
        --out    ${engine}_rescore_input.tsv

    ms2rescore \\
        -t tsv \\
        -s scans \\
        -p ${engine}_rescore_input.tsv \\
        -c ${config_json} \\
        -o ${engine} \\
        > ${engine}_pipeline_log.txt 2>&1

    # ms2rescore outputs ${engine}.psms.tsv + other files in cwd; normalise name.
    if [ ! -f "${engine}.psms.tsv" ] && [ -f "${engine}.ms2rescore.psms.tsv" ]; then
        mv ${engine}.ms2rescore.psms.tsv ${engine}.psms.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ms2rescore: \$(ms2rescore --version 2>&1 | awk '{print \$NF}' | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${engine}.psms.tsv ${engine}_rescore_input.tsv ${engine}_pipeline_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ms2rescore: "stub"
    END_VERSIONS
    """
}
