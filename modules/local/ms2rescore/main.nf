process MS2RESCORE {
    tag "${meta.id}_${engine}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ms2rescore:3.1.5--pyhdfd78af_0' :
        'biocontainers/ms2rescore:3.1.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(pin_file)
    val(engine)
    path(scans_pkls, stageAs: 'scans/*')
    path(mgf_files,  stageAs: 'scans/*')
    path(config_json)

    output:
    tuple val(meta), path("${engine}.psms.tsv"),   emit: psms, optional: true
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
        --pin    ${pin_file} \\
        --scans  ${scans_pkls.join(' ')} \\
        --mod    ${mod} \\
        --out    ${engine}_rescore_input.tsv

    # MS2Rescore's rescoring engine (mokapot/percolator) needs decoy PSMs to train
    # the FDR model. An engine that reports only confident hits against a small DB
    # can emit a decoy-free PIN, which aborts rescoring. Skip this engine instead
    # (no psms.tsv -> INTEGRATE proceeds with the engines that do have decoys).
    n_decoy=\$(awk -F'\\t' 'NR>1 && \$5=="True"{c++} END{print c+0}' ${engine}_rescore_input.tsv)
    if [ "\$n_decoy" -eq 0 ]; then
        echo "WARN: [MS2RESCORE_${engine}] ${meta.id}: 0 decoy PSMs in the PIN - skipping rescoring for this engine; INTEGRATE will use the remaining engines." | tee ${engine}_pipeline_log.txt >&2
    else
        # Apply user-configurable rescoring settings to the base config: rescoring
        # engine (mokapot|percolator), mokapot internal calibration FDR, and the
        # HTML QC report toggle. These override whatever the JSON asset ships with.
        cp ${config_json} run_config.json
        python3 -c "import json; c=json.load(open('run_config.json')); r=c['ms2rescore']; mk=r.get('rescoring_engine',{}).get('mokapot',{}); mk.update({'train_fdr': ${params.mokapot_train_fdr}, 'test_fdr': ${params.mokapot_test_fdr}}); r['rescoring_engine']=({'percolator':{}} if '${rescore_engine}'=='percolator' else {'mokapot': mk}); r['write_report']=${params.ms2rescore_write_report ? 'True' : 'False'}; json.dump(c, open('run_config.json','w'), indent=4)"

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
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ms2rescore: \$(ms2rescore --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
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
