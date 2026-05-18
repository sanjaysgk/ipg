process SAGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::sage-proteomics=0.14.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sage-proteomics:0.14.7--h92b4e83_0' :
        'biocontainers/sage-proteomics:0.14.7--h92b4e83_0' }"

    input:
    tuple val(meta), path(mzml_files)
    path(fasta)
    path(config_json)
    path(msfragger_log)  // optional; pass file('NO_FILE') if absent

    output:
    tuple val(meta), path("*.pin"),     emit: pin
    tuple val(meta), path("*.sage.tsv"), emit: tsv, optional: true
    path("sage_search_log.txt"),        emit: log
    path("used_params.json"),           emit: used_params
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def cpus = Math.max(1, (task.cpus ?: 2).intdiv(2))
    """
    # If an MSFragger search_log is present, re-use its calibrated fragment
    # tolerance and topN peaks (matches core.py run_Sage behaviour).
    ms2_tol=""
    max_peaks=""
    if [ -s "${msfragger_log}" ] && [ "${msfragger_log}" != "NO_FILE" ]; then
        ms2_tol=\$(grep 'New fragment_mass_tolerance' ${msfragger_log} | awk '{print \$4}' | tail -1 || true)
        max_peaks=\$(grep 'New use_topN_peaks'        ${msfragger_log} | awk '{print \$4}' | tail -1 || true)
    fi

    if [ -n "\$ms2_tol" ] && [ -n "\$max_peaks" ]; then
        python3 - <<PY
import json, sys
with open("${config_json}") as f:
    cfg = json.load(f)
cfg["fragment_tol"] = {"ppm": [-float("\$ms2_tol"), float("\$ms2_tol")]}
cfg["max_peaks"] = int("\$max_peaks")
with open("used_params.json", "w") as f:
    json.dump(cfg, f, indent=2)
PY
    else
        cp ${config_json} used_params.json
    fi

    export RUST_MIN_STACK=4194304
    sage \\
        --batch-size ${cpus} \\
        --write-pin \\
        -o . \\
        -f ${fasta} \\
        used_params.json \\
        ${mzml_files} \\
        > sage_search_log.txt 2>&1

    # Sage writes results.sage.pin — keep a clean name for mokapot
    if [ -f results.sage.pin ]; then
        mv results.sage.pin ${meta.id}.sage.pin
    fi

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    sage: \$(sage --version 2>&1 | awk '{print \$NF}' | head -1)
END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sage.pin
    touch sage_search_log.txt
    echo '{}' > used_params.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: "stub"
    END_VERSIONS
    """
}
