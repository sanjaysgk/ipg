process INSTANOVO_PREDICT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    // De novo inference needs InstaNovo + torch (GPU-preferred). Under -profile pixi
    // the tool is provided by the isolated `denovo` pixi env (nextflow.config beforeScript).
    // TODO: GPU singularity image for -profile singularity / real-data runs.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(spectra)

    output:
    tuple val(meta), path("${prefix}.denovo.csv"), emit: predictions
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    // instanovo_predict.py is a thin shim around `instanovo` that forces
    // torch.load(weights_only=False) (instanovo hardcodes weights_only=True, which
    // torch>=2.6 rejects for the v1.2.0 checkpoint's defaultdict). InstaNovo reads
    // mzML/mzXML/MGF directly; output CSV carries the predictions_beam_* columns
    // Winnow consumes downstream.
    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    instanovo_predict.py transformer predict \\
        -d ${spectra} \\
        -o ${prefix}.denovo.csv \\
        --denovo ${args}

    ver=\$(python -c "import instanovo; print(instanovo.__version__)" 2>/dev/null || true)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        instanovo: \${ver:-unknown}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'experiment_name,scan_number,spectrum_id,precursor_mz,precursor_charge,predictions,log_probs,predictions_beam_0,predictions_log_probability_beam_0\\nstub,1,stub:1,500.0,2,PEPTIDEK,-0.1,PEPTIDEK,-0.1\\n' > ${prefix}.denovo.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        instanovo: "stub"
    END_VERSIONS
    """
}
