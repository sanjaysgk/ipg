process PEPQUERY {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pepquery:2.0.2--hdfd78af_0' :
        'biocontainers/pepquery:2.0.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(peptides), path(mgf, stageAs: 'spectra/*'), path(ref_fasta)

    output:
    tuple val(meta), path("pepquery_out/psm_rank.txt"), emit: psm_rank, optional: true
    tuple val(meta), path("pepquery_out"),              emit: outdir
    path "versions.yml",                                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Non-specific digestion (-e 0 -c 0) is required for immunopeptides; the
    // tolerance/modification defaults match the orbitrap config and can be
    // overridden per run via task.ext.args (e.g. instrument-specific tols).
    def args = task.ext.args ?: '-e 0 -c 0 -fixMod 0 -varMod 2 -tol 10 -tolu ppm -itol 0.02 -fast'
    def mem  = task.memory ? "-Xmx${task.memory.toGiga()}g" : '-Xmx8g'
    """
    # PepQuery reads MGF reliably (mzML triggers SQLite index errors); concat the
    # per-sample MGFs from CONVERT_MZML into one file for a single -ms input.
    cat spectra/*.mgf > combined.mgf

    pepquery ${mem} \\
        -i ${peptides} \\
        -t peptide \\
        -s 1 \\
        -ms combined.mgf \\
        -db ${ref_fasta} \\
        -cpu ${task.cpus} \\
        ${args} \\
        -o pepquery_out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pepquery: "2.0.2"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p pepquery_out
    touch pepquery_out/psm_rank.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pepquery: "stub"
    END_VERSIONS
    """
}
