process CONVERT_PEAKS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-bd0ce31f08f73d42bb0edebd8f5c7dcd9b1eefb0:1a671ad7f1b9f7bae925b3fe7baf49e06e59e70e-0' :
        'biocontainers/mulled-v2-bd0ce31f08f73d42bb0edebd8f5c7dcd9b1eefb0:1a671ad7f1b9f7bae925b3fe7baf49e06e59e70e-0' }"

    input:
    tuple val(meta), path(peaks_csv)
    path(index2scan_pkls, stageAs: 'idx/*')

    output:
    tuple val(meta), path("peaks.pin"), emit: pin
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    min_match = params.peaks_min_match_fraction ?: 0.98
    // Coerce to a list: with a single index2scan pickle Nextflow passes a scalar Path, and
    // .join(' ') would iterate its path components ('idx', 'name.pkl') instead of the file.
    idx_pkls = index2scan_pkls instanceof List ? index2scan_pkls : [index2scan_pkls]
    template 'convert_peaks_to_pin.py'

    stub:
    """
    touch peaks.pin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
