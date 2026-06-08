process CONVERT_MZML {
    tag "${meta.id}_${mzml.baseName}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pymzml:2.5.2--pyhdfd78af_0' :
        'biocontainers/pymzml:2.5.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(mzml)

    output:
    tuple val(meta), path("*.mgf"),             emit: mgf
    tuple val(meta), path("*.scans.pkl"),       emit: scans
    tuple val(meta), path("*.index2scan.pkl"),  emit: index2scan
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    run = mzml.baseName
    template 'convert_mzml.py'

    stub:
    def run = mzml.baseName
    """
    touch ${run}.mgf ${run}.scans.pkl ${run}.index2scan.pkl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymzml: "stub"
    END_VERSIONS
    """
}
