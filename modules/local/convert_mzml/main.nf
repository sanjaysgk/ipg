process CONVERT_MZML {
    tag "${meta.id}_${mzml.baseName}"
    label 'process_low'

    conda "bioconda::pymzml=2.5.2"
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
    def run = mzml.baseName
    """
    convert_mzml.py --run ${run} ${mzml}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pymzml: \$(python3 -c "import pymzml; print(pymzml.__version__)")
    END_VERSIONS
    """

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
