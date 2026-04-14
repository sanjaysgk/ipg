process FLASHLFQ {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::flashlfq=2.1.4 conda-forge::python>=3.9 conda-forge::pandas>=1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flashlfq:2.1.4--ha8f3691_0' :
        'biocontainers/flashlfq:2.1.4--ha8f3691_0' }"

    input:
    tuple val(meta), path(peptides_tsv), path(ms_files)

    output:
    tuple val(meta), path("quant/QuantifiedPeptides.tsv"),    emit: peptides
    tuple val(meta), path("quant/QuantifiedPeaks.tsv"),       emit: peaks, optional: true
    tuple val(meta), path("quant/flashlfq_input.tsv"),        emit: input
    path("quant/quant_log.txt"),                              emit: log
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mbr = (ms_files instanceof List && ms_files.size() > 1) ? 'true' : 'false'
    def cpus = task.cpus ?: 4
    """
    mkdir -p quant
    prepare_flashlfq_input.py \\
        --peptides ${peptides_tsv} \\
        --out quant/flashlfq_input.tsv

    FlashLFQ \\
        --idt quant/flashlfq_input.tsv \\
        --thr ${cpus} \\
        --chg --ppm 5 --ath \\
        --rep quant \\
        --mbr ${mbr} \\
        > quant/quant_log.txt 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flashlfq: \$(FlashLFQ --help 2>&1 | grep -oE 'FlashLFQ v[^ ]+' | head -1 | sed 's/FlashLFQ v//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p quant
    touch quant/QuantifiedPeptides.tsv quant/flashlfq_input.tsv quant/quant_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flashlfq: "stub"
    END_VERSIONS
    """
}
