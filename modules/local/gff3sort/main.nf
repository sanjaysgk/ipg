process GFF3SORT {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::gff3sort=0.1.a1a2bc9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gff3sort:0.1.a1a2bc9--pl5321hdfd78af_1' :
        'biocontainers/gff3sort:0.1.a1a2bc9--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*_sorted.gtf"), emit: gtf
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when.every { it }

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${gtf.baseName}"
    def args   = task.ext.args   ?: "--chr_order original"
    """
    gff3sort.pl ${args} ${gtf} > ${prefix}_sorted.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gff3sort: \$(gff3sort.pl --version 2>&1 | grep -oP '[0-9][^ ]*' | head -1 || echo "0.1.a1a2bc9")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${gtf.baseName}"
    """
    touch ${prefix}_sorted.gtf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gff3sort: stub
    END_VERSIONS
    """
}
