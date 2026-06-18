process PICKED_GROUP_FDR {
    tag "${meta.id}_${engine}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), val(engine), path(mokapot_psms), path(mokapot_decoy), path(cryptic_fasta)

    output:
    tuple val(meta), val(engine), path("*_protein_groups.tsv"), emit: protein_groups, optional: true
    path "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix             = "${meta.id}_${engine}"
    fdr                = params.protein_group_fdr ?: '0.01'
    fdr_class          = params.protein_group_fdr_class ?: 'global'
    canonical_prefixes = params.canonical_protein_prefixes ?: 'sp|,tr|'
    template 'picked_group_fdr.py'

    stub:
    prefix = "${meta.id}_${engine}"
    """
    touch ${prefix}_protein_groups.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picked_group_fdr: "stub"
    END_VERSIONS
    """
}
