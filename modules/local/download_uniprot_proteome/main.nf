process DOWNLOAD_UNIPROT_PROTEOME {
    tag "${url.tokenize('/')[-1]}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    // Persistent cache: the proteome is fetched at most once and reused across
    // runs/work-dir resets. Needs network (same as conda/container pulls).
    storeDir "${params.outdir}/canonical_db"

    input:
    val(url)

    output:
    path("*.fasta"),     emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def gz   = url.tokenize('/')[-1]
    def base = gz.endsWith('.gz') ? gz[0..-4] : gz
    """
    wget --no-check-certificate -q "${url}" -O "${gz}"
    gunzip -f "${gz}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -1 | grep -oE '[0-9]+\\.[0-9]+(\\.[0-9]+)?' | head -1)
    END_VERSIONS
    """

    stub:
    def gz   = url.tokenize('/')[-1]
    def base = gz.endsWith('.gz') ? gz[0..-4] : gz
    """
    printf '>sp|STUB1|STUB_HUMAN stub canonical\\nPEPTIDESTUBKR\\n' > "${base}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: "stub"
    END_VERSIONS
    """
}
