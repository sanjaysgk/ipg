process NETMHCIIPAN {
    tag "${meta.id}"
    label 'process_low'

    // netMHCIIpan-4.3 is closed-source, academic-license. Binary supplied by user.
    conda "conda-forge::python>=3.9 conda-forge::pandas>=1.5"

    input:
    tuple val(meta), path(peptides_tsv)
    val(hla_alleles)         // comma-separated class-II alleles
    path(netmhciipan_bin)    // user-provided path

    output:
    tuple val(meta), path("netmhciipan_output.txt"), emit: raw
    tuple val(meta), path("netmhciipan_parsed.tsv"), emit: parsed
    tuple val(meta), path("netmhciipan_best.tsv"),   emit: best
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Class-II lengths are 13-18.
    python3 - <<PY
import pandas as pd
df = pd.read_csv("${peptides_tsv}", sep="\\t")
df = df[df["length"].between(13, 18)]
df["peptide"].to_csv("input.txt", header=False, index=False)
PY

    if [ ! -s input.txt ]; then
        echo "no class-II length peptides; emitting empty output" >&2
        : > netmhciipan_output.txt
    else
        ${netmhciipan_bin} -inptype 1 -f input.txt -s \\
            -a '${hla_alleles}' > netmhciipan_output.txt
    fi

    parse_netmhcpan.py \\
        --output netmhciipan_output.txt \\
        --tool   netMHCIIpan \\
        --out-parsed netmhciipan_parsed.tsv \\
        --out-best   netmhciipan_best.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        netmhciipan: "4.3 (user-provided binary)"
    END_VERSIONS
    """

    stub:
    """
    touch netmhciipan_output.txt netmhciipan_parsed.tsv netmhciipan_best.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        netmhciipan: "stub"
    END_VERSIONS
    """
}
