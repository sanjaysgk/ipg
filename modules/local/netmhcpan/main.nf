process NETMHCPAN {
    tag "${meta.id}"
    label 'process_low'

    // netMHCpan-4.1 is closed-source, academic-license. Binary supplied by user.
    conda "conda-forge::python>=3.9 conda-forge::pandas>=1.5"

    input:
    tuple val(meta), path(peptides_tsv)
    val(hla_alleles)       // comma-separated list e.g. "HLA-A01:01,HLA-B07:02"
    path(netmhcpan_bin)    // user-provided path (staged into workdir)

    output:
    tuple val(meta), path("netmhcpan_output.txt"),    emit: raw
    tuple val(meta), path("netmhcpan_parsed.tsv"),    emit: parsed
    tuple val(meta), path("netmhcpan_best.tsv"),      emit: best
    path "versions.yml",                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Restrict to class-I lengths (8-12) before calling the binary.
    python3 - <<PY
import pandas as pd
df = pd.read_csv("${peptides_tsv}", sep="\\t")
df = df[df["length"].between(8, 12)]
df["peptide"].to_csv("input.txt", header=False, index=False)
PY

    if [ ! -s input.txt ]; then
        echo "no class-I length peptides; emitting empty output" >&2
        : > netmhcpan_output.txt
    else
        ${netmhcpan_bin} -p input.txt -l 8,9,10,11,12 -s \\
            -a '${hla_alleles}' > netmhcpan_output.txt
    fi

    parse_netmhcpan.py \\
        --output netmhcpan_output.txt \\
        --tool   netMHCpan \\
        --out-parsed netmhcpan_parsed.tsv \\
        --out-best   netmhcpan_best.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        netmhcpan: "4.1 (user-provided binary)"
    END_VERSIONS
    """

    stub:
    """
    touch netmhcpan_output.txt netmhcpan_parsed.tsv netmhcpan_best.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        netmhcpan: "stub"
    END_VERSIONS
    """
}
