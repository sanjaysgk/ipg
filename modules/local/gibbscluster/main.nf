process GIBBSCLUSTER {
    tag "${meta.id}"
    label 'process_low'

    conda "conda-forge::python>=3.9 conda-forge::pandas>=1.5 conda-forge::perl"

    input:
    tuple val(meta), path(peptides_tsv)
    path(gibbscluster_pl)   // user-provided path to GibbsCluster-2.0e_SA.pl

    output:
    tuple val(meta), path("gibbs_clusters.tsv"),                 emit: clusters
    tuple val(meta), path("gibbs_output/images/gibbs.KLDvsClusters.tab"), emit: kld
    tuple val(meta), path("gibbs_output/**"),                    emit: raw
    path "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def clusters = (params.gibbs_clusters == 'auto') ? '[2-5]' : params.gibbs_clusters
    """
    # Filter to immunopeptides only — mirrors core.py run_Gibbs L2028.
    python3 - <<PY
import pandas as pd
df = pd.read_csv("${peptides_tsv}", sep="\\t")
immuno = df[df.get("immuno", False).astype(bool)]
if "peptide" not in immuno.columns or len(immuno) == 0:
    raise SystemExit("no immunopeptides in input")
immuno["peptide"].to_csv("peptide_input.txt", header=False, index=False)
# Record length range so the shell below knows whether to add -C/-D/-I.
lens = sorted(set(immuno["length"]))
with open("_lens.txt", "w") as f:
    f.write(f"{min(lens)} {max(lens)} {len(lens)}\\n")
PY

    mkdir -p gibbs_output
    read min_len max_len n_lens < _lens.txt

    extra=""
    if [ "\$max_len" -lt 13 ]; then
        extra="-C"
        [ "\$n_lens" -gt 1 ] && extra="\$extra -D 4 -I 1"
    fi

    perl ${gibbscluster_pl} \\
        -H \$(which R) -T -j 2 -l 9 -k 32 -S 5 \\
        -P GibbsCluster2 \\
        -g '${clusters}' \\
        -f peptide_input.txt \\
        -R gibbs_output \\
        \$extra

    # GibbsCluster appends a trailing integer to its result directory.
    gibbs_subdir=\$(find gibbs_output -maxdepth 1 -type d -name 'GibbsCluster2*' | head -1)
    parse_gibbs.py \\
        --kld gibbs_output/images/gibbs.KLDvsClusters.tab \\
        --gibbs-dir "\$gibbs_subdir" \\
        --out gibbs_clusters.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gibbscluster: "2.0 (user-provided)"
        perl: \$(perl --version | awk '/\\(v/ {gsub(/[()v]/,""); print \$NF}' | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p gibbs_output/images gibbs_output/res
    touch gibbs_output/images/gibbs.KLDvsClusters.tab
    touch gibbs_clusters.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gibbscluster: "stub"
    END_VERSIONS
    """
}
