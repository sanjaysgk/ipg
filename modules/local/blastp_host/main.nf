process BLASTP_HOST {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::blast=2.15.0 conda-forge::python>=3.9 conda-forge::pandas>=1.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-bd0ce31f08f73d42bb0edebd8f5c7dcd9b1eefb0:1a671ad7f1b9f7bae925b3fe7baf49e06e59e70e-0' :
        'biocontainers/mulled-v2-bd0ce31f08f73d42bb0edebd8f5c7dcd9b1eefb0:1a671ad7f1b9f7bae925b3fe7baf49e06e59e70e-0' }"

    input:
    tuple val(meta), path(peptides_tsv)
    path(blast_db_files, stageAs: 'db/*')   // every file that shares --blast_db prefix
    val(db_prefix)                           // basename of the db, e.g. 'human'
    val(host_species)

    output:
    tuple val(meta), path("peptides.blastp.tsv"), emit: peptides
    tuple val(meta), path("blastp_out.tsv"),      emit: raw
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 - <<PY
import pandas as pd
df = pd.read_csv("${peptides_tsv}", sep="\\t")
mask = (~df["species"].astype(str).str.contains("${host_species}"))             \\
        & (~df["peptide"].astype(str).str.contains("contaminant"))
subset = df.loc[mask, "peptide"].dropna().unique()
with open("blastp_in.fasta", "w") as fh:
    for pep in subset:
        fh.write(f">{pep}\\n{str(pep).replace('I', 'L')}\\n")
PY

    if [ ! -s blastp_in.fasta ]; then
        echo "no non-host peptides; emitting unchanged table" >&2
        cp ${peptides_tsv} peptides.blastp.tsv
        : > blastp_out.tsv
    else
        blastp -task blastp-short \\
            -db db/${db_prefix} \\
            -query blastp_in.fasta \\
            -out blastp_out.tsv \\
            -outfmt "6 qacc sacc qlen nident sseq"

        python3 - <<PY
import pandas as pd
df = pd.read_csv("${peptides_tsv}", sep="\\t")
try:
    hits = pd.read_csv("blastp_out.tsv", sep="\\t",
                        names=["peptide", "BLASTP_match", "qlen", "nident", "BLASTP_matchedSeq"])
    hits["BLASTP_ident%"] = (hits["nident"] / hits["qlen"] * 100).round(2)
    top = (
        hits[["peptide", "BLASTP_ident%", "BLASTP_match", "BLASTP_matchedSeq"]]
        .sort_values("BLASTP_ident%", ascending=False)
        .groupby("peptide").first().reset_index()
    )
    merged = df.merge(top, on="peptide", how="left")
    merged["BLASTP_ident%"].fillna(0, inplace=True)
    merged["BLASTP_match"].fillna("NA", inplace=True)
    merged["BLASTP_matchedSeq"].fillna("NA", inplace=True)
    merged["species"] = merged.apply(
        lambda r: ";".join(sorted(set(str(r["species"]).split(";") +
                                        (["${host_species}"] if r["BLASTP_ident%"] == 100 else [])))),
        axis=1,
    )
except pd.errors.EmptyDataError:
    merged = df
merged.to_csv("peptides.blastp.tsv", sep="\\t", index=False)
PY
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastp -version | head -1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    cp ${peptides_tsv} peptides.blastp.tsv
    : > blastp_out.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: "stub"
    END_VERSIONS
    """
}
