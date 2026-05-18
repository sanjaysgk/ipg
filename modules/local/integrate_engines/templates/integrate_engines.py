#!/usr/bin/env python3
"""Integrate rescored PSMs across MS search engines at 1% peptide-level FDR.

Extracted from immunopeptidomics/core.py:
  - read_psms() (L1313)
  - count_chimera() (L1408)
  - read_results() (L1451) — aggregation + chimeric PSM resolution
  - pepform2seq() (L1276), sort_protein_list() (L1293)

Emits:
  integrated_psms.tsv      unique-peptide PSMs (post-chimera removal)
  integrated_peptides.tsv  per-peptide aggregated identifications
  chimeric_PSMs.txt        scans with multiple peptide assignments
  chimera_only_peptides.txt peptides only seen in chimeric PSMs
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO


OUT_COLUMNS = [
    "scan_id", "run", "peptide", "length", "immuno", "peptidoform",
    "protein_list", "gene", "species", "description leftmost",
    "engine", "psm_qval", "psm_score", "psm_PEP",
    "peptide_qval", "peptide_score", "peptide_PEP",
    "precursor_mz", "retention_time", "ion_mobility",
    "rescoring:spec_pearson", "rescoring:rt_diff_best",
    "rescoring:predicted_retention_time", "rescoring:observed_retention_time",
]

PEP_COLUMNS = [
    "peptide", "length", "immuno", "peptidoform", "protein_list",
    "gene", "species", "description leftmost", "engine",
    "peptide_qval", "peptide_score", "peptide_PEP",
    "rescoring:spec_pearson", "rescoring:rt_diff_best",
]


def pepform2seq(pepform: str) -> str:
    return re.sub(r"[^A-Z]", "", re.sub(r"\[.*?\]", "", pepform))


def sort_protein_list(protein_list: str) -> str:
    s = str(protein_list).strip("[]").replace("'", "").replace(", ", ";")
    return ";".join(sorted(p for p in s.split(";") if not p.startswith("rev_")))


def parse_fasta_prot_info(fasta: Path) -> dict:
    """Extract per-protein gene/species/description from FASTA headers."""
    info: dict = {}
    for rec in SeqIO.parse(str(fasta), "fasta"):
        header = rec.description
        gene = ""
        species = ""
        if "GN=" in header:
            m = re.search(r"GN=(\S+)", header)
            gene = m.group(1) if m else ""
        if "OS=" in header:
            m = re.search(r"OS=([^=]+?)(?=\s+[A-Z]{2}=|$)", header)
            species = m.group(1).strip() if m else ""
        info[rec.id] = {
            "gene": gene,
            "species": species,
            "description": header,
        }
    return info


def read_engine_psms(tsv: Path, engine: str, prot_info: dict,
                     immuno_lengths: list[int],
                     fdr: float = 0.01) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(tsv, sep="\t")
    # PSM- and peptide-level FDR, target only.
    df = df[
        (df["qvalue"] <= fdr)
        & (df["meta:peptide_qvalue"] <= fdr)
        & (~df["is_decoy"].astype(bool))
    ].copy()
    if df.empty:
        return (pd.DataFrame(columns=OUT_COLUMNS),
                pd.DataFrame(columns=PEP_COLUMNS))

    df["peptide"] = df["peptidoform"].apply(pepform2seq)
    df["length"] = df["peptide"].str.len()
    df["engine"] = engine
    df["scan_id"] = df["run"] + "_" + df["spectrum_id"].astype(str)

    df = df.sort_values(by=["qvalue", "pep"], ascending=[True, True])
    df = df.rename(columns={
        "qvalue": "psm_qval", "score": "psm_score", "pep": "psm_PEP",
        "meta:peptide_qvalue": "peptide_qval",
        "meta:peptide_score": "peptide_score",
        "meta:peptide_pep": "peptide_PEP",
    })

    df["protein_list"] = df["protein_list"].apply(sort_protein_list)
    df["gene"] = df["protein_list"].apply(
        lambda x: ";".join(prot_info.get(p, {}).get("gene", "") for p in x.split(";"))
    )
    df["species"] = df["protein_list"].apply(
        lambda x: ";".join(sorted(set(prot_info.get(p, {}).get("species", "") for p in x.split(";"))))
    )
    df["description leftmost"] = df["protein_list"].apply(
        lambda x: prot_info.get(x.split(";")[0], {}).get("description", "")
    )
    df["immuno"] = df["length"].isin(immuno_lengths)

    psms = df[[c for c in OUT_COLUMNS if c in df.columns]].copy()

    agg = {
        "length": "first", "immuno": "first",
        "peptidoform": lambda x: ";".join(pd.Series.unique(x)),
        "protein_list": "first", "gene": "first", "species": "first",
        "description leftmost": "first", "engine": "first",
        "peptide_qval": "min", "peptide_score": "max", "peptide_PEP": "min",
        "rescoring:spec_pearson": "max", "rescoring:rt_diff_best": "min",
    }
    agg = {k: v for k, v in agg.items() if k in psms.columns}
    peptides = psms.groupby("peptide").agg(agg).reset_index()
    return psms, peptides


def count_chimera(chimera_psms: pd.DataFrame) -> tuple[dict, dict]:
    count: dict = {}
    copep: dict = {}
    for peptides in chimera_psms["peptide"]:
        pep_list = peptides.split(";")
        for peptide in pep_list:
            count[peptide] = count.get(peptide, 0) + 1
            copep.setdefault(peptide, {})
            for p2 in pep_list:
                if p2 == peptide:
                    continue
                copep[peptide][p2] = copep[peptide].get(p2, 0) + 1
    for p in copep:
        copep[p] = [f"{q}({n})" for q, n in copep[p].items()]
    return count, copep


def integrate(engine_tsvs: dict[str, Path], fasta: Path,
              immuno_lengths: list[int], outdir: Path,
              fdr: float = 0.01) -> None:
    prot_info = parse_fasta_prot_info(fasta)
    psms_all = pd.DataFrame(columns=OUT_COLUMNS)
    peptides_all = pd.DataFrame(columns=PEP_COLUMNS)

    for engine, tsv in engine_tsvs.items():
        psms, peps = read_engine_psms(tsv, engine, prot_info, immuno_lengths, fdr=fdr)
        psms_all = pd.concat([psms_all, psms], ignore_index=True)
        peptides_all = pd.concat([peptides_all, peps], ignore_index=True)

    if psms_all.empty:
        print(f"[integrate_engines] no PSMs passed {fdr*100:.0f}% FDR in any engine",
              file=sys.stderr)
        psms_all.to_csv(outdir / "integrated_psms.tsv", sep="\t", index=False)
        peptides_all.to_csv(outdir / "integrated_peptides.tsv", sep="\t",
                            index=False)
        return

    # Aggregate PSMs across engines on scan_id.
    agg_psm = {
        "run": "first",
        "peptide": lambda x: ";".join(pd.Series.unique(x)),
        "length": "first", "immuno": "first",
        "peptidoform": list, "protein_list": "first",
        "gene": "first", "species": "first",
        "description leftmost": "first",
        "engine": list,
        "psm_qval": list, "psm_score": list, "psm_PEP": list,
        "peptide_qval": list, "peptide_score": list, "peptide_PEP": list,
        "precursor_mz": "first", "retention_time": "first",
        "ion_mobility": "first",
        "rescoring:spec_pearson": list, "rescoring:rt_diff_best": list,
        "rescoring:predicted_retention_time": list,
        "rescoring:observed_retention_time": list,
    }
    agg_psm = {k: v for k, v in agg_psm.items() if k in psms_all.columns}
    psm_agg = psms_all.groupby("scan_id").agg(agg_psm).reset_index()

    chimera = psm_agg[psm_agg["peptide"].str.contains(";")]
    chimera_count, chimera_pep = count_chimera(chimera)
    chimera.to_csv(outdir / "chimeric_PSMs.txt", sep="\t", index=False)
    psms = psm_agg[~psm_agg["peptide"].str.contains(";")].copy()
    print(f"[integrate_engines] unique PSMs: {len(psms)}  chimeric: {len(chimera)}",
          file=sys.stderr)

    # Aggregate peptides across engines.
    agg_pep = {
        "length": "first", "immuno": "first",
        "protein_list": "first", "gene": "first", "species": "first",
        "description leftmost": "first",
        "peptidoform": lambda x: ";".join(pd.Series.unique(x)),
        "engine": list,
        "peptide_qval": list, "peptide_score": list, "peptide_PEP": list,
        "rescoring:spec_pearson": "max", "rescoring:rt_diff_best": "min",
    }
    agg_pep = {k: v for k, v in agg_pep.items() if k in peptides_all.columns}
    peptides = peptides_all.groupby("peptide").agg(agg_pep).reset_index()

    peptides["peptidoform"] = peptides["peptidoform"].apply(
        lambda x: ";".join(set(x.split(";")))
    )
    peptides["#peptidoform"] = peptides["peptidoform"].str.count(";") + 1
    peptides["#engine"] = peptides["engine"].apply(len)
    peptides["peptide_qval_lowest"] = peptides["peptide_qval"].apply(min)
    peptides["peptide_PEP_lowest"] = peptides["peptide_PEP"].apply(min)
    peptides["peptide_score_highest"] = peptides["peptide_score"].apply(max)
    peptides["#chimeric_PSMs"] = (
        peptides["peptide"].map(chimera_count).fillna(0).astype(int)
    )
    peptides["#chimeric_peptides"] = (
        peptides["peptide"].map(chimera_pep).fillna("NA").astype(str)
    )

    psms_total = psms.groupby("peptide").size().reset_index(name="PSMs_total")
    peptides = peptides.merge(psms_total, on="peptide", how="left")
    psms_run = psms.groupby(["peptide", "run"]).size().unstack(fill_value=0)
    psms_run.columns = [f"PSMs_run_{c}" for c in psms_run.columns]
    peptides = peptides.merge(psms_run, left_on="peptide", right_index=True,
                              how="left")

    chim_only = peptides[
        (peptides["PSMs_total"] <= 0) | (peptides["PSMs_total"].isna())
    ]
    chim_only.to_csv(outdir / "chimera_only_peptides.txt", sep="\t", index=False)
    peptides = peptides[peptides["PSMs_total"] > 0]
    print(f"[integrate_engines] unique peptides: {len(peptides)}",
          file=sys.stderr)

    psms.to_csv(outdir / "integrated_psms.tsv", sep="\t", index=False)
    peptides.to_csv(outdir / "integrated_peptides.tsv", sep="\t", index=False)


def parse_engine_tsv_arg(values: list[str]) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for v in values:
        if "=" not in v:
            raise SystemExit(f"--engine-tsv expects name=path, got: {v}")
        name, path = v.split("=", 1)
        out[name] = Path(path)
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--engine-tsv", nargs="+", required=True,
                   help="one or more name=path pairs, e.g. msfragger=...tsv")
    p.add_argument("--fasta", required=True, type=Path,
                   help="search FASTA database (for protein info)")
    p.add_argument("--peptide-length", default="9",
                   help="single length or range like 8-11")
    p.add_argument("--fdr", type=float, default=0.01,
                   help="PSM and peptide FDR threshold (default: 0.01)")
    p.add_argument("--outdir", default=".", type=Path)
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    if "-" in args.peptide_length:
        a, b = args.peptide_length.split("-")
        immuno_lengths = list(range(int(a), int(b) + 1))
    else:
        immuno_lengths = [int(args.peptide_length)]

    engine_tsvs = parse_engine_tsv_arg(args.engine_tsv)
    integrate(engine_tsvs, args.fasta, immuno_lengths, args.outdir, fdr=args.fdr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
