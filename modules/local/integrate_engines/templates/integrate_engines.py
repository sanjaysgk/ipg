#!/usr/bin/env python3
'''Integrate rescored PSMs across MS search engines at 1% peptide-level FDR.

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
'''
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO


OUT_COLUMNS = [
    "scan_id", "run", "peptide", "length", "immuno", "peptidoform",
    "protein_list", "class", "class_detail", "gene", "species",
    "description leftmost",
    "engine", "psm_qval", "psm_score", "psm_PEP",
    "peptide_qval", "peptide_score", "peptide_PEP",
    "precursor_mz", "retention_time", "ion_mobility",
    "rescoring:spec_pearson", "rescoring:rt_diff_best",
    "rescoring:predicted_retention_time", "rescoring:observed_retention_time",
]

PEP_COLUMNS = [
    "peptide", "length", "immuno", "peptidoform", "protein_list", "class",
    "class_detail", "gene", "species", "description leftmost", "engine",
    "peptide_qval", "peptide_score", "peptide_PEP",
    "rescoring:spec_pearson", "rescoring:rt_diff_best",
]

# Per-class target-decoy q-values need enough target PSMs to be meaningful;
# below this the per-class FDR estimate is unreliable. Warned per engine per
# class, not enforced.
MIN_CLASS_TARGETS = 100


def pepform2seq(pepform: str) -> str:
    return re.sub(r"[^A-Z]", "", re.sub(r"\\[.*?\\]", "", pepform))


def sort_protein_list(protein_list: str) -> str:
    s = str(protein_list).strip("[]").replace("'", "").replace(", ", ";")
    return ";".join(sorted(p for p in s.split(";") if not p.startswith("rev_")))


def classify_peptide(protein_list: str, canonical_prefixes: list[str]) -> str:
    '''canonical if ANY mapped protein is canonical (UniProt-style accession);
    cryptic only when every protein is non-canonical. A peptide shared with a
    canonical protein is NOT a cryptic discovery, so the canonical label wins.
    Strips the rev_ decoy prefix first so a decoy inherits its target's class
    (rev_sp|.. -> canonical, rev_TCONS.. -> cryptic), which separated-class FDR
    relies on. Accepts the raw "['sp|..']" repr or a ;-joined accession list.'''
    s = str(protein_list).strip("[]").replace("'", "").replace(", ", ";")
    proteins = [p[4:] if p.startswith("rev_") else p for p in s.split(";") if p]
    if any(p.startswith(pre) for p in proteins for pre in canonical_prefixes):
        return "canonical"
    return "cryptic"


def classify_peptide_detail(protein_list: str, canonical_prefixes: list[str]) -> str:
    '''Three-way REPORT label refining the FDR `class` (which stays canonical):
      canonical_only  every mapped protein is canonical
      overlap         maps to a canonical AND a cryptic protein
      cryptic         maps only to cryptic proteins (a true cryptic discovery)
    A peptide shared with a canonical protein is reported as overlap, never as a
    cryptic discovery, so the cryptic headline counts only novel sequences.
    Mirrors classify_peptide (rev_-stripped) so it never disagrees on class.'''
    s = str(protein_list).strip("[]").replace("'", "").replace(", ", ";")
    proteins = [p[4:] if p.startswith("rev_") else p for p in s.split(";") if p]
    is_canon = [any(p.startswith(pre) for pre in canonical_prefixes) for p in proteins]
    if any(is_canon) and not all(is_canon):
        return "overlap"
    if any(is_canon):
        return "canonical_only"
    return "cryptic"


def _tdc_q(score, is_decoy):
    '''Monotone target-decoy q-values for one competition (high score = better).
    q(s) = min over t<=s of (#decoys>=t)/(#targets>=t); matches mokapot semantics.'''
    s = np.asarray(score, dtype=float)
    d = np.asarray(is_decoy, dtype=bool)
    order = np.argsort(-s, kind="mergesort")
    d_sorted = d[order]
    fdr = np.cumsum(d_sorted) / np.maximum(np.cumsum(~d_sorted), 1)
    q_sorted = np.minimum.accumulate(fdr[::-1])[::-1]
    out = np.empty_like(q_sorted)
    out[order] = q_sorted
    return out


def class_separated_psm_q(df: pd.DataFrame, score_col: str) -> pd.Series:
    '''PSM-level TDC q estimated independently within each class group, so
    cryptic targets compete only against cryptic decoys (never canonical).'''
    q = pd.Series(1.0, index=df.index, dtype=float)
    for _, idx in df.groupby("class").groups.items():
        sub = df.loc[idx]
        q.loc[idx] = _tdc_q(sub[score_col].to_numpy(), sub["_decoy"].to_numpy())
    return q


def class_separated_peptide_q(df: pd.DataFrame):
    '''Peptide-level TDC q within class: best peptide_score per
    (class, peptide, decoy), competed within class, mapped back to every PSM.'''
    rep = df.groupby(["class", "peptide", "_decoy"], as_index=False)["peptide_score"].max()
    rep["pq"] = 1.0
    for _, idx in rep.groupby("class").groups.items():
        sub = rep.loc[idx]
        rep.loc[idx, "pq"] = _tdc_q(sub["peptide_score"].to_numpy(), sub["_decoy"].to_numpy())
    merged = df.merge(rep[["class", "peptide", "_decoy", "pq"]],
                      on=["class", "peptide", "_decoy"], how="left")
    return merged["pq"].to_numpy()


def parse_fasta_prot_info(fasta: Path) -> dict:
    '''Extract per-protein gene/species/description from FASTA headers.'''
    info: dict = {}
    for rec in SeqIO.parse(str(fasta), "fasta"):
        header = rec.description
        gene = ""
        species = ""
        if "GN=" in header:
            m = re.search(r"GN=(\\S+)", header)
            gene = m.group(1) if m else ""
        if "OS=" in header:
            m = re.search(r"OS=([^=]+?)(?=\\s+[A-Z]{2}=|\$)", header)
            species = m.group(1).strip() if m else ""
        info[rec.id] = {
            "gene": gene,
            "species": species,
            "description": header,
        }
    return info


def read_engine_psms(tsv: Path, engine: str, prot_info: dict,
                     immuno_lengths: list[int], canonical_prefixes: list[str],
                     fdr: float = 0.01) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(tsv, sep="\\t")
    if df.empty:
        return (pd.DataFrame(columns=OUT_COLUMNS),
                pd.DataFrame(columns=PEP_COLUMNS))

    # Classify EVERY PSM (targets + decoys); a decoy inherits its target's class
    # via rev_-stripping so each class gets its own target-decoy competition.
    df["class"] = df["protein_list"].apply(
        lambda x: classify_peptide(x, canonical_prefixes)
    )
    df["class_detail"] = df["protein_list"].apply(
        lambda x: classify_peptide_detail(x, canonical_prefixes)
    )
    df["_decoy"] = df["is_decoy"].astype(bool)
    df["peptide"] = df["peptidoform"].apply(pepform2seq)
    df = df.rename(columns={
        "score": "psm_score", "pep": "psm_PEP",
        "meta:peptide_score": "peptide_score",
        "meta:peptide_pep": "peptide_PEP",
    })

    # Warn when a class is too small for a reliable per-class TDC estimate.
    target_class = df.loc[~df["_decoy"], "class"].value_counts()
    for _cls in ("canonical", "cryptic"):
        _n = int(target_class.get(_cls, 0))
        if 0 < _n < MIN_CLASS_TARGETS:
            print(f"[integrate_engines] WARNING: {engine} {_cls} class has {_n} "
                  f"target PSMs (<{MIN_CLASS_TARGETS}); per-class {fdr*100:.0f}% FDR "
                  f"may be unreliable", file=sys.stderr)

    # Separated-class FDR: cryptic targets compete only against cryptic decoys,
    # canonical only against canonical. Replaces the canonical-dominated global q
    # so the cryptic subset is not swamped by the canonical majority.
    df["psm_qval"] = class_separated_psm_q(df, "psm_score")
    df["peptide_qval"] = class_separated_peptide_q(df)

    # Class-specific PSM- and peptide-level FDR cut, target only.
    df = df[
        (df["psm_qval"] <= fdr)
        & (df["peptide_qval"] <= fdr)
        & (~df["_decoy"])
    ].copy()
    if df.empty:
        return (pd.DataFrame(columns=OUT_COLUMNS),
                pd.DataFrame(columns=PEP_COLUMNS))

    df["length"] = df["peptide"].str.len()
    df["engine"] = engine
    df["scan_id"] = df["run"] + "_" + df["spectrum_id"].astype(str)
    df = df.sort_values(by=["psm_qval", "psm_PEP"], ascending=[True, True])

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
        "protein_list": "first", "class": "first", "class_detail": "first",
        "gene": "first", "species": "first",
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


def cryptic_spectral_status(row, min_spec_pearson, max_rt_diff) -> str:
    '''Class-specific orthogonal-evidence gate. Cryptic peptides lack database
    corroboration, so MS2PIP fragment correlation and DeepLC RT agreement are
    independent evidence the PSM is real. Canonical peptides are not gated; when
    no threshold is configured cryptic peptides are flagged 'not_filtered'.'''
    if str(row.get("class", "")) != "cryptic":
        return ""
    if min_spec_pearson is None and max_rt_diff is None:
        return "not_filtered"
    sp = row.get("rescoring:spec_pearson")
    rt = row.get("rescoring:rt_diff_best")
    if min_spec_pearson is not None and (pd.isna(sp) or sp < min_spec_pearson):
        return "fail"
    if max_rt_diff is not None and (pd.isna(rt) or rt > max_rt_diff):
        return "fail"
    return "pass"


def integrate(engine_tsvs: dict[str, Path], fasta: Path,
              immuno_lengths: list[int], canonical_prefixes: list[str],
              outdir: Path, fdr: float = 0.01,
              min_spec_pearson=None, max_rt_diff=None) -> None:
    prot_info = parse_fasta_prot_info(fasta)
    psms_all = pd.DataFrame(columns=OUT_COLUMNS)
    peptides_all = pd.DataFrame(columns=PEP_COLUMNS)

    for engine, tsv in engine_tsvs.items():
        psms, peps = read_engine_psms(tsv, engine, prot_info, immuno_lengths,
                                      canonical_prefixes, fdr=fdr)
        psms_all = pd.concat([psms_all, psms], ignore_index=True)
        peptides_all = pd.concat([peptides_all, peps], ignore_index=True)

    if psms_all.empty:
        print(f"[integrate_engines] no PSMs passed {fdr*100:.0f}% FDR in any engine",
              file=sys.stderr)
        psms_all.to_csv(outdir / "integrated_psms.tsv", sep="\\t", index=False)
        peptides_all.to_csv(outdir / "integrated_peptides.tsv", sep="\\t",
                            index=False)
        return

    # Aggregate PSMs across engines on scan_id.
    agg_psm = {
        "run": "first",
        "peptide": lambda x: ";".join(pd.Series.unique(x)),
        "length": "first", "immuno": "first",
        "peptidoform": list, "protein_list": "first", "class": "first",
        "class_detail": "first",
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
    chimera.to_csv(outdir / "chimeric_PSMs.txt", sep="\\t", index=False)
    psms = psm_agg[~psm_agg["peptide"].str.contains(";")].copy()
    print(f"[integrate_engines] unique PSMs: {len(psms)}  chimeric: {len(chimera)}",
          file=sys.stderr)

    # Aggregate peptides across engines.
    agg_pep = {
        "length": "first", "immuno": "first",
        "protein_list": "first", "class": "first", "class_detail": "first",
        "gene": "first", "species": "first",
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
    peptides["spectral_status"] = peptides.apply(
        lambda r: cryptic_spectral_status(r, min_spec_pearson, max_rt_diff),
        axis=1,
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
    chim_only.to_csv(outdir / "chimera_only_peptides.txt", sep="\\t", index=False)
    peptides = peptides[peptides["PSMs_total"] > 0]
    print(f"[integrate_engines] unique peptides: {len(peptides)}",
          file=sys.stderr)
    print(f"[integrate_engines] peptides by class_detail: "
          f"{peptides['class_detail'].value_counts().to_dict()}", file=sys.stderr)

    psms.to_csv(outdir / "integrated_psms.tsv", sep="\\t", index=False)
    peptides.to_csv(outdir / "integrated_peptides.tsv", sep="\\t", index=False)


def parse_engine_tsv_arg(values: list[str]) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for v in values:
        if "=" not in v:
            raise SystemExit(f"--engine-tsv expects name=path, got: {v}")
        name, path = v.split("=", 1)
        out[name] = Path(path)
    return out


# --- Nextflow template entry point ---------------------------------
# Variables interpolated by Nextflow at process execution time.
import pandas as pd  # already imported above; ensures versions block has it
import Bio

ENGINE_TSV_PAIRS = "${pairs}"   # space-separated "name=path" tokens
FASTA = "${fasta}"
PEPTIDE_LENGTH = "${len_arg}"
FDR = float("${fdr_arg}")
CANONICAL_PREFIXES = [p for p in "${canonical_prefixes}".split(",") if p]
MIN_SPEC_PEARSON = float("${min_spec_pearson}") if "${min_spec_pearson}" else None
MAX_RT_DIFF = float("${max_rt_diff}") if "${max_rt_diff}" else None
PROCESS_NAME = "${task.process}"

if "-" in PEPTIDE_LENGTH:
    _a, _b = PEPTIDE_LENGTH.split("-")
    immuno_lengths = list(range(int(_a), int(_b) + 1))
else:
    immuno_lengths = [int(PEPTIDE_LENGTH)]

engine_tsvs = parse_engine_tsv_arg(ENGINE_TSV_PAIRS.split())
integrate(engine_tsvs, Path(FASTA), immuno_lengths, CANONICAL_PREFIXES,
          Path("."), fdr=FDR,
          min_spec_pearson=MIN_SPEC_PEARSON, max_rt_diff=MAX_RT_DIFF)

with open("versions.yml", "w") as _f:
    _f.write(f'"{PROCESS_NAME}":\\n')
    _f.write(f"    python: {sys.version.split()[0]}\\n")
    _f.write(f"    pandas: {pd.__version__}\\n")
    _f.write(f"    biopython: {Bio.__version__}\\n")
