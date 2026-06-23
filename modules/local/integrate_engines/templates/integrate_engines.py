#!/usr/bin/env python3
'''Integrate rescored PSMs across MS search engines at a configurable peptide-level FDR (default 1%).

Emits:
  integrated_psms.tsv      unique-peptide PSMs (post-chimera removal)
  integrated_peptides.tsv  per-peptide aggregated identifications
  chimeric_PSMs.txt        scans with multiple peptide assignments
  chimera_only_peptides.txt peptides only seen in chimeric PSMs
'''
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO


OUT_COLUMNS = [
    "scan_id", "run", "peptide", "length", "peptidoform",
    "protein_list", "protein_list_by_engine", "class", "class_detail",
    "gene", "species", "description leftmost",
    "engine", "psm_qval", "psm_score", "psm_PEP",
    "peptide_qval", "peptide_score", "peptide_PEP",
    "precursor_mz", "retention_time", "ion_mobility",
    "rescoring:spec_pearson", "rescoring:rt_diff_best",
    "rescoring:predicted_retention_time", "rescoring:observed_retention_time",
]

PEP_COLUMNS = [
    "peptide", "length", "peptidoform", "protein_list",
    "protein_list_by_engine", "class",
    "class_detail", "gene", "species", "description leftmost", "engine",
    "peptide_qval", "peptide_score", "peptide_PEP",
    "rescoring:spec_pearson", "rescoring:rt_diff_best",
]

# Per-class target-decoy q-values need enough target PSMs to be meaningful;
# below this the per-class FDR estimate is unreliable. Warned per engine per
# class, not enforced.
MIN_CLASS_TARGETS = 100

# Header-parsing patterns and options, configurable via params. Compiled once at
# the entry point (fail-fast) and read by the classify / parse helpers. Defaults
# reproduce the previous hardcoded behaviour. CRYPTIC_RE keys on the ORF suffix
# db_construct always appends (<transcript>_f<frame>p<orf>_<copy>), not a TCONS_
# prefix, so it holds for any input transcriptome. A protein matching neither
# pattern is class 'other' (contaminants, iRT) rather than silently 'cryptic'.
CANONICAL_RE = None
CRYPTIC_RE = None
GENE_RE = None
SPECIES_RE = None
FILTER_OTHER = False


def compile_regex(pattern: str, param: str):
    try:
        return re.compile(pattern)
    except re.error as exc:
        raise SystemExit(f"[integrate_engines] invalid {param} regex '{pattern}': {exc}")


def pepform2seq(pepform: str) -> str:
    return re.sub(r"[^A-Z]", "", re.sub(r"\\[.*?\\]", "", pepform))


def sort_protein_list(protein_list: str) -> str:
    s = str(protein_list).strip("[]").replace("'", "").replace(", ", ";")
    return ";".join(sorted(p for p in s.split(";") if not p.startswith("rev_")))


def _proteins_of(protein_list: str) -> list[str]:
    # Normalise the raw "['sp|..']" repr or a ;-joined list, drop the rev_ decoy
    # prefix so a decoy inherits its target's class (separated-class FDR needs it).
    s = str(protein_list).strip("[]").replace("'", "").replace(", ", ";")
    return [p[4:] if p.startswith("rev_") else p for p in s.split(";") if p]


def classify_protein(protein: str) -> str:
    # canonical if it matches the canonical pattern; else cryptic if it matches
    # the cryptic pattern; else 'other' (contaminant / iRT / unrecognised).
    if CANONICAL_RE.search(protein):
        return "canonical"
    if CRYPTIC_RE.search(protein):
        return "cryptic"
    return "other"


def classify_peptide(protein_list: str) -> str:
    '''FDR class with precedence canonical > cryptic > other. A peptide shared
    with a canonical protein is NOT a cryptic discovery, so canonical wins; a
    peptide on a cryptic ORF outranks an 'other' (contaminant) co-assignment.'''
    classes = {classify_protein(p) for p in _proteins_of(protein_list)}
    if "canonical" in classes:
        return "canonical"
    if "cryptic" in classes:
        return "cryptic"
    return "other"


def classify_peptide_detail(protein_list: str) -> str:
    '''REPORT label refining `class`:
      canonical_only  every mapped protein is canonical
      overlap         maps to a canonical AND a non-canonical protein
      cryptic         maps only to cryptic proteins (a true cryptic discovery)
      other           maps only to unrecognised proteins (contaminant / iRT)
    Mirrors classify_peptide (rev_-stripped) so it never disagrees on class.'''
    classes = {classify_protein(p) for p in _proteins_of(protein_list)}
    if "canonical" in classes:
        return "overlap" if len(classes) > 1 else "canonical_only"
    if "cryptic" in classes:
        return "cryptic"
    return "other"


def union_protein_lists(protein_lists) -> str:
    # Merge the ;-joined accessions from several rows into one sorted,
    # de-duplicated list. Used when engines (or PSMs of one peptide) report
    # different protein sets for the same identification, so the class is
    # derived from every protein any engine assigned rather than one engine's
    # row. A union only adds proteins, and classify_peptide is canonical-if-ANY,
    # so re-deriving on the union can move a peptide off the cryptic list but
    # never onto it. Non-string cells (e.g. a missing list) are skipped so they
    # never inject a bogus 'nan' accession.
    proteins = set()
    for value in protein_lists:
        if isinstance(value, str):
            proteins.update(p for p in value.split(";") if p)
    return ";".join(sorted(proteins))


def rederive_protein_fields(df, prot_info) -> None:
    # Recompute every protein-derived column from protein_list, in place. Run
    # after protein_list is merged across engines/PSMs so class, class_detail,
    # gene, species and the leftmost description reflect the merged protein set
    # rather than one engine's row. Reuses the same expressions as the per-PSM
    # pass, so only the protein set changes, not the formatting.
    df["class"] = df["protein_list"].apply(classify_peptide)
    df["class_detail"] = df["protein_list"].apply(classify_peptide_detail)
    # Drop blank lookups (cryptic ORFs carry no UniProt gene symbol) and
    # de-duplicate, so isoforms of one gene collapse to a single symbol and the
    # union-widened join is not padded with empty ';' separators.
    df["gene"] = df["protein_list"].apply(
        lambda x: ";".join(sorted({
            g for g in (prot_info.get(p, {}).get("gene", "") for p in x.split(";")) if g})))
    df["species"] = df["protein_list"].apply(
        lambda x: ";".join(sorted({
            s for s in (prot_info.get(p, {}).get("species", "") for p in x.split(";")) if s})))
    df["description leftmost"] = df["protein_list"].apply(
        lambda x: prot_info.get(x.split(";")[0], {}).get("description", "") if x else "")


def build_protein_list_by_engine(engines, protein_lists) -> str:
    # JSON {engine: [sorted proteins]} preserving each engine's own protein list
    # for a merged identification, so a cross-engine disagreement stays auditable
    # after the union. An engine contributing several rows is unioned (dict of
    # sets) so nothing is dropped; keys and lists are sorted for deterministic
    # output. The engine key is kept even when its list is missing; a non-string
    # value just contributes no proteins (never a bogus 'nan').
    by_engine: dict = {}
    for engine, value in zip(engines, protein_lists):
        proteins = by_engine.setdefault(str(engine), set())
        if isinstance(value, str):
            proteins.update(p for p in value.split(";") if p)
    return json.dumps(
        {e: sorted(v) for e, v in sorted(by_engine.items())},
        separators=(",", ":"),
    )


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
    '''Extract per-protein gene/species/description from FASTA headers using the
    configurable GENE_RE / SPECIES_RE patterns (capture group 1).'''
    info: dict = {}
    for rec in SeqIO.parse(str(fasta), "fasta"):
        header = rec.description
        gm = GENE_RE.search(header)
        sm = SPECIES_RE.search(header)
        info[rec.id] = {
            "gene": gm.group(1) if gm else "",
            "species": sm.group(1).strip() if sm else "",
            "description": header,
        }
    return info


def read_engine_psms(tsv: Path, engine: str, prot_info: dict,
                     fdr: float = 0.01) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(tsv, sep="\\t")
    if df.empty:
        return (pd.DataFrame(columns=OUT_COLUMNS),
                pd.DataFrame(columns=PEP_COLUMNS))

    # Classify EVERY PSM (targets + decoys); a decoy inherits its target's class
    # via rev_-stripping so each class gets its own target-decoy competition.
    df["class"] = df["protein_list"].apply(classify_peptide)
    df["class_detail"] = df["protein_list"].apply(classify_peptide_detail)
    # 'other' = contaminant / iRT / unrecognised. Dropped before FDR only when
    # --filter_other is set; default keeps them (annotate-only) so a #CONTAM DB
    # carrying real hits is never silently removed.
    if FILTER_OTHER:
        df = df[df["class"] != "other"].copy()
        if df.empty:
            return (pd.DataFrame(columns=OUT_COLUMNS),
                    pd.DataFrame(columns=PEP_COLUMNS))
    df["_decoy"] = df["is_decoy"].astype(bool)
    df["peptide"] = df["peptidoform"].apply(pepform2seq)
    df = df.rename(columns={
        "score": "psm_score", "pep": "psm_PEP",
        "meta:peptide_score": "peptide_score",
        "meta:peptide_pep": "peptide_PEP",
    })

    # Warn when a class is too small for a reliable per-class TDC estimate.
    target_class = df.loc[~df["_decoy"], "class"].value_counts()
    for _cls in ("canonical", "cryptic", "other"):
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
    psms = df[[c for c in OUT_COLUMNS if c in df.columns]].copy()

    # Within one engine a peptide can appear in several PSMs that disagree on the
    # reported protein set; union them so the per-engine protein_list is complete
    # before engines are merged. Protein-derived fields are recomputed from the
    # union (rederive_protein_fields), so they are not aggregated here.
    agg = {
        "length": "first",
        "peptidoform": lambda x: ";".join(pd.Series.unique(x)),
        "protein_list": union_protein_lists, "engine": "first",
        "peptide_qval": "min", "peptide_score": "max", "peptide_PEP": "min",
        "rescoring:spec_pearson": "max", "rescoring:rt_diff_best": "min",
    }
    agg = {k: v for k, v in agg.items() if k in psms.columns}
    peptides = psms.groupby("peptide").agg(agg).reset_index()
    rederive_protein_fields(peptides, prot_info)
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
              outdir: Path, fdr: float = 0.01,
              min_spec_pearson=None, max_rt_diff=None) -> None:
    prot_info = parse_fasta_prot_info(fasta)
    psms_all = pd.DataFrame(columns=OUT_COLUMNS)
    peptides_all = pd.DataFrame(columns=PEP_COLUMNS)

    for engine, tsv in engine_tsvs.items():
        psms, peps = read_engine_psms(tsv, engine, prot_info, fdr=fdr)
        psms_all = pd.concat([psms_all, psms], ignore_index=True)
        peptides_all = pd.concat([peptides_all, peps], ignore_index=True)

    if psms_all.empty:
        print(f"[integrate_engines] no PSMs passed {fdr*100:.0f}% FDR in any engine",
              file=sys.stderr)
        psms_all.to_csv(outdir / "integrated_psms.tsv", sep="\\t", index=False)
        peptides_all.to_csv(outdir / "integrated_peptides.tsv", sep="\\t",
                            index=False)
        return

    # Aggregate PSMs across engines on scan_id. protein_list is unioned across
    # engines; _pl_raw keeps each engine's own list (aligned with engine: list)
    # so protein_list_by_engine can preserve the per-engine breakdown. Protein-
    # derived fields are recomputed from the union afterwards.
    psms_all["_pl_raw"] = psms_all["protein_list"]
    agg_psm = {
        "run": "first",
        "peptide": lambda x: ";".join(pd.Series.unique(x)),
        "length": "first",
        "peptidoform": list, "protein_list": union_protein_lists,
        "engine": list, "_pl_raw": list,
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
    rederive_protein_fields(psm_agg, prot_info)
    psm_agg["protein_list_by_engine"] = [
        build_protein_list_by_engine(eng, raw)
        for eng, raw in zip(psm_agg["engine"], psm_agg["_pl_raw"])
    ]
    psm_agg = psm_agg.drop(columns="_pl_raw")

    chimera = psm_agg[psm_agg["peptide"].str.contains(";")]
    chimera_count, chimera_pep = count_chimera(chimera)
    chimera.to_csv(outdir / "chimeric_PSMs.txt", sep="\\t", index=False)
    psms = psm_agg[~psm_agg["peptide"].str.contains(";")].copy()
    print(f"[integrate_engines] unique PSMs: {len(psms)}  chimeric: {len(chimera)}",
          file=sys.stderr)

    # Aggregate peptides across engines. protein_list is unioned; _pl_raw keeps
    # each engine's own list (aligned with engine: list). Protein-derived fields
    # are recomputed from the union, and protein_list_by_engine preserves the raw
    # per-engine lists.
    peptides_all["_pl_raw"] = peptides_all["protein_list"]
    agg_pep = {
        "length": "first",
        "protein_list": union_protein_lists,
        "peptidoform": lambda x: ";".join(pd.Series.unique(x)),
        "engine": list, "_pl_raw": list,
        "peptide_qval": list, "peptide_score": list, "peptide_PEP": list,
        "rescoring:spec_pearson": "max", "rescoring:rt_diff_best": "min",
    }
    agg_pep = {k: v for k, v in agg_pep.items() if k in peptides_all.columns}
    peptides = peptides_all.groupby("peptide").agg(agg_pep).reset_index()
    rederive_protein_fields(peptides, prot_info)
    peptides["protein_list_by_engine"] = [
        build_protein_list_by_engine(eng, raw)
        for eng, raw in zip(peptides["engine"], peptides["_pl_raw"])
    ]
    peptides = peptides.drop(columns="_pl_raw")

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
FDR = float("${fdr_arg}")
MIN_SPEC_PEARSON = float("${min_spec_pearson}") if "${min_spec_pearson}" else None
MAX_RT_DIFF = float("${max_rt_diff}") if "${max_rt_diff}" else None
FILTER_OTHER = "${filter_other}" == "true"
PROCESS_NAME = "${task.process}"

# Header-parsing patterns. An explicit param value (interpolated) wins; otherwise
# the default raw-string regex is used. Defaults are written with doubled
# backslashes / escaped '\$' so the Nextflow template engine yields the intended
# pattern. Canonical: explicit regex wins, else compile the prefix list (sugar).
_canon_regex = "${canonical_id_regex}"
if _canon_regex:
    CANONICAL_RE = compile_regex(_canon_regex, "canonical_id_regex")
else:
    _prefixes = [re.escape(p) for p in "${canonical_prefixes}".split(",") if p]
    _canon_default = ("^(?:" + "|".join(_prefixes) + ")") if _prefixes else r"(?!x)x"
    CANONICAL_RE = compile_regex(_canon_default, "canonical_protein_prefixes")
CRYPTIC_RE = compile_regex("${cryptic_id_regex}" or r"_f\\d+p\\d+_\\d+\$", "cryptic_id_regex")
GENE_RE = compile_regex("${gene_regex}" or r"GN=(\\S+)", "gene_regex")
SPECIES_RE = compile_regex("${species_regex}" or r"OS=([^=]+?)(?=\\s+[A-Z]{2}=|\$)", "species_regex")

engine_tsvs = parse_engine_tsv_arg(ENGINE_TSV_PAIRS.split())
integrate(engine_tsvs, Path(FASTA), Path("."), fdr=FDR,
          min_spec_pearson=MIN_SPEC_PEARSON, max_rt_diff=MAX_RT_DIFF)

with open("versions.yml", "w") as _f:
    _f.write(f'"{PROCESS_NAME}":\\n')
    _f.write(f"    python: {sys.version.split()[0]}\\n")
    _f.write(f"    pandas: {pd.__version__}\\n")
    _f.write(f"    biopython: {Bio.__version__}\\n")
