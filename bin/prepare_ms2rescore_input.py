#!/usr/bin/env python3
"""Build an MS2Rescore TSV input from mokapot PSMs + PIN features + scans pickles.

Extracted from immunopeptidomics/core.py:
  - prepare_tsv() (L809)
  - get_peptidoform() (L749)
  - rescore_MSFragger / _Comet / _Sage / _PEAKS (L951-1117)

Engine-specific SpecId parsing is encapsulated in `_extract_run_specid()`.
"""
from __future__ import annotations

import argparse
import pickle
import re
import sys
from pathlib import Path

import pandas as pd

MOD_MAPPING = {
    "[+57.0215]": "[Carbamidomethyl]", "+57.021465": "Carbamidomethyl",
    "[+15.9949]": "[Oxidation]", "[+42.010567]": "[Acetyl]",
    "[+119.0041]": "[Cysteinylation]",
    "n[304.2072]": "[TMTpro]-", "n[+304.20715]": "[TMTpro]-",
    "[304.2072]": "[TMTpro]", "[+304.20715]": "[TMTpro]",
    "n[229.1629]": "[TMT6plex]-", "[229.1629]": "[TMT6plex]",
    "[+229.16293]": "[TMT6plex]", "(+229.16)": "[TMT6plex]",
    "n[+229.1629]": "[TMT6plex]-", "[+229.1629]": "[TMT6plex]",
    "[15.9949]": "[Oxidation]", "[57.0215]": "[Carbamidomethyl]",
    "[119.0041]": "[Cysteinylation]",
    "[-17.0265]": "[Gln->pyro-Glu]", "[-18.0106]": "[Glu->pyro-Glu]",
    "n[42.0106]": "[Acetyl]-",
    "[-17.026548]": "[Gln->pyro-Glu]", "[-18.010565]": "[Glu->pyro-Glu]",
    "(+15.99)": "[Oxidation]", "(+42.01)": "[Acetyl]",
    "(-17.03)": "[Gln->pyro-Glu]", "(-18.01)": "[Glu->pyro-Glu]",
    "(+119.00)": "[Cysteinylation]",
}


def get_peptidoform(peptide: str) -> str:
    for mod, psimod in MOD_MAPPING.items():
        peptide = peptide.replace(mod, psimod)
    for nt_mod in ("TMT6plex", "TMTpro", "Acetyl"):
        peptide = re.sub(rf"^([A-Z])\[({nt_mod})\]", r"[\2]-\1", peptide)
    return peptide


def _extract_run_specid(df: pd.DataFrame, engine: str) -> pd.DataFrame:
    """Engine-aware extraction of `run`, `specid`, and stripped `Peptide`."""
    if engine == "msfragger":
        df["run"] = df["SpecId"].apply(lambda x: ".".join(x.split(".")[:-3]))
        df["specid"] = df["run"] + "_" + df["ScanNr"].astype(str)
        df["Peptide"] = df["Peptide"].apply(lambda x: x[2:-3])
    elif engine == "comet":
        df["run"] = df["SpecId"].apply(
            lambda x: "_".join(Path(x).name.split("_")[:-3])
        )
        df["specid"] = df["run"] + "_" + df["ScanNr"].astype(str)
        df["Peptide"] = df["Peptide"].apply(lambda x: x[2:-2])
    elif engine == "sage":
        df["run"] = df["FileName"].apply(lambda x: Path(x).stem)
        df["specid"] = df["run"] + "_" + df["ScanNr"].astype(str)
    elif engine == "peaks":
        df["run"] = df["SpecId"].apply(lambda x: "_".join(x.split("_")[:-1]))
        df["specid"] = df["SpecId"]
    else:
        raise ValueError(f"unknown engine {engine}")
    return df


def prepare_tsv(psms: pd.DataFrame, scans: dict) -> pd.DataFrame:
    psms["z"] = psms["specid"].apply(lambda x: scans[x]["z"])
    tsv = pd.DataFrame()
    tsv["peptidoform"] = (
        psms["Peptide"].apply(get_peptidoform) + "/" + psms["z"]
    )
    tsv["peptide"] = psms["Peptide"].apply(lambda x: re.sub(r"[^A-Z]+", "", x))
    tsv["spectrum_id"] = (
        psms["specid"].apply(lambda x: scans[x]["rescore"]) + tsv["peptide"]
    )
    tsv["run"] = psms["run"]
    tsv["is_decoy"] = ~psms["Label"].astype(bool)
    tsv["score"] = psms["mokapot score"]
    tsv["qvalue"] = psms["mokapot q-value"]
    tsv["pep"] = psms["mokapot PEP"]
    tsv["ion_mobility"] = psms["specid"].apply(lambda x: scans[x]["im"])
    tsv["precursor_mz"] = psms["specid"].apply(lambda x: scans[x]["mz"])
    tsv["retention_time"] = psms["specid"].apply(lambda x: scans[x]["rt"])
    tsv["protein_list"] = psms["Proteins"].apply(lambda x: str(x).split(";"))
    tsv["rank"] = 1
    return tsv


def build_features(pin_path: Path, engine: str, scans: dict) -> pd.DataFrame:
    f = pd.read_csv(pin_path, sep="\t", dtype=str)
    f = _extract_run_specid(f, engine)
    # _extract_run_specid already stripped flanking residues for
    # msfragger/comet/peaks — just normalise to A-Z here. Sage leaves
    # Peptide raw, so the regex is the only stripping it gets.
    f["Peptide"] = f["Peptide"].apply(lambda x: re.sub(r"[^A-Z]+", "", x))
    f["spectrum_id"] = (
        f["specid"].apply(lambda x: scans[x]["rescore"]) + f["Peptide"]
    )
    drop_cols = ["SpecId", "Label", "Peptide", "Proteins", "run", "specid", "ScanNr"]
    if engine == "sage":
        drop_cols.append("FileName")
    if engine == "peaks" and "Delta RT" in f.columns:
        drop_cols.append("Delta RT")
    f = f.drop(columns=[c for c in drop_cols if c in f.columns])
    return f


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--engine", required=True,
                   choices=["msfragger", "comet", "sage", "peaks"])
    p.add_argument("--target", required=True, help="mokapot target PSMs TSV")
    p.add_argument("--decoy", required=True, help="mokapot decoy PSMs TSV")
    p.add_argument("--pin", required=True, help="combined PIN used by mokapot")
    p.add_argument("--scans", required=True, nargs="+",
                   help="scans.pkl file(s) from CONVERT_MZML")
    p.add_argument("--mod", default="mod",
                   help="modification profile (mod, nomod, TMT10, TMT16, mhcii, lowres)")
    p.add_argument("--out", required=True, help="output rescore_input.tsv")
    args = p.parse_args()

    # Merge scans dicts from all runs. Accept a mix of *.scans.pkl files
    # and directories containing such files — Nextflow stageAs sometimes
    # hands us one of each depending on channel composition.
    from pathlib import Path as _P
    pkl_paths: list[str] = []
    for s in args.scans:
        p = _P(s)
        if p.is_dir():
            pkl_paths.extend(str(x) for x in p.glob("*.scans.pkl"))
        elif p.is_file():
            pkl_paths.append(str(p))
    scans: dict = {}
    for pkl in pkl_paths:
        with open(pkl, "rb") as fh:
            scans.update(pickle.load(fh))
    print(f"[prepare_ms2rescore] merged {len(scans)} MS2 scans from "
          f"{len(pkl_paths)} pickle(s)", file=sys.stderr)

    target = pd.read_csv(args.target, sep="\t")
    decoy = pd.read_csv(args.decoy, sep="\t")
    psms = pd.concat([target, decoy], ignore_index=True)
    psms = _extract_run_specid(psms, args.engine)

    tsv = prepare_tsv(psms, scans)

    # Comet loses fixed TMT mods; reinstate them.
    if args.engine == "comet" and args.mod == "TMT16":
        tsv["peptidoform"] = tsv["peptidoform"].apply(
            lambda x: "[TMTpro]-" + x.replace("K", "K[TMTpro]")
        )
    elif args.engine == "comet" and args.mod == "TMT10":
        tsv["peptidoform"] = tsv["peptidoform"].apply(
            lambda x: "[TMT6plex]-" + x.replace("K", "K[TMT6plex]")
        )

    # Join percolator features (prefixed rescoring:).
    features = build_features(Path(args.pin), args.engine, scans)
    features = features.set_index("spectrum_id").add_prefix("rescoring:")
    if "rescoring:rank" in features.columns:
        features = features.drop(columns=["rescoring:rank"])

    tsv = tsv.set_index("spectrum_id", drop=False).join(features, how="inner")
    tsv.to_csv(args.out, sep="\t", index=False)
    print(f"[prepare_ms2rescore] wrote {len(tsv)} rows → {args.out}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
