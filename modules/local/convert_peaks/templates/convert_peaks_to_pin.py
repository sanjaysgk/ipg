#!/usr/bin/env python3
"""Convert a PEAKS Studio db.psms.csv export into a Percolator input (PIN).

Extracted from immunopeptidomics/core.py run_PEAKS(). PEAKS reports
spectrum *indices*, not raw scan numbers — so one or more index2scan.pkl
pickles from CONVERT_MZML are required to resolve true ScanNr values.

Invoked as a Nextflow template from modules/local/convert_peaks/main.nf.
PEAKS_CSV / INDEX2SCAN_PKLS / MIN_MATCH / PROCESS_NAME are interpolated by Nextflow.
"""
from __future__ import annotations

import pickle
import shlex
import sys
from pathlib import Path

import pandas as pd

PEAKS_CSV = "${peaks_csv}"
INDEX2SCAN_PKLS_STR = "${index2scan_pkls.join(' ')}"
MIN_MATCH = float("${min_match}")
PROCESS_NAME = "${task.process}"


def load_index2scan(paths: list[Path]) -> dict:
    merged: dict = {}
    for p in paths:
        with open(p, "rb") as fh:
            merged.update(pickle.load(fh))
    return merged


def convert(csv_path: Path, index2scan: dict, out_path: Path,
            min_match_fraction: float = 0.98) -> None:
    df = pd.read_csv(csv_path)
    df["Run"] = df["Source File"].apply(lambda x: Path(x).stem)
    df["scan_id"] = df["Run"] + "_" + df["Scan"].astype(str)
    df["ScanNr"] = df["scan_id"].map(index2scan)

    missing = df["ScanNr"].isna().sum()
    resolved = df["ScanNr"].notna().sum()
    frac = (resolved / len(df)) if len(df) else 0.0
    if frac < min_match_fraction:
        raise SystemExit(
            f"Only {resolved}/{len(df)} ({frac:.1%}) PEAKS spectrum indices "
            f"resolved to scan numbers. Did you search the same mzMLs that "
            f"CONVERT_MZML parsed? Below {min_match_fraction:.0%} threshold."
        )
    if missing:
        print(f"[convert_peaks] warning: {missing}/{len(df)} PEAKS rows "
              f"dropped (no matching scan)", file=sys.stderr)
    df = df[df["ScanNr"].notna()].copy()

    df["SpecId"] = df["Run"] + "_" + df["ScanNr"].astype(str)
    df["decoy"] = df["decoy"].astype(str).str.lower().map({"true": True, "false": False}).fillna(False)
    df["Label"] = df["decoy"].map({True: -1, False: 1})
    df["Proteins"] = df["Accession"].str.replace("#DECOY#", "rev_", regex=False)

    cols = ["SpecId", "Label", "ScanNr", "-10LgP", "Tag Length", "Delta RT",
            "MS2 Correlation", "ppm", "Peptide", "Proteins"]
    if "Delta 1/k0" in df.columns:
        cols.insert(-2, "Delta 1/k0")
    df[cols].to_csv(out_path, sep="\\t", index=False)
    print(f"[convert_peaks] wrote {len(df)} rows -> {out_path}", file=sys.stderr)


index2scan_paths = [Path(p) for p in shlex.split(INDEX2SCAN_PKLS_STR)]
index2scan = load_index2scan(index2scan_paths)
print(f"[convert_peaks] loaded {len(index2scan)} index->scan entries", file=sys.stderr)
convert(Path(PEAKS_CSV), index2scan, Path("peaks.pin"), MIN_MATCH)

with open("versions.yml", "w") as f:
    f.write(f'"{PROCESS_NAME}":\\n')
    f.write(f"    python: {sys.version.split()[0]}\\n")
    f.write(f"    pandas: {pd.__version__}\\n")
