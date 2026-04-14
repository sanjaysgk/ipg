#!/usr/bin/env python3
"""Convert a PEAKS Studio db.psms.csv export into a Percolator input (PIN).

Extracted from immunopeptidomics/core.py run_PEAKS() (L555). PEAKS reports
spectrum *indices*, not raw scan numbers — so one or more index2scan.pkl
pickles from CONVERT_MZML are required to resolve true ScanNr values.

The emitted PIN is directly consumable by mokapot via the existing MOKAPOT
module (engine == 'peaks' triggers no extra PIN cleanup there).

Only validated against PEAKS 12 exports with decoys present.
"""
from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import pandas as pd


def load_index2scan(paths: list[Path]) -> dict:
    merged: dict = {}
    for p in paths:
        with open(p, "rb") as fh:
            merged.update(pickle.load(fh))
    return merged


def convert(csv_path: Path, index2scan: dict, out_path: Path,
            min_match_fraction: float = 0.98) -> None:
    df = pd.read_csv(csv_path)
    df["Run"] = df["Source File"].str.replace(".mzML", "", regex=False)
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
    df["Label"] = df["decoy"].map({True: -1, False: 1})
    df["Proteins"] = df["Accession"].str.replace("#DECOY#", "rev_", regex=False)

    cols = ["SpecId", "Label", "ScanNr", "-10LgP", "Tag Length", "Delta RT",
            "MS2 Correlation", "ppm", "Peptide", "Proteins"]
    if "Delta 1/k0" in df.columns:          # timsTOF ion-mobility feature
        cols.insert(-2, "Delta 1/k0")
    df[cols].to_csv(out_path, sep="\t", index=False)
    print(f"[convert_peaks] wrote {len(df)} rows → {out_path}", file=sys.stderr)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv", required=True, type=Path,
                   help="PEAKS db.psms.csv export")
    p.add_argument("--index2scan", required=True, nargs="+", type=Path,
                   help="index2scan.pkl file(s) from CONVERT_MZML")
    p.add_argument("--out", required=True, type=Path,
                   help="output PIN file")
    p.add_argument("--min-match", type=float, default=0.98,
                   help="minimum fraction of PEAKS rows that must resolve (default: 0.98)")
    args = p.parse_args()
    index2scan = load_index2scan(args.index2scan)
    print(f"[convert_peaks] loaded {len(index2scan)} index→scan entries",
          file=sys.stderr)
    convert(args.csv, index2scan, args.out, args.min_match)
    return 0


if __name__ == "__main__":
    sys.exit(main())
