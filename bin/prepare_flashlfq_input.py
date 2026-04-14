#!/usr/bin/env python3
"""Build a FlashLFQ input TSV from the integrated peptides table.

Extracted from immunopeptidomics/core.py run_FlashLFQ() (L1119). The
upstream implementation reads the post-rescore mokapot peptide files
directly; after INTEGRATE_ENGINES we already have a per-peptide table
with run, charge, retention_time and protein_list — this script reshapes
it into FlashLFQ's expected column set.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

FLASHLFQ_COLS = [
    "File Name", "Scan Retention Time", "Precursor Charge",
    "Base Sequence", "Full Sequence", "Peptide Monoisotopic Mass",
    "Protein Accession",
]


def explode_runs(df: pd.DataFrame) -> pd.DataFrame:
    """Turn PSMs_run_<name> columns into a long form keyed by run."""
    run_cols = [c for c in df.columns if c.startswith("PSMs_run_")]
    if not run_cols:
        raise SystemExit("no PSMs_run_* columns in integrated peptides input")
    long = df.melt(
        id_vars=[c for c in df.columns if c not in run_cols],
        value_vars=run_cols,
        var_name="run_col",
        value_name="psm_count",
    )
    long = long[long["psm_count"].fillna(0) > 0].copy()
    long["File Name"] = long["run_col"].str.replace("PSMs_run_", "", regex=False)
    return long


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--peptides", required=True, type=Path,
                   help="integrated_peptides.tsv from INTEGRATE_ENGINES")
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()

    df = pd.read_csv(args.peptides, sep="\t")
    long = explode_runs(df)

    long["Full Sequence"] = long["peptidoform"].fillna(long["peptide"])
    long["Base Sequence"] = long["Full Sequence"].apply(
        lambda x: re.sub(r"\[.*?\]", "", str(x)).replace("-", "")
    )
    long["Peptide Monoisotopic Mass"] = long.get("calcmass", "")
    long["Protein Accession"] = long["protein_list"].fillna("")
    long["Precursor Charge"] = long.get("charge", 2)
    long["Scan Retention Time"] = long.get("retention_time", 0.0)

    out = long[FLASHLFQ_COLS].drop_duplicates(
        subset=["Full Sequence", "Precursor Charge", "File Name"]
    )
    out.to_csv(args.out, sep="\t", index=False)
    print(f"[flashlfq_input] wrote {len(out)} rows → {args.out}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
