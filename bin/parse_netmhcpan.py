#!/usr/bin/env python3
"""Parse netMHCpan / netMHCIIpan stdout output into a per-peptide best-binder TSV.

Extracted from immunopeptidomics/core.py netMHCpan() (L1802) and
get_best_binder() (L1780). Ranks by percentile rank (ascending) and keeps
the best MHC allele per peptide.
"""
from __future__ import annotations

import argparse
import sys
from io import StringIO
from pathlib import Path

import pandas as pd

NETMHC_I_COLS = ["Pos", "mhc", "peptide", "core", "Of", "Gp", "Gl",
                 "Ip", "Il", "Icore", "Identity", "score_EL", "rank",
                 "arrow", "bind_level"]
NETMHC_II_COLS = ["Pos", "mhc", "peptide", "Of", "core", "Core_Rel",
                  "Inverted", "Identity", "score_EL", "rank", "Exp_Bind",
                  "arrow", "bind_level"]


def parse(output_file: Path, tool: str) -> pd.DataFrame:
    lines = output_file.read_text().splitlines()
    if tool == "netMHCpan":
        hits = [ln for ln in lines if "    PEPLIST" in ln]
        cols = NETMHC_I_COLS
    else:
        hits = [ln for ln in lines if "    Sequence  " in ln]
        cols = NETMHC_II_COLS
    if not hits:
        return pd.DataFrame(columns=cols)
    df = pd.read_csv(StringIO("\n".join(hits)), sep=r"\s+",
                     header=None, names=cols)
    return df


def best_binder(df: pd.DataFrame, tool: str) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=["peptide", "mhc", "core", "rank",
                                     "bind_level", "tool"])
    df = df.sort_values(by="rank", ascending=True)
    best = (
        df.groupby("peptide")
        .agg({"mhc": "first", "core": "first",
              "rank": "first", "bind_level": "first"})
        .reset_index()
        .fillna("NB")
    )
    version = "-4.3" if "II" in tool else "-4.1"
    best["tool"] = tool + version
    return best


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path,
                   help="netMHCpan or netMHCIIpan stdout capture")
    p.add_argument("--tool", required=True,
                   choices=["netMHCpan", "netMHCIIpan"])
    p.add_argument("--out-parsed", required=True, type=Path,
                   help="all parsed lines TSV")
    p.add_argument("--out-best", required=True, type=Path,
                   help="best-binder-per-peptide TSV")
    args = p.parse_args()

    parsed = parse(args.output, args.tool)
    parsed.to_csv(args.out_parsed, sep="\t", index=False)
    best = best_binder(parsed, args.tool)
    best.to_csv(args.out_best, sep="\t", index=False)
    print(f"[parse_netmhcpan] {args.tool}: {len(parsed)} rows → "
          f"{len(best)} unique peptides", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
