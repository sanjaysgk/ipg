#!/usr/bin/env python3
"""Pick the winning GibbsCluster-2.0 cluster by max KLD sum and dump peptide→cluster.

Extracted from immunopeptidomics/core.py run_Gibbs() (L2021 — post-run parsing).
"""
from __future__ import annotations

import argparse
import sys
from io import StringIO
from pathlib import Path

import pandas as pd

CLUSTER_COLS = [
    "G", "cluster", "nr0", "peptide", "gibbs_core",
    "o", "nr1", "ip", "nr2", "il", "nr3", "dp", "nr4", "dl",
    "nrx", "peplist", "sS", "nr5", "bgG", "nr6", "bgS", "nr7", "cS", "nr8",
]


def pick_cluster(kld_path: Path) -> int:
    kld = pd.read_csv(kld_path, sep="\t")
    sums = kld.sum(axis=1)
    return int(sums.idxmax())


def parse_cluster_output(out_path: Path) -> pd.DataFrame:
    lines = out_path.read_text().splitlines()
    hits = [ln for ln in lines if "    Peplist" in ln]
    if not hits:
        return pd.DataFrame(columns=["peptide", "cluster", "gibbs_core"])
    df = pd.read_csv(StringIO("\n".join(hits)), sep=r"\s+",
                     header=None, names=CLUSTER_COLS)
    return df[["peptide", "cluster", "gibbs_core"]]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--kld", required=True, type=Path,
                   help="gibbs.KLDvsClusters.tab from GibbsCluster output")
    p.add_argument("--gibbs-dir", required=True, type=Path,
                   help="GibbsCluster output directory with res/ subfolder")
    p.add_argument("--out", required=True, type=Path,
                   help="per-peptide cluster assignment TSV")
    args = p.parse_args()

    winner = pick_cluster(args.kld)
    res_file = args.gibbs_dir / "res" / f"gibbs.{winner}g.out"
    if not res_file.exists():
        raise SystemExit(f"cluster result file missing: {res_file}")
    df = parse_cluster_output(res_file)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"[parse_gibbs] winning cluster count: {winner} ({len(df)} peptides)",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
