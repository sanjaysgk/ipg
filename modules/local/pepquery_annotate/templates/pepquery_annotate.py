#!/usr/bin/env python3
'''Annotate the integrated peptide table with PepQuery2 validation status.

A cryptic peptide is 'confident' if it has at least one PSM in psm_rank.txt
with confident == 'Yes'. Canonical peptides are not tested (PepQuery validates
the cryptic subset only), so they carry an empty status.
'''
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

PEPTIDES_TSV = "${peptides_tsv}"
PSM_RANK = "${psm_rank}"
PROCESS_NAME = "${task.process}"

pep = pd.read_csv(PEPTIDES_TSV, sep="\\t")

# Build per-peptide verdict from PepQuery. Missing/empty result => no confident
# hits; every tested cryptic peptide then resolves to 'not_confident'.
confident_peps: set = set()
best_pvalue: dict = {}
pq_path = Path(PSM_RANK)
if pq_path.exists() and pq_path.stat().st_size > 0:
    pq = pd.read_csv(pq_path, sep="\\t")
    if not pq.empty and "peptide" in pq.columns and "confident" in pq.columns:
        yes = pq[pq["confident"].astype(str).str.strip().str.lower() == "yes"]
        confident_peps = set(yes["peptide"].astype(str))
        if "pvalue" in pq.columns:
            real = pq[pq["pvalue"] < 100]   # 100 is PepQuery's "not scored" sentinel
            best_pvalue = (real.groupby("peptide")["pvalue"].min().to_dict())

def status(row) -> str:
    if str(row.get("class", "")) != "cryptic":
        return ""
    return "confident" if row["peptide"] in confident_peps else "not_confident"

pep["pepquery_status"] = pep.apply(status, axis=1)
pep["pepquery_pvalue"] = pep["peptide"].map(best_pvalue)

n_conf = int((pep["pepquery_status"] == "confident").sum())
n_cryp = int((pep.get("class", pd.Series(dtype=str)) == "cryptic").sum())
print(f"[pepquery_annotate] cryptic peptides confident: {n_conf}/{n_cryp}",
      file=sys.stderr)

pep.to_csv("integrated_peptides_validated.tsv", sep="\\t", index=False)

with open("versions.yml", "w") as f:
    f.write(f'"{PROCESS_NAME}":\\n')
    f.write(f"    python: {sys.version.split()[0]}\\n")
    f.write(f"    pandas: {pd.__version__}\\n")
