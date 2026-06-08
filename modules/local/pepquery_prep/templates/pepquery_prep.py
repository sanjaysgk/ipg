#!/usr/bin/env python3
'''Prepare PepQuery2 inputs for cryptic-peptide validation.

From the integrated peptide table, pull the cryptic-class peptide sequences
into a one-per-line list, and build a canonical-only reference FASTA (the
proteome PepQuery competes each cryptic candidate against). Decoys and
non-canonical (cryptic) entries are excluded from the reference.
'''
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def is_canonical(header_id: str, prefixes: list[str]) -> bool:
    hid = header_id[4:] if header_id.startswith("rev_") else header_id
    return any(hid.startswith(pre) for pre in prefixes)


PEPTIDES_TSV = "${peptides_tsv}"
SEARCH_FASTA = "${search_fasta}"
CANONICAL_PREFIXES = [p for p in "${canonical_prefixes}".split(",") if p]
PROCESS_NAME = "${task.process}"

# 1. Cryptic peptide list (unique sequences, class == cryptic).
df = pd.read_csv(PEPTIDES_TSV, sep="\\t")
if "class" in df.columns:
    cryptic = df[df["class"] == "cryptic"]["peptide"].dropna().unique()
else:
    cryptic = []
with open("cryptic_peptides.txt", "w") as fh:
    for pep in cryptic:
        fh.write(f"{pep}\\n")
print(f"[pepquery_prep] cryptic peptides: {len(cryptic)}", file=sys.stderr)

# 2. Canonical-only reference FASTA.
n_ref = 0
with open("canonical_ref.fasta", "w") as out:
    for rec in SeqIO.parse(SEARCH_FASTA, "fasta"):
        if is_canonical(rec.id, CANONICAL_PREFIXES):
            SeqIO.write(rec, out, "fasta")
            n_ref += 1
print(f"[pepquery_prep] canonical reference proteins: {n_ref}", file=sys.stderr)

with open("versions.yml", "w") as f:
    f.write(f'"{PROCESS_NAME}":\\n')
    f.write(f"    python: {sys.version.split()[0]}\\n")
    f.write(f"    pandas: {pd.__version__}\\n")
    import Bio
    f.write(f"    biopython: {Bio.__version__}\\n")
