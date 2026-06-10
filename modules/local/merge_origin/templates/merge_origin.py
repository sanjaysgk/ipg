#!/usr/bin/env python3
'''Append origins(Ensembl) region calls to the integrated peptide table.

Parses the kescull origins *_origins_rna.csv (two header rows, then one row per
peptide-transcript match) and grafts the normal-transcriptome region/frame/biotype/
ENSG/coordinate fields onto integrated_peptides_origin.tsv, joined by peptide
sequence. A peptide spanning several transcript matches has its distinct non-empty
values joined with ';'. Canonical and non-matched peptides are left blank.

origins_rna.csv normal-group columns (origins.c): 0 Num, 1 Sequence, 2 Conventional,
3 Transcript, 4 Mutations, 5 Metadata, 6 Category, 7 Chromosome, 8 PeptideCoords,
9 Strand, 10 ReferenceId, 11 AltGenomeCoords, 12 EnsemblCoords, 13 Biotype, 14 Gene.
'''
from __future__ import annotations

import csv
import sys

import pandas as pd

PEPTIDES_TSV = "${peptides_tsv}"
ORIGINS_CSV = "${origins_csv}"
PROCESS_NAME = "${task.process}"

# origins_rna.csv column index -> output column name (normal transcriptome group).
COLMAP = {
    2: "ensembl_conventional",
    3: "ensembl_transcript",
    4: "ensembl_variants",
    5: "ensembl_metadata",
    6: "ensembl_category",
    7: "ensembl_chrom",
    8: "ensembl_pep_coords",
    9: "ensembl_strand",
    10: "ensembl_ref_id",
    11: "ensembl_alt_coords",
    12: "ensembl_coords",
    13: "ensembl_biotype",
    14: "ensembl_gene",
}
OUT_COLS = list(COLMAP.values())
MIN_FIELDS = max(COLMAP) + 1   # a usable normal-group row needs >= 15 fields


def joinu(vals) -> str:
    '''Distinct, order-preserving, non-empty join with ';'.'''
    return ";".join(dict.fromkeys(v for v in vals if v))


def parse_origins(path: str) -> dict:
    '''sequence -> {out_col: aggregated value} from the rna.csv normal group.'''
    acc: dict = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        rows = list(reader)
    # Drop the two header rows; tolerate a short/empty file.
    for row in rows[2:]:
        if len(row) < MIN_FIELDS:
            continue
        seq = row[1].strip()
        if not seq:
            continue
        bucket = acc.setdefault(seq, {c: [] for c in OUT_COLS})
        for idx, col in COLMAP.items():
            bucket[col].append(row[idx].strip())
    return {seq: {c: joinu(v) for c, v in cols.items()} for seq, cols in acc.items()}


df = pd.read_csv(PEPTIDES_TSV, sep="\\t")
for c in OUT_COLS:
    df[c] = ""

origins = parse_origins(ORIGINS_CSV)

n_annot = 0
if origins and "peptide" in df.columns:
    for i in df.index:
        rec = origins.get(str(df.at[i, "peptide"]))
        if not rec:
            continue
        n_annot += 1
        for c in OUT_COLS:
            df.at[i, c] = rec[c]

print(f"[merge_origin] peptides annotated with Ensembl region calls: "
      f"{n_annot}/{len(df)}", file=sys.stderr)

df.to_csv("integrated_peptides_origin_region.tsv", sep="\\t", index=False)

with open("versions.yml", "w") as f:
    f.write(f'"{PROCESS_NAME}":\\n')
    f.write(f"    python: {sys.version.split()[0]}\\n")
    f.write(f"    pandas: {pd.__version__}\\n")
