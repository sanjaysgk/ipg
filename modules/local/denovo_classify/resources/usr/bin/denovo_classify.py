#!/usr/bin/env python3
"""Classify de novo peptide predictions against a cryptic + canonical database.

De novo engines (InstaNovo/Casanovo) predict peptide sequences directly from spectra
- this is NOT a database search. Each predicted peptide is looked up (substring match)
in the cryptic ORF database and, optionally, the canonical proteome to assign an
origin class:
  - canonical : present in the canonical proteome (known)
  - cryptic   : present in the cryptic DB but not canonical (corroborates DB search)
  - novel     : present in neither (de novo's unique contribution)

Confidence is the raw de novo log-probability (Phase 1a; Winnow calibration/FDR is a
separate, later step). Predictions below --min-log-prob or shorter than --min-length
are dropped.

I and L are isobaric (identical mass, indistinguishable by MS), so peptides and the
database are normalised I->L before matching.

Usage:
    denovo_classify.py --predictions preds.csv --cryptic-fasta cryptic.fasta \\
        [--canonical-fasta canonical.fasta] --min-log-prob -0.5 --out classified.tsv
"""

import argparse
import csv
import re
import sys


def fasta_concat_il(path):
    """Concatenate all sequences into one I/L-normalised, upper-cased string.

    Sequences are joined with '|' so a peptide (which never contains '|') cannot
    match across a protein boundary.
    """
    seqs, cur = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line)
        if cur:
            seqs.append("".join(cur))
    return "|".join(seqs).upper().replace("I", "L")


def il(pep):
    return pep.upper().replace("I", "L")


def clean_peptide(seq):
    """Keep only amino-acid letters (drop modification annotations, charges, etc.)."""
    return re.sub(r"[^A-Za-z]", "", seq or "").upper()


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--predictions", required=True, help="De novo predictions CSV.")
    p.add_argument("--cryptic-fasta", required=True, help="Cryptic ORF database FASTA.")
    p.add_argument("--canonical-fasta", help="Canonical proteome FASTA (optional).")
    p.add_argument("--peptide-col", default="predictions", help="Predicted-peptide column.")
    p.add_argument("--score-col", default="log_probs", help="De novo log-probability column.")
    p.add_argument("--id-col", default="spectrum_id", help="Spectrum-id column.")
    p.add_argument("--engine", default="instanovo", help="De novo engine name (tag in output).")
    p.add_argument("--min-log-prob", type=float, default=float("-inf"),
                   help="Drop predictions with log-prob below this.")
    p.add_argument("--min-length", type=int, default=7, help="Minimum peptide length to keep.")
    p.add_argument("--out", required=True, help="Output classified TSV.")
    args = p.parse_args()

    cryptic = fasta_concat_il(args.cryptic_fasta)
    canonical = fasta_concat_il(args.canonical_fasta) if args.canonical_fasta else ""

    n_in = n_kept = 0
    counts = {"canonical": 0, "cryptic": 0, "novel": 0}
    with open(args.predictions, newline="") as fh, open(args.out, "w", newline="") as out:
        reader = csv.DictReader(fh)
        for col in (args.peptide_col, args.score_col):
            if col not in (reader.fieldnames or []):
                sys.exit(f"ERROR: column '{col}' not in predictions CSV "
                         f"(have: {reader.fieldnames})")
        w = csv.writer(out, delimiter="\t")
        w.writerow(["spectrum_id", "peptide", "length", "log_prob", "class", "engine"])
        for row in reader:
            n_in += 1
            pep = clean_peptide(row.get(args.peptide_col))
            if len(pep) < args.min_length:
                continue
            try:
                score = float(row.get(args.score_col) or "nan")
            except ValueError:
                score = float("nan")
            if not (score >= args.min_log_prob):  # also drops NaN
                continue
            pep_il = il(pep)
            if canonical and pep_il in canonical:
                cls = "canonical"
            elif pep_il in cryptic:
                cls = "cryptic"
            else:
                cls = "novel"
            counts[cls] += 1
            n_kept += 1
            w.writerow([row.get(args.id_col, ""), pep, len(pep), f"{score:.6g}", cls, args.engine])

    sys.stderr.write(
        f"[denovo_classify] engine={args.engine} in={n_in} kept={n_kept} "
        f"canonical={counts['canonical']} cryptic={counts['cryptic']} novel={counts['novel']}\n"
    )


if __name__ == "__main__":
    main()
