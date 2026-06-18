#!/usr/bin/env python3
"""Protein/ORF-group FDR via PickedGroupFDR on MS2Rescore/mokapot output.

Runs the picked target-decoy protein-group FDR using the kusterlab example
workflow for MS2Rescore input (`picked_protein_group_no_remap`). picked's
ms2rescore reader keys on `spectrum_id` and reads `protein_list` / `mokapot
score` / `mokapot q-value` / `mokapot pep` directly, and recognises `rev_`
decoys natively, so no column renaming is needed.

Scope is set by `fdr_class`:
  global  : full target+decoy list -> one overall protein-group FDR (default).
  cryptic : restrict to the cryptic class first (peptides whose mapped proteins
            are all non-canonical; decoys inherit class via rev_-stripping, the
            same rule as INTEGRATE_ENGINES) -> class-specific cryptic FDR.

Invoked as a Nextflow `template` from modules/local/picked_group_fdr/main.nf.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

MOKAPOT_PSMS  = "${mokapot_psms}"
MOKAPOT_DECOY = "${mokapot_decoy}"
DB_FASTA      = "${cryptic_fasta}"
PREFIX        = "${prefix}"
FDR           = float("${fdr}")
FDR_CLASS     = "${fdr_class}"
CANON_PREFIXES = [p for p in "${canonical_prefixes}".split(",") if p]
PROCESS_NAME  = "${task.process}"

PG_OUT = f"{PREFIX}_protein_groups.tsv"


def is_cryptic(protein_field: str) -> bool:
    """Cryptic iff every mapped protein is non-canonical. Decoys inherit their
    target's class: strip a leading rev_ before testing the prefix."""
    s = protein_field.strip().strip("[]").replace("'", "").replace('"', "")
    s = s.replace(", ", ";").replace(",", ";").replace(" ", ";")
    proteins = [p[4:] if p.startswith("rev_") else p for p in s.split(";") if p]
    if not proteins:
        return False
    return not any(p.startswith(pre) for p in proteins for pre in CANON_PREFIXES)


def filter_cryptic(in_path: str, out_path: str) -> int:
    """Keep header + rows whose protein_list classifies cryptic. Lines are kept
    byte-for-byte (no reformatting) so picked's reader sees its native format."""
    kept = 0
    with open(in_path) as fh, open(out_path, "w") as out:
        header = fh.readline()
        out.write(header)
        cols = [c.strip().lower() for c in header.rstrip("\\n").split("\\t")]
        try:
            pcol = cols.index("protein_list")
        except ValueError:
            raise SystemExit(
                f"[picked_group_fdr] '{in_path}' has no 'protein_list' column "
                f"(headers: {cols}); expected MS2Rescore/mokapot output."
            )
        for line in fh:
            fields = line.rstrip("\\n").split("\\t")
            if len(fields) > pcol and is_cryptic(fields[pcol]):
                out.write(line)
                kept += 1
    return kept


too_few = False
if FDR_CLASS == "cryptic":
    n_target = filter_cryptic(MOKAPOT_PSMS, "fdr.psms.txt")
    n_decoy  = filter_cryptic(MOKAPOT_DECOY, "fdr.decoy.psms.txt")
    print(f"[picked_group_fdr] cryptic PSMs: {n_target} target, {n_decoy} decoy",
          file=sys.stderr)
    psms_in, decoy_in = "fdr.psms.txt", "fdr.decoy.psms.txt"
    too_few = n_target == 0 or n_decoy == 0
else:
    # global: full target+decoy, i.e. the picked MS2Rescore example as-is.
    psms_in, decoy_in = MOKAPOT_PSMS, MOKAPOT_DECOY

if too_few:
    # No cryptic competition possible (e.g. canonical-only DB). Emit an empty
    # result + header so downstream joins don't break; warn loudly.
    print(f"[picked_group_fdr] WARNING: cryptic class has {n_target} target / "
          f"{n_decoy} decoy PSMs — too few for a protein-group FDR; writing empty "
          f"{PG_OUT}.", file=sys.stderr)
    Path(PG_OUT).write_text("Protein IDs\\tProtein group FDR\\n")
else:
    cmd = [
        sys.executable, "-m", "picked_group_fdr",
        "--fasta", DB_FASTA,
        "--perc_evidence", psms_in, decoy_in,
        "--protein_groups_out", PG_OUT,
        "--output_format", "fragpipe",
        "--methods", "picked_protein_group_no_remap",
        "--protein_group_fdr_threshold", str(FDR),
    ]
    print("[picked_group_fdr] " + " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)

try:
    import importlib.metadata as im
    pgf_ver = im.version("picked_group_fdr")
except Exception:
    pgf_ver = "unknown"

with open("versions.yml", "w") as f:
    f.write(f'"{PROCESS_NAME}":\\n')
    f.write(f"    python: {sys.version.split()[0]}\\n")
    f.write(f"    picked_group_fdr: {pgf_ver}\\n")
