#!/usr/bin/env python3
"""Prepare a target-decoy FASTA database for MS search.

Reads a FASTA file and checks for decoy sequences (prefixed with 'rev_').
If none are found, generates reversed decoy sequences and writes a
concatenated target-decoy FASTA. If decoys already exist, copies the
input unchanged.

Invoked as a Nextflow `template` from modules/local/prepare_fasta/main.nf.
INPUT_FASTA / OUTPUT_FASTA / CONTAMINANT are interpolated by Nextflow.
"""

import os
import platform
import shutil
import sys

INPUT_FASTA = "${fasta}"
OUTPUT_FASTA = "${prefix}_tda.fasta"
CONTAMINANT = "CON__"
PROCESS_NAME = "${task.process}"


def prepare_fasta(input_fasta, output_fasta, contaminant_prefix):
    target_count = 0
    decoy_count = 0
    target_seq = {}
    headers = {}

    with open(input_fasta, "r") as fh:
        prot_id = None
        for line in fh:
            line = line.rstrip()
            if line.startswith(">rev_"):
                decoy_count += 1
                prot_id = None
            elif line.startswith(">"):
                target_count += 1
                prot_id = line.split()[0][1:]
                headers[prot_id] = line[1:]
                target_seq[prot_id] = ""
            elif prot_id:
                target_seq[prot_id] += line

    print(
        f"Input FASTA contains {target_count} target and "
        f"{decoy_count} decoy sequences.",
        file=sys.stderr,
    )

    if decoy_count > 0:
        print("Decoys already present, copying input.", file=sys.stderr)
        shutil.copy2(input_fasta, output_fasta)
    else:
        print("No decoys found, generating reverse target-decoy database.", file=sys.stderr)
        with open(output_fasta, "w") as out:
            for prot_id, seq in target_seq.items():
                rev_seq = seq[::-1]
                out.write(f">{headers[prot_id]}\\n{seq}\\n")
                out.write(f">rev_{prot_id} decoy_{prot_id}\\n{rev_seq}\\n")
        print(
            f"Wrote {len(target_seq)} targets + {len(target_seq)} decoys "
            f"to {os.path.basename(output_fasta)}",
            file=sys.stderr,
        )


prepare_fasta(INPUT_FASTA, OUTPUT_FASTA, CONTAMINANT)

with open("versions.yml", "w") as f:
    f.write(f'"{PROCESS_NAME}":\\n')
    f.write(f"    python: {platform.python_version()}\\n")
    f.write('    prepare_fasta: "sanjaysgk/ipg"\\n')

