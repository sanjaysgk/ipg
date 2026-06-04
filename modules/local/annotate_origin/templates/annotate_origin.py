#!/usr/bin/env python3
'''Backtrack cryptic peptides to their RNA-seq genomic origin — no network.

For each cryptic-class peptide: find the exact source ORF by substring-matching
against the unsquished 3-frame translation (resolving the squish concatenation),
decode the ORF id into transcript/frame/orf, then enrich from the StringTie/
gffcompare artifacts:
  - combined.gtf  -> gene symbol, genomic locus (chrom:start-end:strand),
                     reference transcript, novelty class code, TSS id
  - .tmap         -> transcript expression (FPKM, TPM) + length
Canonical peptides are left blank. Reuses artifacts DB_CONSTRUCT already emits;
no Ensembl/REST (use the optional ORIGINS module for 5'/3'UTR region calls).
'''
from __future__ import annotations

import re
import sys

import pandas as pd
from Bio import SeqIO

ORF_RE = re.compile(r"(TCONS_\\d+)_f(\\d+)p(\\d+)")
_ATTR = {
    "gene_id": re.compile(r'gene_id "([^"]+)"'),
    "gene_name": re.compile(r'gene_name "([^"]+)"'),
    "oid": re.compile(r'oId "([^"]+)"'),
    "nearest_ref": re.compile(r'cmp_ref "([^"]+)"'),
    "novelty_class": re.compile(r'class_code "([^"]+)"'),
    "tss_id": re.compile(r'tss_id "([^"]+)"'),
}
_TID = re.compile(r'transcript_id "([^"]+)"')

PEPTIDES_TSV = "${peptides_tsv}"
ORF_FASTA = "${orf_fasta}"
COMBINED_GTF = "${combined_gtf}"
TMAP = "${tmap}"
PROCESS_NAME = "${task.process}"


def parse_gtf(path: str) -> dict:
    '''transcript_id -> per-transcript provenance from the combined GTF.'''
    out: dict = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\\n").split("\\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            m = _TID.search(cols[8])
            if not m:
                continue
            rec = {
                "chrom": cols[0], "start": cols[3], "end": cols[4],
                "strand": cols[6],
            }
            for key, rx in _ATTR.items():
                mm = rx.search(cols[8])
                rec[key] = mm.group(1) if mm else ""
            out[m.group(1)] = rec
    return out


def parse_tmap(path: str) -> dict:
    '''gffcompare qry_id (StringTie id) -> {FPKM, TPM, len}.'''
    out: dict = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\\n").split("\\t")
        idx = {c: i for i, c in enumerate(header)}
        need = ("qry_id", "FPKM", "TPM", "len")
        if not all(k in idx for k in need):
            return out
        for line in fh:
            p = line.rstrip("\\n").split("\\t")
            if len(p) <= idx["qry_id"]:
                continue
            out[p[idx["qry_id"]]] = {
                "FPKM": p[idx["FPKM"]], "TPM": p[idx["TPM"]], "len": p[idx["len"]],
            }
    return out


orfs = [(rec.id, str(rec.seq)) for rec in SeqIO.parse(ORF_FASTA, "fasta")]
gtf = parse_gtf(COMBINED_GTF)
tmap = parse_tmap(TMAP)

df = pd.read_csv(PEPTIDES_TSV, sep="\\t")
is_cryptic = df.get("class", pd.Series([""] * len(df))) == "cryptic"

cols = ["source_orf", "source_transcript", "frame", "orf", "gene_name",
        "gene_locus", "genomic_locus", "strand", "nearest_ref",
        "novelty_class", "tss_id", "FPKM", "TPM"]
for c in cols:
    df[c] = ""


def joinu(vals) -> str:
    return ";".join(dict.fromkeys(v for v in vals if v))


n_located = 0
for i in df.index[is_cryptic]:
    pep = str(df.at[i, "peptide"])
    hits = [oid for oid, seq in orfs if pep in seq]
    if not hits:
        continue
    n_located += 1
    df.at[i, "source_orf"] = ";".join(hits)
    acc = {c: [] for c in cols if c not in ("source_orf",)}
    for oid in hits:
        m = ORF_RE.search(oid)
        if not m:
            continue
        tcons, frame, orf_n = m.group(1), m.group(2), m.group(3)
        g = gtf.get(tcons, {})
        t = tmap.get(g.get("oid", ""), {})
        acc["source_transcript"].append(tcons)
        acc["frame"].append(frame)
        acc["orf"].append(orf_n)
        acc["gene_name"].append(g.get("gene_name", ""))
        acc["gene_locus"].append(g.get("gene_id", ""))
        if g:
            acc["genomic_locus"].append(f"{g['chrom']}:{g['start']}-{g['end']}")
        acc["strand"].append(g.get("strand", ""))
        acc["nearest_ref"].append(g.get("nearest_ref", ""))
        acc["novelty_class"].append(g.get("novelty_class", ""))
        acc["tss_id"].append(g.get("tss_id", ""))
        acc["FPKM"].append(t.get("FPKM", ""))
        acc["TPM"].append(t.get("TPM", ""))
    for c, vals in acc.items():
        df.at[i, c] = joinu(vals)

print(f"[annotate_origin] cryptic peptides located in ORF set: "
      f"{n_located}/{int(is_cryptic.sum())}", file=sys.stderr)

df.to_csv("integrated_peptides_origin.tsv", sep="\\t", index=False)

with open("versions.yml", "w") as f:
    f.write(f'"{PROCESS_NAME}":\\n')
    f.write(f"    python: {sys.version.split()[0]}\\n")
    f.write(f"    pandas: {pd.__version__}\\n")
    import Bio
    f.write(f"    biopython: {Bio.__version__}\\n")
