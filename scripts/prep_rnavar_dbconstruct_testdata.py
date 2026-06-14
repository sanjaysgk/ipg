#!/usr/bin/env python3
"""Prepare the rnavar chr22 db_construct e2e test data (provided-ref topology).

The nf-core/rnavar chr22 bundle ships fasta/fai/dict/gtf, a packed STAR index,
and dbsnp/mills VCFs — i.e. references are PROVIDED (the real-run channel
topology), unlike the synthetic spike which builds them. This script fills the
three gaps and fans the single RNA sample into a multi-sample / multi-rep matrix:

  1. germline_af_only.vcf.gz  — a handful of biallelic chr22 SNP sites with an
        AF INFO field (Mutect2 GetPileupSummaries / CalculateContamination need
        a gnomAD-style germline resource; the bundle has none).
  2. genes.bed                — BED12 gene models from genome.gtf (RSeQC
        infer_experiment requires a BED; prepare_genome does not derive it).
  3. star/                    — unpacked from star.tar.gz (provide the dir).
  4. reads/rnavar_rep{1..6}   — the 4159 real read pairs split round-robin into
        six distinct rep files (genuine per-rep data, a real merge adds depth).
  5. two samplesheets (sample,fastq_1,fastq_2,strandedness,rep,condition):
        merged  : 2 conditions x 3 reps, shared id/condition -> 2 cryptic DBs
        per-rep : 6 distinct ids                              -> 6 cryptic DBs

known_indels reuses the bundle's mills VCF (Mills IS a known-indels set). All
outputs land under <out>; nothing here is committed — the config points at them.

    python scripts/prep_rnavar_dbconstruct_testdata.py \\
        --bundle /fs04/.../test-datasets/rnavar \\
        --out    /fs04/.../test-datasets/rnavar/dbconstruct
"""
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
from collections import defaultdict

N_REPS = 6
CONDS = ("RVcondA", "RVcondB")   # 3 reps each in the merged sheet


def read_fasta_one(path: str) -> tuple[str, str]:
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
            else:
                seq.append(line.strip())
    return name, "".join(seq)


def exons_by_transcript(gtf: str):
    """transcript_id -> (chrom, strand, [(start,end), ...])  (1-based, inclusive)."""
    tx = defaultdict(lambda: [None, None, []])
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            attrs = f[8]
            tid = attrs.split('transcript_id "', 1)[1].split('"', 1)[0]
            rec = tx[tid]
            rec[0], rec[1] = f[0], f[6]
            rec[2].append((int(f[3]), int(f[4])))
    return tx


def write_germline_vcf(contig: str, seq: str, exons, out_vcf: str) -> int:
    """A few biallelic SNP sites inside expressed exons, with AF INFO."""
    # take the first exon of the first few transcripts; offset into the exon so
    # the site is comfortably inside a covered region.
    sites = []
    for _tid, (_chrom, _strand, exs) in list(exons.items()):
        if not exs:
            continue
        s, e = sorted(exs)[0]
        pos = min(e, s + 40)
        if 1 <= pos <= len(seq):
            sites.append(pos)
        if len(sites) >= 6:
            break
    sites = sorted(set(sites))
    with open(out_vcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID={contig},length={len(seq)}>\n")
        f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for p in sites:
            ref = seq[p - 1].upper()
            if ref not in "ACGT":
                continue
            alt = next(b for b in "ACGT" if b != ref)
            f.write(f"{contig}\t{p}\t.\t{ref}\t{alt}\t.\t.\tAF=0.5\n")
    return len(sites)


def write_bed12(exons, out_bed: str) -> int:
    n = 0
    with open(out_bed, "w") as f:
        for tid, (chrom, strand, exs) in exons.items():
            if not chrom or not exs:
                continue
            exs = sorted(exs)
            tx_start = exs[0][0] - 1            # BED 0-based
            tx_end = exs[-1][1]
            sizes = [str(e - s + 1) for s, e in exs]
            starts = [str((s - 1) - tx_start) for s, _ in exs]
            f.write("\t".join([
                chrom, str(tx_start), str(tx_end), tid, "0", strand,
                str(tx_start), str(tx_end), "0",
                str(len(exs)), ",".join(sizes) + ",", ",".join(starts) + ",",
            ]) + "\n")
            n += 1
    return n


def split_fastq(r1: str, r2: str, out_dir: str, n_reps: int):
    def read_recs(path):
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rt") as fh:
            buf = []
            for line in fh:
                buf.append(line)
                if len(buf) == 4:
                    yield "".join(buf)
                    buf = []
    os.makedirs(out_dir, exist_ok=True)
    outs1 = [gzip.open(f"{out_dir}/rnavar_rep{i+1}_R1.fastq.gz", "wt") for i in range(n_reps)]
    outs2 = [gzip.open(f"{out_dir}/rnavar_rep{i+1}_R2.fastq.gz", "wt") for i in range(n_reps)]
    n = 0
    for j, (a, b) in enumerate(zip(read_recs(r1), read_recs(r2))):
        k = j % n_reps
        outs1[k].write(a)
        outs2[k].write(b)
        n += 1
    for fh in outs1 + outs2:
        fh.close()
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bundle", required=True, help="rnavar bundle dir")
    ap.add_argument("--out", required=True, help="staging output dir (scratch2)")
    args = ap.parse_args()

    b = os.path.abspath(args.bundle)
    d = os.path.abspath(args.out)
    g = f"{b}/genome"
    os.makedirs(d, exist_ok=True)
    os.makedirs(f"{d}/reads", exist_ok=True)

    # 3. unpack STAR index
    subprocess.run(["tar", "-xzf", f"{g}/star/star.tar.gz", "-C", d], check=True)
    star_dir = f"{d}/star"
    if not os.path.isdir(star_dir):
        # tar may extract into a named subdir; find it
        cand = [p for p in os.listdir(d) if os.path.isdir(f"{d}/{p}") and "star" in p.lower()]
        star_dir = f"{d}/{cand[0]}" if cand else d

    # genome + exon model
    contig, seq = read_fasta_one(f"{g}/genome.fasta")
    exons = exons_by_transcript(f"{g}/genome.gtf")

    # 1. germline af-only VCF -> bgzip + tabix
    n_sites = write_germline_vcf(contig, seq, exons, f"{d}/germline_af_only.vcf")
    subprocess.run(["bgzip", "-f", f"{d}/germline_af_only.vcf"], check=True)
    subprocess.run(["tabix", "-p", "vcf", f"{d}/germline_af_only.vcf.gz"], check=True)

    # known_indels: a DISTINCT copy of mills (Mills IS a known-indels set),
    # renamed so the BQSR known_sites list [dbsnp, known_indels, mills] has three
    # distinct filenames — reusing mills verbatim causes a staging name-collision.
    shutil.copy(f"{g}/vcf/mills_and_1000G.indels.vcf.gz", f"{d}/known_indels.vcf.gz")
    shutil.copy(f"{g}/vcf/mills_and_1000G.indels.vcf.gz.tbi", f"{d}/known_indels.vcf.gz.tbi")

    # 2. BED12 gene models
    n_tx = write_bed12(exons, f"{d}/genes.bed")

    # 4. split RNA reads into 6 distinct reps
    n_pairs = split_fastq(f"{b}/fastq/test_rnaseq_1.fastq.gz",
                          f"{b}/fastq/test_rnaseq_2.fastq.gz", f"{d}/reads", N_REPS)

    # 5. samplesheets
    def row(sample, rep, cond, i):
        r1 = f"{d}/reads/rnavar_rep{i}_R1.fastq.gz"
        r2 = f"{d}/reads/rnavar_rep{i}_R2.fastq.gz"
        return f"{sample},{r1},{r2},forward,{rep},{cond}\n"

    with open(f"{d}/samplesheet_rnavar_merged.csv", "w") as f:
        f.write("sample,fastq_1,fastq_2,strandedness,rep,condition\n")
        i = 0
        for cond in CONDS:
            for rep in (1, 2, 3):
                i += 1
                f.write(row(cond, rep, cond, i))

    with open(f"{d}/samplesheet_rnavar_perrep.csv", "w") as f:
        f.write("sample,fastq_1,fastq_2,strandedness,rep,condition\n")
        i = 0
        for ci, cond in enumerate(CONDS):
            for rep in (1, 2, 3):
                i += 1
                f.write(row(f"RV_rep{i}", 1, cond, i))

    print(f"[prep_rnavar] contig {contig} {len(seq)} bp; {n_tx} transcripts -> genes.bed")
    print(f"[prep_rnavar] germline_af_only.vcf.gz: {n_sites} AF sites; STAR index: {star_dir}")
    print(f"[prep_rnavar] split {n_pairs} read pairs into {N_REPS} reps")
    print(f"[prep_rnavar] wrote merged (2 DBs) + per-rep (6 DBs) samplesheets to {d}")
    print(f"[prep_rnavar] known_indels: reuse {g}/vcf/mills_and_1000G.indels.vcf.gz")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
