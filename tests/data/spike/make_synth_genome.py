#!/usr/bin/env python3
"""Generate a tiny synthetic genome + RNA-seq reads for the db_construct positive control.

A ~3 kb contig carries TWO single-exon transcripts:
  * gene 1 (CRYPTICALLY) — its exon embeds, in a fixed reading frame, a stop-flanked block of
    stop-free codons spelling the planted cryptic peptide CRYPTICALLY.
  * gene 2 (sentinel) — plain DNA, downstream of gene 1.

Why two genes: triple_translate.c silently DROPS the last transcript in the transcriptome (it
translates but never writes the final FASTA record — fine at scale, fatal for a 1-transcript
fixture). Coordinate-sorted, gene 1 is emitted first (written) and the sentinel is last
(absorbs the drop). A real transcriptome has many transcripts, so this is also more faithful.

db_construct: STAR align → StringTie assemble → gffread extract → triple_translate (3-frame,
split at stops, >=8 aa) → squish. Frame (block_offset mod 3) of gene 1 emits exactly
`CRYPTICALLY`, which must appear in <sample>_cryptic.fasta.

Pure stdlib, reproducible (fixed seed). Emits genome.fa, genes.gtf, genes.bed, paired
reads_{R1,R2}.fastq.gz, four VCFs (germline carries AF + real sites for GATK
GetPileupSummaries; the others are header-only known-sites for BQSR), and samplesheet_rnaseq.csv.
bgzip + tabix the VCFs afterwards (see the spike README).

    python tests/data/spike/make_synth_genome.py --out-dir tests/data/spike/genome
"""
from __future__ import annotations

import argparse
import gzip
import os
import random

CRYPTIC = "CRYPTICALLY"
CONTIG = "chrSPIKE"

AA_CODON = {
    "C": "TGC", "R": "CGC", "Y": "TAC", "P": "CCC", "T": "ACC",
    "I": "ATC", "A": "GCC", "L": "CTC",
}
STOP = "TAA"


def cryptic_block() -> str:
    """TAA + codons(CRYPTICALLY) + TAA — a self-framed, stop-flanked block."""
    return STOP + "".join(AA_CODON[a] for a in CRYPTIC) + STOP


def random_dna(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def reads_for(rng, seq, ridx, rl, fl, cov, tag="SPIKE"):
    """Paired FR reads tiling one exon (sense on +). Returns (r1_list, r2_list, next_idx).
    `tag` namespaces read IDs so per-rep fastqs carry unique names once merged."""
    ln = len(seq)
    n = max(1, cov * ln // (2 * rl))
    span = max(0, ln - fl)
    r1, r2 = [], []
    for _ in range(n):
        off = rng.randint(0, span)
        frag = seq[off:off + fl] if ln >= fl else seq
        a = frag[:rl]
        b = revcomp(frag[-rl:])
        r1.append(f"@{tag}_{ridx}/1\n{a}\n+\n{'I' * len(a)}\n")
        r2.append(f"@{tag}_{ridx}/2\n{b}\n+\n{'I' * len(b)}\n")
        ridx += 1
    return r1, r2, ridx


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--contig-len", type=int, default=3000)
    ap.add_argument("--exon-len", type=int, default=450)
    ap.add_argument("--read-len", type=int, default=100)
    ap.add_argument("--frag-len", type=int, default=200)
    ap.add_argument("--coverage", type=int, default=120)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.out_dir, exist_ok=True)

    # ---- two exons on one contig --------------------------------------------------------
    block = cryptic_block()
    g1_len, g2_len = args.exon_len, 300
    g1s = 1200                                  # gene1 exon (CRYPTICALLY), 1-based start
    g2s = g1s + g1_len + 350                    # gene2 (sentinel), downstream
    pad = (g1_len - len(block)) // 2
    g1_seq = random_dna(rng, pad) + block + random_dna(rng, g1_len - pad - len(block))
    g2_seq = random_dna(rng, g2_len)

    contig = list(random_dna(rng, args.contig_len))
    contig[g1s - 1:g1s - 1 + g1_len] = list(g1_seq)
    contig[g2s - 1:g2s - 1 + g2_len] = list(g2_seq)
    contig = "".join(contig)

    genes = [
        ("SPIKE_G1", "SPIKE_T1", g1s, g1s + g1_len - 1, g1_seq),   # CRYPTICALLY (written)
        ("SPIKE_G2", "SPIKE_T2", g2s, g2s + g2_len - 1, g2_seq),   # sentinel (TT drops last)
    ]
    d = os.path.abspath(args.out_dir)

    # ---- genome.fa ----------------------------------------------------------------------
    with open(f"{d}/genome.fa", "w") as f:
        f.write(f">{CONTIG}\n")
        for i in range(0, len(contig), 70):
            f.write(contig[i:i + 70] + "\n")

    # ---- GTF + BED12 --------------------------------------------------------------------
    with open(f"{d}/genes.gtf", "w") as gtf, open(f"{d}/genes.bed", "w") as bed:
        for gid, tid, s, e, _ in genes:
            attr = f'gene_id "{gid}"; transcript_id "{tid}";'
            for feat in ("gene", "transcript", "exon"):
                a = attr if feat != "exon" else attr + ' exon_number "1";'
                gtf.write(f"{CONTIG}\tsynth\t{feat}\t{s}\t{e}\t.\t+\t.\t{a}\n")
            bed.write(f"{CONTIG}\t{s - 1}\t{e}\t{tid}\t0\t+\t{s - 1}\t{e}\t0\t1\t{e - s + 1},\t0\n")

    # ---- reads tiling both exons --------------------------------------------------------
    r1_recs, r2_recs, ridx = [], [], 0
    for _, _, _, _, seq in genes:
        a, b, ridx = reads_for(rng, seq, ridx, args.read_len, args.frag_len, args.coverage)
        r1_recs += a
        r2_recs += b
    with gzip.open(f"{d}/reads_R1.fastq.gz", "wt") as f:
        f.write("".join(r1_recs))
    with gzip.open(f"{d}/reads_R2.fastq.gz", "wt") as f:
        f.write("".join(r2_recs))

    # ---- VCFs ---------------------------------------------------------------------------
    base = ("##fileformat=VCFv4.2\n"
            f"##contig=<ID={CONTIG},length={args.contig_len}>\n")
    cols = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    for name in ("dbsnp", "known_indels", "mills"):           # BQSR known-sites: header only
        with open(f"{d}/{name}.vcf", "w") as f:
            f.write(base + cols)
    # germline resource: GATK GetPileupSummaries needs an AF INFO field + biallelic SNP sites.
    with open(f"{d}/germline_resource.vcf", "w") as f:
        f.write(base + '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n' + cols)
        for p in (g1s + 60, g1s + 180, g1s + 300):
            ref = contig[p - 1]
            alt = next(b for b in "ACGT" if b != ref)
            f.write(f"{CONTIG}\t{p}\t.\t{ref}\t{alt}\t.\t.\tAF=0.5\n")

    # ---- RNA-seq samplesheet (absolute fastq paths) -------------------------------------
    with open(f"{d}/samplesheet_rnaseq.csv", "w") as f:
        f.write("sample,fastq_1,fastq_2,strandedness\n")
        f.write(f"SPIKE,{d}/reads_R1.fastq.gz,{d}/reads_R2.fastq.gz,forward\n")

    # ---- multi-sample / multi-rep fixtures (DISTINCT reads per rep) ----------------------
    # Two conditions (A, B) x three reps. Each rep draws from its own RNG stream and
    # carries a unique read-name tag, so a merge adds real depth (not pure duplicates)
    # and read names stay unique once reps are combined at MarkDuplicates. Drives both
    # styles WITHOUT hard-coding any DB count:
    #   samplesheet_rnaseq_multirep.csv : same sample id per condition -> reps merge -> 2 DBs
    #   samplesheet_rnaseq_perrep.csv   : distinct id per rep          -> 6 DBs
    conds = ("A", "B")
    reps = (1, 2, 3)
    rep_files = {}
    for ci, cond in enumerate(conds, 1):
        for rep in reps:
            rrng = random.Random(args.seed + 1000 * ci + rep)
            r1r, r2r, ri = [], [], 0
            for _, _, _, _, seq in genes:
                a, b, ri = reads_for(rrng, seq, ri, args.read_len, args.frag_len,
                                     args.coverage, tag=f"SPIKE{cond}r{rep}")
                r1r += a
                r2r += b
            p1 = f"{d}/reads_{cond}_rep{rep}_R1.fastq.gz"
            p2 = f"{d}/reads_{cond}_rep{rep}_R2.fastq.gz"
            with gzip.open(p1, "wt") as f:
                f.write("".join(r1r))
            with gzip.open(p2, "wt") as f:
                f.write("".join(r2r))
            rep_files[(cond, rep)] = (p1, p2)

    with open(f"{d}/samplesheet_rnaseq_multirep.csv", "w") as f:
        f.write("sample,fastq_1,fastq_2,strandedness\n")
        for cond in conds:
            for rep in reps:
                p1, p2 = rep_files[(cond, rep)]
                f.write(f"SPIKE_{cond},{p1},{p2},forward\n")

    with open(f"{d}/samplesheet_rnaseq_perrep.csv", "w") as f:
        f.write("sample,fastq_1,fastq_2,strandedness\n")
        for cond in conds:
            for rep in reps:
                p1, p2 = rep_files[(cond, rep)]
                f.write(f"SPIKE_{cond}_rep{rep},{p1},{p2},forward\n")

    print(f"[make_synth_genome] contig {CONTIG} {args.contig_len} bp")
    print(f"[make_synth_genome] gene1 CRYPTICALLY exon {g1s}-{g1s + g1_len - 1} "
          f"(block frame {pad % 3}); gene2 sentinel exon {g2s}-{g2s + g2_len - 1}")
    print(f"[make_synth_genome] {len(r1_recs)} read pairs; germline VCF has 3 AF sites")
    print(f"[make_synth_genome] wrote fixtures to {d}; planted peptide: {CRYPTIC}")
    print("[make_synth_genome] next: bgzip + tabix the 4 .vcf files (see README)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
