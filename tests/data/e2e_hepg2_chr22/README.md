# E2E HepG2 chr22 test fixture

End-to-end test data for the full immunopeptidogenomics pipeline
(db_construct → ms_search) using the same HepG2 cell line for both
RNA-seq and MS/MS immunopeptidome data.

## Source

| Component | Source | Accession |
|-----------|--------|-----------|
| RNA-seq | ENCODE polyA+ RNA-seq, Graveley lab (UConn) | ENCSR985KAT rep1 |
| MS data | nf-core/mhcquant HepG2 immunopeptidome | PXD009752 |
| Reference | GENCODE v44, GRCh38 primary assembly, chr22 only | — |

## Contents (~65 MB)

```
e2e_hepg2_chr22/
├── samplesheet_rnaseq.csv    # db_construct input (chr22 paired reads)
├── samplesheet_ms.csv        # ms_search input (HepG2 clean mzML)
├── fastq/                    # 42 MB — chr22-only paired reads
│   ├── HepG2_chr22_R1.fastq.gz
│   └── HepG2_chr22_R2.fastq.gz
└── README.md
```

Reference + variant calling + MS data files are in sibling directories
(`ipg-db-constructions/` and `ms_search/`).

## Build script

```bash
bash bin/build_e2e_hepg2_chr22.sh
```

Downloads ENCODE BAM (ENCFF916YZY, GRCh38-aligned), extracts chr22
properly-paired reads to FASTQ. Requires samtools.

## Run the E2E test

```bash
# Phase 1: db_construct (~5 min)
pixi run nextflow run . -profile pixi,test_e2e_hepg2_chr22 \
    --step db_construct \
    --outdir test_e2e_hepg2/db_construct

# Phase 2: ms_search (~10 min)
# Combine cryptic FASTA with Swiss-Prot for adequate mokapot training
cat test_e2e_hepg2/db_construct/squish/HepG2_chr22_cryptic.fasta \
    tests/data/ms_search/human_swissprot_mini.fasta \
    > test_e2e_hepg2/combined.fasta

pixi run nextflow run . -profile pixi,test_e2e_hepg2_chr22 \
    --step ms_search \
    --search_fasta test_e2e_hepg2/combined.fasta \
    --outdir test_e2e_hepg2/ms_search
```

## Expected results (validated 2026-04-20)

| Phase | Output | Count |
|-------|--------|-------|
| db_construct | Cryptic ORFs (squish/) | 3,354 |
| ms_search — Comet | Target PSMs | 744 (135 cryptic) |
| ms_search — Sage | Target PSMs | 623 (131 cryptic) |

~18-21% of identified peptides are cryptic (non-Swiss-Prot), confirming
the pipeline detects real immunopeptidogenomic signal from chr22 variants
and novel transcripts in the HepG2 immunopeptidome.
