# db_construct test fixtures — chr22 subset

Pre-built chr22-only bundle for exercising `--step db_construct`
end-to-end without needing the full GRCh38 + GATK resource tree. Derived
from a real patient sample so the pipeline produces realistic cryptic
peptide output at every stage — not synthetic placeholders.

## Source

| | |
|---|---|
| Reference genome | GENCODE v44, GRCh38 primary assembly, **chr22 only** |
| Patient sample | D100_liver — first-batch Monash ATLANTIS cohort |
| RNA-seq subsetting | chr22 reads extracted with `samtools view` from the full D100 BAM, then re-exported to paired FASTQ |
| Build script | `bin/build_test_bundle.sh` (regenerates from the full GRCh38 tree) |
| Origin path on M3 | Copied from `/fs04/scratch2/xy86/sanjay/ipg-test-data/` on 2026-04-15 |

## Contents (~545 MB total)

```
ipg-db-constructions/
├── samplesheet_test.csv        # nf-core ipg input (D100_liver_chr22)
├── fastq/                      # 36 MB — chr22-only paired reads
│   ├── test_R1.fastq.gz
│   └── test_R2.fastq.gz
├── reference/                  # 486 MB — chr22 genome + annotation + STAR index
│   ├── GRCh38.chr22.fa
│   ├── GRCh38.chr22.fa.fai
│   ├── GRCh38.chr22.dict
│   ├── gencode.v44.chr22.gtf
│   ├── gencode_assembly.chr22.bed     # RSeQC strandedness BED
│   └── star_index_chr22/              # pre-built STAR 2.7.11b index
├── variant_calling/            # 20 MB — GATK4 resource chr22 subsets
│   ├── dbsnp138.chr22.vcf.gz (+ .tbi)
│   ├── known_indels.chr22.vcf.gz (+ .tbi)
│   ├── Mills.chr22.vcf.gz (+ .tbi)
│   └── small_exac_common_3.chr22.vcf.gz (+ .tbi)   # Mutect2 germline resource
└── _prefilter/                 # 2.7 MB — intermediate artefacts from the
                                #          bundle build (safe to ignore)
```

## Not committed to git

These files together weigh ~545 MB which exceeds GitHub limits. `.gitignore`
excludes `*.fa`, `*.fastq.gz`, `*.vcf.gz*`, and the whole STAR index. To
populate this directory on a fresh clone, either:

1. Copy from Monash M3 reference location:
   ```
   cp -r /fs04/scratch2/xy86/sanjay/ipg-test-data/* \
         tests/data/ipg-db-constructions/
   ```
2. Regenerate from the full GRCh38 tree:
   ```
   bash bin/build_test_bundle.sh
   ```

## Launch the test

```bash
pixi run nextflow run . -profile pixi,test \
    --outdir test_db_construct_results
```

Runtime: ~15 minutes on a laptop, ~5 minutes on a Monash M3 comp node.

## Expected outputs (for nf-test assertions)

After a successful run, the pipeline produces:

| Output | Where | Rows |
|---|---|---|
| Calibrated chr22 BAM | `samtools/` | — |
| Novel transcript GTF | `stringtie/` | ≥1 |
| gffcompare combined GTF | `gffcompare/` | ≥1 |
| Mutect2 filtered VCF | `gatk4/filtermutectcalls/` | ≥1 PASS variant |
| Curated indel + unmasked VCFs | `curate/` | 2 files |
| 3-frame translated FASTAs (ref + indel + unmasked) | `tt/` | 3 files, each ≥100 proteins |
| **Final cryptic FASTA** (after squish dedupe) | `squish/` | ≥5000 unique cryptic ORFs |

The test fixture was validated in ipg SLURM job 54796415 (2026-04-14):
D100_liver_chr22 db_construct completed in 12 min producing
`D100_liver_chr22_cryptic.fasta` with ~6000 ORFs.

## Why chr22 (not the full genome)?

- Full genome db_construct takes 2-3 hours on HPC; not viable in CI.
- chr22 has enough variant density to hit every pipeline module with
  non-empty output (the whole point of the fixture).
- STAR index for chr22 is ~400 MB vs ~30 GB full genome — manageable.
- Matches the nf-core convention of shipping chromosome-subset test data.
