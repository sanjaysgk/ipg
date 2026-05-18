# Testing `sanjaysgk/ipg`

## Quick test — any host with Nextflow + a container engine

```bash
git clone https://github.com/sanjaysgk/ipg.git
cd ipg
bash bin/fetch_test_bundle.sh
nextflow run . -profile test,docker --outdir results_test
```

Runs the full pipeline on a chr22-subset of a real human RNAseq sample. Completes in ~15 min on a laptop. Outputs land in `results_test/`.

### Engine alternatives

- `-profile test,singularity` — on HPC
- `-profile test,conda` — if no container engine

### HPC login-node pitfall

If `bin/fetch_test_bundle.sh` or pipeline activation fails with `StorageFull` on `/tmp`, your HPC login node has a tiny `/tmp`. Point it at scratch:

```bash
export TMPDIR=/path/to/your/scratch/tmp
mkdir -p "$TMPDIR"
```

Then re-run. Or submit the test as a small sbatch job (4 CPU, 14 GB, 1 h).

## Regenerating the bundle (maintainer only)

This is rarely needed — only when references update (GENCODE version bump, new GATK bundle, etc.).

```bash
export REFERENCE_DIR=/path/to/grch38_refs
export VARIANT_DIR=/path/to/gatk_bundle
export FASTQ_R1=/path/to/sample_R1.fastq.gz
export FASTQ_R2=/path/to/sample_R2.fastq.gz
export TEST_BUNDLE_DIR=/path/to/output/bundle

bash bin/build_test_bundle.sh
```

Produces `ipg-test-bundle-chr22.tar.gz` adjacent to `$TEST_BUNDLE_DIR`. Upload:

```bash
gh release create v0.1.X-test-data ipg-test-bundle-chr22.tar.gz \
    --title 'Test data bundle vX (chr22)' \
    --notes 'chr22-subset reference + 200k-pair FASTQ for -profile test'
```

Then update `bin/fetch_test_bundle.sh:URL` default and `conf/test.config` description if the bundle layout changes.

## What's in the bundle

```
.test-bundle/
├── samplesheet_test.csv
├── fastq/
│   ├── test_R1.fastq.gz             200 k read pairs (seqtk -s 42)
│   └── test_R2.fastq.gz
├── reference/
│   ├── GRCh38.chr22.fa + .fai + .dict
│   ├── gencode.v44.chr22.gtf
│   ├── gencode_assembly.chr22.bed   for RSeQC
│   └── star_index_chr22/            STAR pre-built
└── variant_calling/
    ├── dbsnp138.chr22.vcf.gz + .tbi
    ├── known_indels.chr22.vcf.gz + .tbi
    ├── Mills.chr22.vcf.gz + .tbi
    └── small_exac_common_3.chr22.vcf.gz + .tbi
```

Tarball is ~200 MB compressed. Source pipeline `bin/build_test_bundle.sh` reproduces it byte-for-byte from sources (seqtk seed = 42).
