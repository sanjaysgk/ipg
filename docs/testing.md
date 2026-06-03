# Testing `sanjaysgk/ipg`

There are two layers of tests:

1. **Component tests** (`nf-test`) — fast, offline, no data download. Stub-execute every
   local module and subworkflow to verify wiring and channel shapes.
2. **Pipeline test datasets** (`-profile` runs) — execute the real pipeline end-to-end on
   small chr22-subset reference data.

> **Monash M3 / pixi:** prefix all tooling with `pixi run` (e.g. `pixi run nextflow …`,
> `pixi run nf-test …`) — the binaries are provided by the pixi environment, not on `PATH`.
> Elsewhere, drop the `pixi run` prefix and call the tools directly.

---

## 1. Component tests (`nf-test`)

Every local module (`modules/local/*`) and subworkflow (`subworkflows/local/*`) ships a
test under `tests/main.nf.test`. They run as **stub** executions, so they need **no
reference data and no container images** — the stub scripts are trivial shell and the
snapshots are container-independent.

```bash
# run the entire suite
pixi run nf-test test

# run one component by tag
pixi run nf-test test --tag prepare_genome
pixi run nf-test test --tag db_construct

# run a single test file
pixi run nf-test test subworkflows/local/ms_search/tests/main.nf.test

# only the stub tests
pixi run nf-test test --tag stub
```

> **Run stub tests with local execution — do _not_ add `--profile singularity`.**
> `nf-test -stub` still instantiates the declared container; on an offline compute node an
> uncached image fails with `FATAL: the requested image was not found`. The default
> profile (configured in `nf-test.config`) disables all container engines, so the stub
> shell scripts run locally and the snapshots stay portable.

Updating snapshots after an intentional change:

```bash
pixi run nf-test test --tag <name> --update-snapshot
```

`nf-test.config` ignores `modules/nf-core/**` and `subworkflows/nf-core/**` (those are
tested upstream in nf-core/modules) and re-runs tests automatically when
`nextflow.config`, `nf-test.config`, or `conf/test.config` change.

---

## 2. Pipeline test datasets

### 2a. Quick test — `db_construct` on a chr22 bundle

Runs the full RNAseq → cryptic-FASTA path on a chr22 subset of a real human sample.
Completes in ~15 min on a laptop.

```bash
git clone https://github.com/sanjaysgk/ipg.git
cd ipg
bash bin/fetch_test_bundle.sh                       # one-time, downloads the bundle
nextflow run . -profile test,docker --outdir results_test
```

Engine alternatives: `-profile test,singularity` (HPC) or `-profile test,conda`.
Override the bundle location for offline/cached runs with `--test_bundle /local/path`.

### 2b. `ms_search` profiles

```bash
# wiring-only smoke (chr22 FASTA as the search DB) — fast sanity check
nextflow run . -profile test_ms_search,docker --outdir results_ms

# real Orbitrap spectra (nf-core/mhcquant HepG2 fixtures) vs a mini SwissProt DB
nextflow run . -profile test_ms_search_hepg2,docker --outdir results_ms_hepg2
```

Both default to `--ms_engines comet,sage` (MSFragger is skipped so no licensed JAR is
needed) and `--skip_ms2rescore true`.

### 2c. End-to-end — HepG2 chr22 (two phases)

Chains `db_construct → ms_search` on ENCODE HepG2 chr22 RNAseq + nf-core/mhcquant HepG2
immunopeptidome MS. Build the fixture once, then run the two phases:

```bash
bash bin/build_e2e_hepg2_chr22.sh                   # one-time fixture build

# Phase 1 — build the cryptic FASTA from RNAseq
nextflow run . -profile pixi,test_e2e_hepg2_chr22 --step db_construct --outdir results_e2e

# Phase 2 — search MS spectra against the Phase 1 cryptic FASTA
nextflow run . -profile pixi,test_e2e_hepg2_chr22 --step ms_search \
    --search_fasta results_e2e/squish/HepG2_chr22_cryptic.fasta --outdir results_e2e
```

Reference figures (validated 2026-04-20): Phase 1 ≈ 31 processes / 5.5 min / 3354 cryptic
ORFs; Phase 2 ≈ 6 processes / 10 min / ~130 cryptic-peptide PSMs per engine.

### 2d. Synthetic spike-in positive control (known-answer)

The fastest, deterministic check of the **whole** pipeline. A planted cryptic peptide
(`CRYPTICALLY`) must be recovered from synthetic RNAseq into the cryptic FASTA
(`STAR → StringTie → triple_translate → squish`) **and** from synthetic spectra into the
integrated table (`COMET/SAGE → MOKAPOT → MS2RESCORE → INTEGRATE_ENGINES`). No download —
fixtures are generated from a fixed RNG seed (`tests/data/spike/`; regenerate per clone, see
its README).

```bash
# MS-search half (~5 min): search the hand-made FASTA
pixi run nextflow run . -profile pixi,test_spike --outdir results_spike

# db_construct half (~3.5 min): build the cryptic FASTA from synthetic RNAseq
pixi run nextflow run . -profile pixi,test_spike --step db_construct --outdir results_db
```

The **full two-phase e2e** (db_construct → search its *own* cryptic FASTA concatenated with a
background proteome — a cryptic-only DB is too small for Mokapot's FDR) is in
`tests/data/spike/README.md`.

Green = the relevant `integrate/integrated_peptides.tsv` contains `CRYPTICALLY` at peptide
q ≤ 0.01 (engines `['comet','sage']`, q ≈ 0.0033), and db_construct's
`squish/SPIKE_cryptic.fasta` contains it. MSFragger is omitted by default (academic-license
JAR); add it with `--ms_engines comet,sage,msfragger --msfragger_jar /path/to/MSFragger-4.2.jar`.

### Profile reference

| Profile | Step(s) | Data | Purpose |
|---|---|---|---|
| `test` | db_construct | chr22 bundle (fetched) | full RNAseq → cryptic FASTA |
| `test_ms_search` | ms_search | chr22 FASTA as DB | wiring smoke test |
| `test_ms_search_hepg2` | ms_search | HepG2 Orbitrap + mini SwissProt | real MS search |
| `test_e2e_hepg2_chr22` | db_construct → ms_search | HepG2 chr22 | full two-phase e2e |
| `test_full` | db_construct | mirrors `test` | placeholder until full-size data is published |
| `test_spike` | db_construct → ms_search | synthetic (planted `CRYPTICALLY`) | known-answer positive control (full pipeline) |

> The `post_ms` step has no published test profile yet — it is exercised via its
> `nf-test` stub (`--tag post_ms_analysis`) and run on real data downstream of `ms_search`.

---

## Test bundle (`-profile test`)

### Contents

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

Tarball is ~200 MB compressed; `bin/build_test_bundle.sh` reproduces it byte-for-byte
(seqtk seed = 42).

### Regenerating the bundle (maintainer only)

Rarely needed — only when references update (GENCODE bump, new GATK bundle, etc.).

```bash
export REFERENCE_DIR=/path/to/grch38_refs
export VARIANT_DIR=/path/to/gatk_bundle
export FASTQ_R1=/path/to/sample_R1.fastq.gz
export FASTQ_R2=/path/to/sample_R2.fastq.gz
export TEST_BUNDLE_DIR=/path/to/output/bundle

bash bin/build_test_bundle.sh

gh release create v0.1.X-test-data ipg-test-bundle-chr22.tar.gz \
    --title 'Test data bundle vX (chr22)' \
    --notes 'chr22-subset reference + 200k-pair FASTQ for -profile test'
```

Then update the `URL` default in `bin/fetch_test_bundle.sh` if the bundle layout changes.

---

## Troubleshooting

**`StorageFull` on `/tmp` (HPC login nodes).** Many HPC login nodes have a tiny `/tmp`.
Point temp at scratch before fetching/running:

```bash
export TMPDIR=/path/to/your/scratch/tmp
mkdir -p "$TMPDIR"
```

Or submit the test as a small sbatch job (4 CPU, 14 GB, 1 h).

**`FATAL: the requested image was not found` during `nf-test`.** You added
`--profile singularity` to a stub test on an offline node. Drop it — stub tests run
locally (see §1).
