# Test fixtures

Each subdirectory here is a **self-contained test dataset** for one
`--step` of the pipeline. Every subdir MUST have its own `README.md`
documenting:

1. **Source** — where the data came from (PRIDE accession, patient ID,
   public repo URL)
2. **Contents** — every file, its size, what it represents
3. **Launch command** — exact `nextflow run ...` invocation
4. **Expected outputs** — what a successful run produces, for nf-test
   regression assertions

## Current fixtures

| Directory | Step | Source | Size | Licensed tools? |
|---|---|---|---|---|
| [`ipg-db-constructions/`](ipg-db-constructions/README.md) | `--step db_construct` | D100_liver chr22 subset, GRCh38 / GENCODE v44 | ~545 MB | none |
| [`ms_search/`](ms_search/README.md) | `--step ms_search` | HepG2 from nf-core/mhcquant (PRIDE PXD009752) | ~160 MB | MSFragger optional |
| [`post_ms/`](post_ms/README.md) | `--step post_ms` | Hand-built synthetic — 10 cryptic + 12 UniProt peptides | ~10 KB | none |

## Adding a new fixture

When adding `tests/data/<newthing>/`:

1. **Create `<newthing>/README.md` first** — describe what the data is,
   where it came from, what it tests. The README should read well in
   isolation; anyone stumbling into the directory should understand its
   purpose in under a minute.
2. **Add a fetch script** under `bin/build_test_<newthing>_bundle.sh` if
   the data is downloadable / regenerable. Real data → fetch on demand,
   never commit binaries >10 MB.
3. **Update `.gitignore`** to exclude large binary files (mzML, FASTQ,
   VCF, STAR indexes, etc.) while keeping the README + samplesheets
   tracked.
4. **Register a test profile** in `conf/test_<newthing>.config` and
   `nextflow.config` so the fixture can be exercised with
   `-profile pixi,test_<newthing>`.
5. **Update this index** adding a row to the table above.

## Why fixtures are per-step, not monolithic

Separate subdirs keep the test matrix tight:

- A PR that only touches `--step ms_search` doesn't need to download
  the 545 MB db_construct bundle.
- CI can run `test_post_ms` (takes 30 s) on every push without pulling
  any large binaries — the synthetic fixture is tiny enough to commit.
- Each fixture evolves independently with the step it tests.

## Launching the full test matrix

```bash
# fast (synthetic only, no downloads)
pixi run nextflow run . -profile pixi,test_post_ms

# medium (downloads 160 MB, ~10 min)
bash bin/build_test_ms_search_bundle.sh
pixi run nextflow run . -profile pixi,test_ms_search_hepg2

# slow (needs chr22 bundle on Monash / local copy, ~15 min)
pixi run nextflow run . -profile pixi,test
```
