# sanjaysgk/ipg: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

`dev` versions follow the nf-core convention: `<X.Y.Z>dev` (no separator)
indicates active development between releases. The first numeric tag will
be `1.0.0` and will be cut when the pipeline has been validated against the
full ATLANTIS RNA-seq cohort and is ready for publication.

## 1.0.0dev — 2026-04-09

First development snapshot of the nf-core port of the legacy 31-step IPG
cryptic peptide pipeline. Marked as a GitHub pre-release because the
pipeline is not yet feature-complete and the test profile only exercises
a chr22 subset of one D100-liver sample.

### Added

- Six typed nf-core subworkflows under `subworkflows/local/` covering all
  31 steps of the legacy `Select_steps_D122_*.sh` bash pipeline:

  - `align_qc` (steps 1–3): STAR two-pass + samtools sort/index + RSeQC
    `infer_experiment.py`
  - `transcript_assembly` (steps 4–5): StringTie + gffcompare with
    `-R -V -C` for the consensus combined GTF
  - `bam_prep` (steps 6–12): samtools queryname-sort, GATK4 FastqToSam,
    MergeBamAlignment, MarkDuplicates, SplitNCigarReads, ValidateSamFile
  - `bqsr` (steps 13–16): GATK4 BaseRecalibrator (×2), ApplyBQSR,
    AnalyzeCovariates
  - `mutect_calling` (steps 17–23): Mutect2 tumour-only, LearnReadOrientation,
    GetPileupSummaries, CalculateContamination, FilterMutectCalls,
    SelectVariants, curate_vcf
  - `db_construct` (steps 24–31): IndexFeatureFile, FastaAlternateReferenceMaker,
    revert_headers, gff3sort, alt_liftover, gffread, triple_translate, squish

- Twenty-three nf-core/modules installed and pinned via `modules.json`
  for STAR, samtools (sort, index), RSeQC inferexperiment, StringTie,
  gffcompare, gffread, FastQC, MultiQC, and fifteen GATK4 tools.

- Eight local modules under `modules/local/` for tools that are either
  IPG-specific or not yet upstream in `nf-core/modules`:
  `curate_vcf`, `revert_headers`, `alt_liftover`, `triple_translate`,
  `squish`, `gff3sort`, `gatk4_validatesamfile`,
  `gatk4_fastaalternatereferencemaker`.

- `containers/ipg-tools/` reproducible Docker build of the five custom
  IPG C tools, with the source `.c` files bundled directly in the
  repository at `containers/ipg-tools/src/`. Two of the bundled sources
  carry production-quality fixes developed locally at the Li/Purcell labs:

  - `curate_vcf.c` (formerly `curate_vcfV2.c`): larger field buffers
    (`chrom 200`, `ref/alt 4096`, `filter 1024`), dynamic INFO field
    allocation, explicit `free()` cleanup of allocated info strings.
  - `revert_headers.c`: optional 3rd positional argument for the output
    FASTA prefix (writes `<prefix>.fasta` directly instead of the legacy
    hardcoded `tmpc.fasta`). Backward compatible with the upstream 2-arg
    form.
    Both improvements have been pushed upstream to
    `sanjaysgk/immunopeptidogenomics@a09a74c` for separate citeable provenance.

- `pixi` developer environment with `pixi.lock` committed for bit-for-bit
  reproducible toolchain (Nextflow 25.10.4, nf-core 3.5.2, nf-test 0.9.5,
  GATK4 4.6.2, STAR 2.7.11b, samtools 1.23.1, bcftools 1.23.1, StringTie 3.0,
  gffcompare 0.12.10, gffread 0.12.7, RSeQC 5.0.4, FastQC 0.12.1,
  MultiQC 1.33, OpenJDK 17, plus the linting toolchain).

- Five Nextflow profiles:

  - `test` — chr22 subset bundle built locally by `bin/build_test_bundle.sh`
  - `pixi` — run every process from the local pixi env, no containers
  - `singularity` / `apptainer` — pull biocontainers via singularity (HPC default)
  - `docker` — pull biocontainers via docker (laptop / cloud / CI default)
  - `monash` — SLURM executor on Monash M3 `comp` partition under the `xy86`
    project, with shared singularity cache on `/fs04/scratch2/xy86/singularity_cache`

- `--include_variant_peptides` parameter (default `false`) to optionally
  include the alt-reference variant-derived peptide branches (unmasked +
  indel) in the final cryptic peptide database. The default behaviour
  matches the legacy Scull et al. 2021 / D122_Lung run, which was verified
  empirically against the legacy `D122_Lung_squish.log` to fall back on
  the reference branch only. **The flag does NOT switch the variant
  caller to matched tumour-normal mode** — variant calling is always
  performed in tumour-only Mutect2 mode against a gnomAD-style germline
  allele-frequency database, regardless of this flag.

- `bin/build_test_bundle.sh` reproducible chr22 test bundle builder.
  Idempotent, env-var-overridable, derives from the real Monash GRCh38
  reference and a real D100-liver FASTQ sample. Default destination
  `/fs04/scratch2/xy86/sanjay/ipg-test-data/`.

- `nf-test` harness (`nf-test.config`, `tests/nextflow.config`,
  `tests/.nftignore`) with eight per-module stub tests for the local
  modules. Snapshots are committed and verified stable across two
  consecutive runs.

- GitHub Actions workflows:

  - `ci.yml` — parse-only check matrix (Nextflow 24.04.2 and 25.10.4
    × docker and singularity profiles)
  - `build-ipg-tools.yml` — builds and publishes
    `ghcr.io/sanjaysgk/ipg-tools` on push to `main` and on tags matching
    `ipg-tools-v*`. The first image will be cut when the maintainer pushes
    a `ipg-tools-v0.1.0` tag.

- Mermaid workflow diagram in the README showing the six subworkflows,
  their inter-dependencies, and the conditional fan-out gated by
  `--include_variant_peptides`.

- `CITATIONS.md` rewritten to cite Scull et al. 2021 prominently and
  list every tool used by the pipeline (STAR, samtools, StringTie,
  gffcompare, gffread, RSeQC, gff3sort, GATK4, Mutect2, FastQC, MultiQC,
  pixi, bioconda, biocontainers, singularity).

- `.claude/` working memory directory (gitignored) with project-local
  ADRs, plan, status, and progress notes.

### Verified

- All 31 cryptic peptide steps execute end-to-end against the chr22 test
  bundle in approximately 2 minutes on a Monash M3 compute node via
  `pixi run nextflow run . -profile test,pixi --outdir results`.

- The final `cryptic.fasta` deliverable carries `>TCONS_NNNNN_fNpN`
  headers consistent with the legacy `D122_Lung_cryptic.fasta` format
  (TCONS IDs from gffcompare's consensus combined GTF).

- All 8 local-module nf-test stub tests pass and remain
  snapshot-stable across consecutive runs.

- `nextflow config -profile <test|docker|singularity>` parses cleanly
  under both Nextflow 24.04.2 and 25.10.4. CI matrix runs both
  combinations on every push.

- 11 real runtime bugs surfaced by the first end-to-end test run have
  been root-caused and fixed (topic-channel versions misuse, STAR
  `--readFilesCommand zcat` ext.args, GATK4 FastqToSam read-group
  naming, ApplyBQSR cram-vs-bam output extension, BaseRecalibrator
  prefix collision, AnalyzeCovariates R+ggplot2 dependency, Mutect2
  `--f1r2-tar-gz` output, CalculateContamination `--tumor-segmentation`
  output, FilterMutect/SelectVariants prefix collision, curate_vcf
  empty-VCF segfault guard, FastaAlternateReferenceMaker GATK Walker
  arg case, gff3sort versions.yml multiline bug, gffread `-w` arg,
  gffcompare `-C` flag for combined GTF).

### Known limitations

- The `test` profile points at the chr22 subset bundle on Monash M3
  scratch (`/fs04/scratch2/xy86/sanjay/ipg-test-data/`). The bundle is
  not committed to the repository or accessible from GitHub Actions
  runners. CI runs only `nextflow config` parse checks; the full
  end-to-end test must be run locally on Monash. Publishing the test
  bundle to a public location (e.g. fork of `nf-core/test-datasets`)
  is a future task.

- `nf-core pipelines lint` is configured with several intentional
  exclusions documented inline in `.nf-core.yml`:

  - `files_exist`: the nf-core template's `nf-test.yml` workflow,
    `actions/get-shards`, `actions/nf-test`, and `tests/default.nf.test`
    are not included (we use a different nf-test harness layout).
  - `files_unchanged`: `.gitattributes`, `.prettierrc.yml`, the issue/PR
    templates, `CONTRIBUTING.md`, and `linting_comment.yml` were updated
    as part of the bulk `nf-core/ipg → sanjaysgk/ipg` identity rename.
  - `readme`: nf-core template badges removed in the rewritten README.
  - `manifest.version`: allowed to contain `dev` until the first
    numeric release.

- The pipeline has not yet been run end-to-end against a full real
  sample (e.g. the D100-liver source the test bundle is derived from).
  This is the next milestone toward the `1.0.0` release.

- The `ghcr.io/sanjaysgk/ipg-tools` container image referenced by the
  five IPG-tool local modules has not been published yet. Pushing a
  `ipg-tools-v0.1.0` git tag will trigger the GHA workflow that publishes
  the first image. Until then, only the `pixi` profile (which uses tools
  from the local PATH) can run end-to-end without modification.

### Pipeline-runtime tool versions pinned in `pixi.lock`

| Tool       | Version |
| ---------- | ------- |
| Nextflow   | 25.10.4 |
| nf-core    | 3.5.2   |
| nf-test    | 0.9.5   |
| OpenJDK    | 17.0.18 |
| STAR       | 2.7.11b |
| samtools   | 1.23.1  |
| bcftools   | 1.23.1  |
| GATK4      | 4.6.2.0 |
| StringTie  | 3.0.0   |
| gffcompare | 0.12.10 |
| gffread    | 0.12.7  |
| RSeQC      | 5.0.4   |
| FastQC     | 0.12.1  |
| MultiQC    | 1.33    |
| seqtk      | 1.5     |

### Authors

- **Sanjay SG Krishna** ([@sanjaysgk](https://github.com/sanjaysgk)) —
  pipeline port; Li Lab, Monash University
- **Kate Scull** — original IPG method, custom C tools
  ([`kescull/immunopeptidogenomics`](https://github.com/kescull/immunopeptidogenomics));
  Purcell Lab, Monash University
- **Chen Li** — supervision; Li Lab, Monash University
- **Anthony W. Purcell** — supervision; Purcell Lab, Monash University

### Citation

If you use this version, please cite the original method paper:

> Scull KE, Pandey K, Ramarathinam SH, Purcell AW.
> _Immunopeptidogenomics: harnessing RNA-seq to illuminate the dark immunopeptidome._
> Mol Cell Proteomics. 2021;20:100143.
> [doi.org/10.1016/j.mcpro.2021.100143](https://doi.org/10.1016/j.mcpro.2021.100143)

A complete reference list for every tool in the pipeline lives in
[`CITATIONS.md`](CITATIONS.md).
