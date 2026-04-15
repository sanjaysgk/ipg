# sanjaysgk/ipg

[![GitHub Actions CI Status](https://github.com/sanjaysgk/ipg/actions/workflows/ci.yml/badge.svg)](https://github.com/sanjaysgk/ipg/actions/workflows/ci.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A525.10.4-23aa62.svg)](https://www.nextflow.io/)
[![pixi](https://img.shields.io/badge/dev_env-pixi-yellow.svg)](https://pixi.sh)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Immunopeptidogenomics — cryptic peptide database construction from RNA-seq.**
> An [nf-core](https://nf-co.re)-style port of the 31-step bash pipeline used by the
> [Li Lab](https://research.monash.edu/en/persons/chen-li) and
> [Purcell Lab](https://research.monash.edu/en/persons/anthony-purcell) at Monash University,
> implementing the cryptic peptide discovery method described in
> [Scull et al. _Mol Cell Proteomics_ 2021](https://doi.org/10.1016/j.mcpro.2021.100143).

## What this pipeline does

`sanjaysgk/ipg` turns paired-end RNA-seq reads into a **cryptic peptide search database**
suitable for downstream MS/MS immunopeptidomics searches. Starting from FASTQs, it:

1. Aligns reads with **STAR** (two-pass) and infers strandedness with **RSeQC**.
2. Assembles novel transcripts with **StringTie** and reconciles them with the reference
   annotation via **gffcompare** (`-R -V -C` for the consensus combined GTF).
3. Builds a variant-calling-ready BAM via the GATK4 RNA-seq Best Practices chain
   (FastqToSam → MergeBamAlignment → MarkDuplicates → SplitNCigarReads).
4. Recalibrates base qualities with two passes of **BaseRecalibrator + ApplyBQSR**.
5. Calls somatic variants with **Mutect2** in tumour-only mode against a gnomAD-style
   germline allele-frequency database, then **CalculateContamination** /
   **FilterMutectCalls** / **SelectVariants** for clean PASS-only sites-only VCFs.
6. Builds the cryptic peptide FASTA database via the IPG custom C tools
   ([`curate_vcf`](https://github.com/sanjaysgk/immunopeptidogenomics),
   [`revert_headers`](https://github.com/sanjaysgk/immunopeptidogenomics),
   [`alt_liftover`](https://github.com/sanjaysgk/immunopeptidogenomics),
   [`triple_translate`](https://github.com/sanjaysgk/immunopeptidogenomics),
   [`squish`](https://github.com/sanjaysgk/immunopeptidogenomics)) plus
   `gff3sort` and `gffread` from bioconda.

Optionally, a **post-MS analysis** step (`--step post_ms`) runs the two-phase
`db_compare` + `origins` workflow on PEAKS (or other search engine) results to
identify and annotate cryptic-only peptides.

An **MS search** step (`--step ms_search`) runs up to four open-source search
engines (**MSFragger**, **Comet**, **Sage**, **PEAKS**) in parallel against the
cryptic FASTA, applies **mokapot** FDR control per engine, rescores PSMs with
**MS2Rescore**, and merges results at 1% peptide-level FDR. Optional
immunoinformatics gates — `--run_netmhcpan`, `--run_netmhciipan`,
`--run_gibbscluster`, `--run_flashlfq`, `--run_blastp_host` — pick up the
integrated peptide table and emit a per-sample HTML report summarising HLA
binding, motif clusters, quantification, and host-background hits. See
[`docs/usage.md`](docs/usage.md) for the full invocation.

The 31 legacy steps are grouped into **seven typed nf-core subworkflows**:

```mermaid
%%{init: {
  "theme": "base",
  "themeVariables": {
    "primaryColor":         "#e0f2fe",
    "primaryTextColor":     "#0c4a6e",
    "primaryBorderColor":   "#0369a1",
    "lineColor":            "#0369a1",
    "secondaryColor":       "#fef3c7",
    "tertiaryColor":        "#dcfce7",
    "fontFamily":           "ui-sans-serif, system-ui, -apple-system, sans-serif",
    "fontSize":             "14px"
  }
}}%%
flowchart LR
    %% ----- nodes -----
    INPUT(["fa:fa-dna <b>Paired-end FASTQ</b><br/>+ samplesheet.csv"]):::input

    subgraph QC1["1. ALIGN_QC &nbsp; <span style='color:#64748b;font-size:11px'>(steps 1–3)</span>"]
        direction TB
        STAR1["STAR 2-pass align"]
        SORT1["samtools sort + index"]
        RSEQC["RSeQC infer_experiment"]
        STAR1 --> SORT1 --> RSEQC
    end

    subgraph TA["2. TRANSCRIPT_ASSEMBLY &nbsp; <span style='color:#64748b;font-size:11px'>(steps 4–5)</span>"]
        direction TB
        STRINGTIE["StringTie"]
        GFFCOMPARE["gffcompare<br/>(-R -V -C)"]
        STRINGTIE --> GFFCOMPARE
    end

    subgraph BP["3. BAM_PREP &nbsp; <span style='color:#64748b;font-size:11px'>(steps 6–12)</span>"]
        direction TB
        SORTQN["samtools sort -n"]
        F2S["GATK4 FastqToSam"]
        MERGE["GATK4 MergeBamAlignment"]
        MARKDUP["GATK4 MarkDuplicates"]
        SPLIT["GATK4 SplitNCigarReads"]
        VAL["GATK4 ValidateSamFile<br/><span style='color:#64748b;font-size:11px'>(audit only)</span>"]
        SORTQN --> F2S --> MERGE --> MARKDUP --> SPLIT --> VAL
    end

    subgraph BQ["4. BQSR &nbsp; <span style='color:#64748b;font-size:11px'>(steps 13–16)</span>"]
        direction TB
        BR1["BaseRecalibrator (1st)"]
        APPLY["ApplyBQSR"]
        BR2["BaseRecalibrator (2nd)"]
        AC["AnalyzeCovariates<br/><span style='color:#64748b;font-size:11px'>(QC plot)</span>"]
        BR1 --> APPLY --> BR2 --> AC
    end

    subgraph MC["5. MUTECT_CALLING &nbsp; <span style='color:#64748b;font-size:11px'>(steps 17–23)</span>"]
        direction TB
        M2["Mutect2 tumour-only"]
        LOM["LearnReadOrientationModel"]
        GPS["GetPileupSummaries"]
        CC["CalculateContamination"]
        FMC["FilterMutectCalls"]
        SV["SelectVariants<br/>(PASS, sites-only)"]
        CV["curate_vcf<br/><span style='color:#64748b;font-size:11px'>(IPG)</span>"]
        M2 --> LOM
        M2 --> GPS --> CC
        LOM --> FMC
        CC --> FMC
        M2 --> FMC --> SV --> CV
    end

    subgraph DC["6. DB_CONSTRUCT &nbsp; <span style='color:#64748b;font-size:11px'>(steps 24–31)</span>"]
        direction TB
        GFF3["gff3sort"]
        GFR["gffread"]
        TT["triple_translate<br/><span style='color:#64748b;font-size:11px'>(IPG)</span>"]
        VAR{{"--include_variant_peptides<br/>true?"}}:::flag
        IFR["IndexFeatureFile + FARM<br/>+ revert_headers + alt_liftover<br/><span style='color:#64748b;font-size:11px'>(unmasked + indel branches)</span>"]
        SQUISH["squish<br/><span style='color:#64748b;font-size:11px'>(IPG)</span>"]
        GFF3 --> GFR --> TT --> SQUISH
        VAR -->|yes| IFR --> SQUISH
    end

    QC2["FastQC"]:::qc
    MQC["fa:fa-chart-bar <b>MultiQC report</b>"]:::report
    OUT(["fa:fa-database <b>cryptic_peptide.fasta</b><br/>(MS/MS-ready DB)"]):::deliverable

    %% ----- edges -----
    INPUT --> QC2
    INPUT --> STAR1
    SORT1 --> STRINGTIE
    SORT1 --> SORTQN
    SPLIT --> BR1
    APPLY --> M2
    GFFCOMPARE --> GFF3
    SV --> IFR
    SQUISH --> OUT

    QC2 --> MQC
    RSEQC --> MQC
    MARKDUP --> MQC
    AC --> MQC
    GFFCOMPARE --> MQC
    FMC --> MQC

    %% ----- classes -----
    classDef input        fill:#fef9c3,stroke:#ca8a04,stroke-width:2px,color:#713f12
    classDef deliverable  fill:#bbf7d0,stroke:#16a34a,stroke-width:3px,color:#14532d
    classDef report       fill:#e9d5ff,stroke:#9333ea,stroke-width:2px,color:#581c87
    classDef qc           fill:#fce7f3,stroke:#db2777,stroke-width:2px,color:#831843
    classDef flag         fill:#fef3c7,stroke:#d97706,stroke-width:2px,color:#92400e

    style QC1 fill:#f0f9ff,stroke:#0284c7,stroke-width:2px
    style TA  fill:#f0fdf4,stroke:#16a34a,stroke-width:2px
    style BP  fill:#fef2f2,stroke:#dc2626,stroke-width:2px
    style BQ  fill:#fefce8,stroke:#ca8a04,stroke-width:2px
    style MC  fill:#faf5ff,stroke:#9333ea,stroke-width:2px
    style DC  fill:#fff7ed,stroke:#ea580c,stroke-width:2px
```

## Quick start

### 1. Install the dev environment

The repository ships a [pixi](https://pixi.sh) project that pins every tool
(nextflow 25.10.4, nf-core 3.5.2, nf-test 0.9.5, GATK4, STAR, samtools, bcftools,
stringtie, gffcompare, gffread, RSeQC, FastQC, MultiQC, OpenJDK 17, plus the
linting toolchain) to bit-for-bit reproducible versions via `pixi.lock`.

```bash
git clone https://github.com/sanjaysgk/ipg.git
cd ipg
pixi install
```

If you don't have pixi: `curl -fsSL https://pixi.sh/install.sh | bash`.

### 2. Build the chr22 test bundle (~5 minutes, one time)

```bash
pixi run bash bin/build_test_bundle.sh
```

By default this writes to `/fs04/scratch2/xy86/sanjay/ipg-test-data/`. Override
on a non-Monash machine:

```bash
TEST_BUNDLE_DIR=/some/path \
REFERENCE_DIR=/path/to/GRCh38 \
VARIANT_DIR=/path/to/variant-resources \
FASTQ_R1=/path/to/sample_R1.fq.gz \
FASTQ_R2=/path/to/sample_R2.fq.gz \
pixi run bash bin/build_test_bundle.sh
```

The bundle is **never committed to git** — only the build script is. See
[`bin/build_test_bundle.sh`](bin/build_test_bundle.sh) for the exact recipe.

### 3. Run the pipeline against the test bundle

| Profile combination                | Container engine             | Use when                                                                                  |
| ---------------------------------- | ---------------------------- | ----------------------------------------------------------------------------------------- |
| `-profile test,pixi`               | none — tools from local PATH | Fastest iteration on a workstation/HPC compute node where pixi can run all tools natively |
| `-profile test,singularity`        | apptainer / singularity      | Standard HPC; pulls biocontainers from `community.wave.seqera.io`                         |
| `-profile test,docker`             | docker                       | Laptops, cloud, GitHub Actions CI                                                         |
| `-profile test,monash,singularity` | singularity                  | Monash M3 SLURM cluster (`comp` partition, `xy86` account)                                |

```bash
# Fastest (no containers, all tools from pixi):
pixi run nextflow run . -profile test,pixi --outdir results
```

The pipeline takes **~2 minutes** end-to-end against the chr22 test bundle on a
single Linux box and produces the cryptic peptide database at
`results/squish/<sample>_cryptic.fasta`.

### 4. Run on real data

Create a samplesheet `samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,strandedness
D100_liver,/path/to/D100-liver_R1.fastq.gz,/path/to/D100-liver_R2.fastq.gz,reverse
D101_pancreas,/path/to/D101-pancreas_R1.fastq.gz,/path/to/D101-pancreas_R2.fastq.gz,reverse
```

Then run with explicit reference paths:

```bash
pixi run nextflow run . \
    -profile monash,singularity \
    --input            samplesheet.csv \
    --outdir           /fs04/scratch2/xy86/sanjay/ipg-results \
    --fasta            /path/to/GRCh38.primary_assembly.genome.fa \
    --fasta_fai        /path/to/GRCh38.primary_assembly.genome.fa.fai \
    --fasta_dict       /path/to/GRCh38.primary_assembly.genome.dict \
    --gtf              /path/to/gencode.v44.primary_assembly.annotation.gtf \
    --star_index       /path/to/human_genome_index_GRCh38 \
    --rseqc_bed        /path/to/gencode_assembly.bed \
    --dbsnp            /path/to/dbsnp138.vcf \
    --dbsnp_tbi        /path/to/dbsnp138.vcf.idx \
    --known_indels     /path/to/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known_indels_tbi /path/to/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi \
    --mills            /path/to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --mills_tbi        /path/to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi \
    --germline_resource     /path/to/small_exac_common_3.hg38.vcf.gz \
    --germline_resource_tbi /path/to/small_exac_common_3.hg38.vcf.gz.tbi
```

## The `--include_variant_peptides` flag

By default the cryptic peptide DB is built from the **reference assembly only**,
matching the legacy Scull et al. 2021 / D122_Lung run (verified empirically by
inspecting the legacy `squish.log`). Pass `--include_variant_peptides true` to
also fold the alt-reference variant-derived peptide branches (unmasked + indel)
into the final database:

```bash
pixi run nextflow run . -profile monash,singularity \
    --input samplesheet.csv \
    --include_variant_peptides true \
    [other reference args]
```

> [!IMPORTANT] > `--include_variant_peptides` does **not** switch the variant caller to matched
> tumour-normal mode. Variant calling is **always** performed in tumour-only
> Mutect2 mode against a gnomAD-style germline allele-frequency database; this
> pipeline does **not** support matched tumour-normal calling. The flag controls
> only whether the discovered variants get folded into the final cryptic peptide
> DB. Set `true` only when the sample is expected to harbour biologically
> meaningful somatic variants (e.g. tumour tissue, hypermutated cell lines,
> MMR-deficient samples). Leave at `false` for normal tissue, cell lines, or any
> sample where variant peptides would mostly add noise. This is a pipeline-level
> flag — to mix modes for a heterogeneous cohort, run the pipeline twice.

## Post-MS analysis (`--step post_ms`)

After running the DB construction pipeline, the cryptic peptide FASTA is searched
against MS/MS data using PEAKS Online (or another search engine such as
MSFragger, Comet, or Sage). The resulting PSM CSVs are then analysed with the
two-phase `db_compare` + `origins` workflow
([Scull et al. 2021](https://doi.org/10.1016/j.mcpro.2021.100143)):

```
Phase 1: db_compare_v2.R  →  cryptic_only.txt
         origins -s        →  origins_discard.txt + origins_unconventional.txt

Phase 2: db_compare_v2.R (with -j discard -u unconventional)
         →  unambiguous_unconventional.txt
         origins (full Ensembl mode)  →  deep origin annotation
```

### Post-MS samplesheet

Create a CSV with one row per sample:

```csv
sample,cryptic_psm_csv,uniprot_psm_csv,cryptic_decoy_score,uniprot_decoy_score
D122_liver,/path/to/D122_Liver_Cryptic_DB.db.psms.csv,/path/to/D122_Liver_Uniprot_DB.db.psms.csv,44.43,36.48
D101_heart,/path/to/D101_Heart_Cryptic_DB.db.psms.csv,/path/to/D101_Heart_Uniprot_DB.db.psms.csv,3.64,3.25
```

The `cryptic_decoy_score` and `uniprot_decoy_score` columns are the `-10lgP`
decoy thresholds from the respective PEAKS searches.

### Running post-MS analysis

```bash
pixi run nextflow run . -profile pixi,monash \
    --step post_ms \
    --post_ms_input       post_ms_samplesheet.csv \
    --uniprot_fasta       /path/to/uniprotkb_human_canonical_isoform.fasta \
    --transcriptome_fasta /path/to/D122_liver_transcriptome.fa \
    --prefix_tracking     /path/to/prefix.tracking \
    --outdir              results_post_ms
```

The `--transcriptome_fasta` and `--prefix_tracking` files are outputs from the
DB construction step (`gffread` and `gffcompare` respectively). You can find them
in the pipeline output directory from the previous run.

### Post-MS output

```
results_post_ms/
├── post_ms/
│   ├── phase1/
│   │   ├── db_compare/     Phase 1 cryptic-only peptide lists + plots
│   │   └── origins/        Phase 1 origins (simple mode) — discard + unconventional lists
│   └── phase2/
│       ├── db_compare/     Phase 2 refined unambiguous unconventional peptides
│       └── origins/        Phase 2 origins (full Ensembl annotation)
└── pipeline_info/
```

## Output

```
results/
├── star/                    STAR alignment BAMs + .Log.final.out
├── samtools/                sorted/indexed BAMs
├── rseqc/                   strandedness inference reports
├── stringtie/               assembled transcript GTFs
├── gffcompare/              prefix.combined.gtf, prefix.tracking, .stats
├── gatk4/                   Mutect2 / BQSR / contamination tables
├── revert/                  alt-reference FASTAs (only when --include_variant_peptides=true)
├── gff3sort/                sorted assembly GTF
├── tt/                      triple_translate per-branch peptide FASTAs
├── squish/
│   └── <sample>_cryptic.fasta   ← THE DELIVERABLE
├── multiqc/
│   └── multiqc_report.html      aggregated QC report
└── pipeline_info/
    ├── execution_report_<timestamp>.html
    ├── execution_timeline_<timestamp>.html
    └── pipeline_dag_<timestamp>.html
```

## Profiles

| Profile                     | Purpose                                                                                                  |
| --------------------------- | -------------------------------------------------------------------------------------------------------- |
| `test`                      | Use the chr22 test bundle (built by `bin/build_test_bundle.sh`)                                          |
| `pixi`                      | Run every process from the local pixi env, no containers                                                 |
| `singularity` / `apptainer` | Pull biocontainers via singularity/apptainer (HPC default)                                               |
| `docker`                    | Pull biocontainers via docker (laptop / cloud / CI default)                                              |
| `monash`                    | SLURM executor on the Monash M3 `comp` partition under the `xy86` project, with shared singularity cache |

## Architecture

```
sanjaysgk/ipg/
├── main.nf                              entry point
├── workflows/ipg.nf                     main workflow — chains the 6 subworkflows
├── subworkflows/local/
│   ├── align_qc/                        steps 1-3
│   ├── transcript_assembly/             steps 4-5
│   ├── bam_prep/                        steps 6-12
│   ├── bqsr/                            steps 13-16
│   ├── mutect_calling/                  steps 17-23
│   ├── db_construct/                    steps 24-31 (branches on --include_variant_peptides)
│   └── post_ms_analysis/               --step post_ms: 2-phase db_compare + origins
├── modules/
│   ├── nf-core/                         23 upstream nf-core modules (STAR, samtools, GATK4, etc.)
│   └── local/                           10 local modules:
│       ├── curate_vcf/                  IPG custom C tool (kescull)
│       ├── revert_headers/              IPG custom C tool (kescull)
│       ├── alt_liftover/                IPG custom C tool (kescull)
│       ├── triple_translate/            IPG custom C tool (kescull)
│       ├── squish/                      IPG custom C tool (kescull)
│       ├── origins/                     IPG custom C tool — peptide origin annotation (kescull)
│       ├── db_compare/                  R script — cryptic vs UniProt PSM comparison (kescull)
│       ├── gff3sort/                    bioconda gff3sort wrapper
│       ├── gatk4_validatesamfile/       missing-from-upstream wrapper
│       └── gatk4_fastaalternatereferencemaker/  missing-from-upstream wrapper
├── containers/
│   └── ipg-tools/                       Reproducible Docker build of the kescull C tools
│                                        from a pinned commit SHA → ghcr.io/sanjaysgk/ipg-tools
├── conf/
│   ├── base.config
│   ├── modules.config                   per-process ext.args / ext.prefix / errorStrategy
│   ├── test.config                      chr22 test bundle paths
│   └── monash.config                    Monash M3 SLURM
├── bin/
│   └── build_test_bundle.sh             reproducible chr22 test bundle builder
├── tests/                               nf-test pipeline-level tests
├── pixi.toml                            dev env definition (committed)
└── pixi.lock                            dev env lockfile (committed, bit-reproducible)
```

## Authors

- **Sanjay SG Krishna** ([@sanjaysgk](https://github.com/sanjaysgk)) — pipeline port,
  Li Lab, Monash University
- **Kate Scull** — original IPG method + custom C tools, Purcell Lab, Monash University
  (see [`kescull/immunopeptidogenomics`](https://github.com/kescull/immunopeptidogenomics))
- **Chen Li** — supervision, Li Lab, Monash University
- **Anthony W. Purcell** — supervision, Purcell Lab, Monash University

## Citation

If you use `sanjaysgk/ipg` in your research, please cite the original method paper:

> **Scull KE, Pandey K, Ramarathinam SH, Purcell AW.** > _Immunopeptidogenomics: harnessing RNA-seq to illuminate the dark immunopeptidome._
> Mol Cell Proteomics. 2021;20:100143.
> doi: [10.1016/j.mcpro.2021.100143](https://doi.org/10.1016/j.mcpro.2021.100143)

A full reference list for every tool in the pipeline lives in [`CITATIONS.md`](CITATIONS.md).

This pipeline is built on top of [Nextflow](https://www.nextflow.io) and the
[nf-core](https://nf-co.re) framework:

> **The nf-core framework for community-curated bioinformatics pipelines.**
> Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU,
> Di Tommaso P, Nahnsen S.
> _Nat Biotechnol._ 2020 Mar;38(3):276-278.
> doi: [10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x).

## License

MIT — see [LICENSE](LICENSE).
