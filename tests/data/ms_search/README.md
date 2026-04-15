# ms_search test fixtures — HepG2 from nf-core/mhcquant

Public LC-MS/MS immunopeptidomics data used to exercise the
`--step ms_search` subworkflow end-to-end.

## Source

| | |
|---|---|
| PRIDE accession | [PXD009752](https://www.ebi.ac.uk/pride/archive/projects/PXD009752) |
| Instrument | Thermo Q Exactive (Orbitrap) |
| Sample | HepG2 hepatocellular carcinoma cell line |
| Curated test dataset | [nf-core/test-datasets @ mhcquant](https://github.com/nf-core/test-datasets/tree/mhcquant) |
| Total size | ~160 MB |

## HepG2 HLA alleles (for netMHCpan-style tests)

**Class I:** `A*02:01;A*24:02;B*35:01;B*51:08;C*04:01;C*16:02`
**Class II:** `HLA-DRB1*13:02;HLA-DRB1*16:02;...` (see `HepG2_allele_sheet.tsv`)

## Files (downloaded by `bin/build_test_ms_search_bundle.sh`)

| File | Purpose |
|---|---|
| `HepG2_rep1_small.mzML` | Replicate 1 MS data (48 MB) |
| `HepG2_rep2_small.mzML` | Replicate 2 (48 MB) |
| `HepG2_rep3_small.mzML` | Replicate 3 (47 MB) |
| `human_swissprot.fasta` | Human Swiss-Prot DB (14 MB) |
| `HepG2_allele_sheet.tsv` | HLA allele reference |
| `HepG2_sample_sheet.tsv` | nf-core/mhcquant original samplesheet |
| `samplesheet.csv` | ipg-format samplesheet for the 3 mzMLs |

**Not committed.** Fetched on demand. `.gitignore` excludes `*.mzML` and
`*.fasta` under this directory.

## Launch

```bash
# one-time fetch
bash bin/build_test_ms_search_bundle.sh

# fast path — Comet + Sage only (no licensed JAR needed)
pixi run nextflow run . -profile pixi,test_ms_search_hepg2 \
    --outdir results/test_ms_search_hepg2

# full path — include MSFragger (requires JAR)
pixi run nextflow run . -profile pixi,test_ms_search_hepg2 \
    --ms_engines msfragger,comet,sage \
    --msfragger_jar /path/to/MSFragger-4.1.jar \
    --outdir results/test_ms_search_hepg2
```

## Expected outputs (non-empty assertions)

After a successful run:

```
results/test_ms_search_hepg2/
├── prepare_fasta/         → target/decoy FASTA, ≥1 sequence
├── comet/                 → HepG2_repN.pin, ≥1 row each
├── sage/                  → results.sage.pin, ≥1 row
├── mokapot/               → *.target.psms.txt, ≥10 rows typically
├── convert_mzml/          → *.mgf, *.scans.pkl, *.index2scan.pkl
├── ms2rescore/            → {engine}.psms.tsv, ≥1 row
└── integrate_engines/     → integrated_peptides.tsv, integrated_psms.tsv
```

## Regression checks

A follow-up nf-test sprint will assert:

1. Each output file exists and is non-zero
2. Mokapot finds ≥10 target PSMs across the 3 replicates
3. INTEGRATE_ENGINES merges engine outputs without error
4. Integrated peptides table has ≥1 shared peptide across engines

Licensed tools (MSFragger, netMHCpan, GibbsCluster) get `tag:
'requires_licensed'` so GitHub Actions skips them; Monash runs them
manually.
