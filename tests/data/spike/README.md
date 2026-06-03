# Synthetic spike-in positive control (`--step ms_search`)

A **known-answer** test for the MS-search half of the pipeline. A single cryptic
peptide, **`CRYPTICALLY`**, is planted in a tiny FASTA and given matching theoretical
spectra. A green run recovers it in the integrated peptide table at ≤1% FDR; a red run
fails fast and localises the broken link
(`PREPARE_FASTA → COMET/SAGE → MOKAPOT → MS2RESCORE → INTEGRATE_ENGINES`).

Unlike the chr22 / HepG2 fixtures this needs **no download** — everything is generated
deterministically (RNG seed 42), so the committed `.mzML`/`.fasta` are regenerable
byte-for-byte.

## Source

Synthetic. `make_synth_spectra.py` emits, for the planted peptide plus ~300 background
peptides, singly-charged theoretical b/y-ion MS2 spectra (mono-isotopic), plus a set of
partial reversed-sequence "decoy-seed" spectra so Mokapot/FDR has a decoy distribution to
fit. Background peptides are ≥10 aa and the marker is an 11-mer so every spectrum clears
Sage's `min_peaks:15` after the `fragment_min_mz:200` cut.

## Contents (~1.3 MB)

```
spike/
├── make_synth_spectra.py   # generator (seed 42, reproducible)
├── synth.mzML              # 451 spectra (301 target + 150 decoy-seed)
├── synth_db.fasta          # 301 sequences incl. >CRYPTIC_SPIKE_CRYPTICALLY
├── samplesheet_ms.csv      # one SPIKE row (absolute ms_file path — regenerate per machine)
└── README.md
```

## Regenerate (per machine)

nf-schema validates the samplesheet's `ms_file` against the launch dir, so the path is
written **absolute** — regenerate after cloning so it points at your checkout:

```bash
pixi run -e ms2rescore python tests/data/spike/make_synth_spectra.py \
    --out-mzml tests/data/spike/synth.mzML \
    --out-fasta tests/data/spike/synth_db.fasta
```

## Run

```bash
pixi run nextflow run . -profile pixi,test_spike --outdir results_spike    # ~5 min
```

Add MSFragger when a JAR is available (its academic-license JAR cannot be committed):

```bash
pixi run nextflow run . -profile pixi,test_spike \
    --ms_engines comet,sage,msfragger \
    --msfragger_jar /path/to/MSFragger-4.2.jar --outdir results_spike
```

## Expected outputs (validated 2026-06-03, Comet 2026.01 + Sage 0.14.7 + MSFragger 4.2)

`results_spike/integrate/integrated_peptides.tsv` contains the planted peptide:

```
peptide        engine                 peptide_qval
CRYPTICALLY    ['comet', 'sage']      [0.0033, 0.0033]      # +'msfragger' when its JAR is supplied
```

Assertion: `CRYPTICALLY` present with `peptide_qval ≤ 0.01`. Per-engine Mokapot tables under
`results_spike/mokapot/` also carry it. `chimeric_PSMs.txt` / `chimera_only_peptides.txt` are
empty (single unambiguous assignment for the planted scan).
