# post_ms test fixtures

Synthetic minimum-viable dataset for exercising the `--step post_ms`
subworkflow (DB_COMPARE Phase 1 → ORIGINS simple → DB_COMPARE Phase 2 →
ORIGINS full). All files are hand-built so the outputs are
deterministic and guaranteed non-empty.

## Files

| File | Purpose |
|---|---|
| `cryptic.peptides.csv` | 10 peptides (5 cryptic-only, 3 shared with UniProt, 2 below threshold) |
| `uniprot.peptides.csv` | 12 peptides (7 UniProt-only, 3 shared, 2 below threshold) |
| `samplesheet.csv` | 1-row post_ms samplesheet, thresholds at 20.0 `-10lgP` |
| `uniprot.fasta` | 10-protein UniProt mini-DB covering the UniProt peptide Accessions |
| `transcriptome.fa` | 8-transcript cryptic mini-transcriptome covering TCONS IDs |
| `prefix.tracking` | gffcompare-style tracking for the 8 TCONS transcripts |

## Expected post_ms outputs

At `--cryptic_decoy_score 20.0 --uniprot_decoy_score 20.0`:

- **cryptic.peptides.csv** filtered: 8 rows (`LOWSCORE_CRYPTIC_A/B` drop)
- **uniprot.peptides.csv** filtered: 10 rows (`LOWSCORE_NOISE_A/B` drop)
- **Phase 1 cryptic-only**: 5 peptides
  - `CRYPTPEPTIDEONE`, `CRYPTPEPTIDETWO`, `CRYPTPEPTIDETHREE`,
    `CRYPTPEPTIDEFOUR`, `CRYPTPEPTIDEFIVE`
- **Shared (discarded from cryptic-only)**: 3 peptides
- **Origins output**: 5 TCONS IDs mapped to XLOC genes

## Launch

```bash
cd $IPG_REPO
pixi run nextflow run . \
    -profile pixi,monash \
    --step post_ms \
    --post_ms_input       tests/data/post_ms/samplesheet.csv \
    --uniprot_fasta       tests/data/post_ms/uniprot.fasta \
    --transcriptome_fasta tests/data/post_ms/transcriptome.fa \
    --prefix_tracking     tests/data/post_ms/prefix.tracking \
    --outdir              test_post_ms_results
```

Runs in ~30 seconds on any laptop; no licensed tools needed.

## Assertions (for nf-test later)

- DB_COMPARE_PHASE1 emits `*_cryptic_only.txt` with ≥1 row
- DB_COMPARE_PHASE2 emits `*_unambiguous_unconventional.txt` with ≥0 rows
- ORIGINS_SIMPLE emits ≥1 origin annotation
- ORIGINS_FULL emits ≥1 origin annotation
- No empty outputs
