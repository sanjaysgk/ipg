# Synthetic spike-in positive control (full pipeline)

A **known-answer** test: a single cryptic peptide, **`CRYPTICALLY`**, is planted so it must be
recovered through the pipeline. Green = every link works; red = fails fast and localises the
break. Two halves, runnable separately or chained end-to-end:

- **db_construct half** — a tiny synthetic genome whose transcript encodes `CRYPTICALLY` in a
  non-canonical frame; `db_construct` must surface it in the cryptic FASTA
  (`STAR → StringTie → gffread → triple_translate → squish`).
- **ms_search half** — theoretical b/y spectra for `CRYPTICALLY` + background; `ms_search` must
  recover it (`PREPARE_FASTA → COMET/SAGE → MOKAPOT → MS2RESCORE → INTEGRATE_ENGINES`).

No download — everything is generated deterministically (RNG seed 42); committed fixtures are
regenerable byte-for-byte.

## Contents (~1.4 MB)

```
spike/
├── make_synth_spectra.py     # MS-half generator (mzML + FASTA + samplesheet)
├── make_synth_genome.py      # db_construct-half generator (genome + reads + VCFs)
├── synth.mzML                # 451 spectra (301 target + 150 decoy-seed)
├── synth_db.fasta            # 301 seqs incl. >CRYPTIC_SPIKE_CRYPTICALLY (MS-half search DB)
├── samplesheet_ms.csv        # ms_search input (absolute ms_file path — regenerate per machine)
├── genome/                   # db_construct fixtures
│   ├── genome.fa             #   3 kb contig: gene1 carries CRYPTICALLY, gene2 is a sentinel*
│   ├── genes.gtf, genes.bed
│   ├── reads_{R1,R2}.fastq.gz
│   ├── {dbsnp,known_indels,mills}.vcf.gz(.tbi)   # BQSR known-sites (empty)
│   ├── germline_resource.vcf.gz(.tbi)            # AF sites for GATK GetPileupSummaries
│   └── samplesheet_rnaseq.csv
└── README.md
```
\* gene2 sentinel: `triple_translate` silently dropped the *last* transcript (kescull bug, fixed
in the fork — see ipg issue #43); a 2nd transcript keeps `CRYPTICALLY` from being last. With the
fixed tool a single transcript also works, but ≥2 is more realistic.

## Regenerate (per machine — paths are written absolute)

```bash
pixi run -e ms2rescore python tests/data/spike/make_synth_spectra.py \
    --out-mzml tests/data/spike/synth.mzML --out-fasta tests/data/spike/synth_db.fasta
python3 tests/data/spike/make_synth_genome.py --out-dir tests/data/spike/genome
# then bgzip + tabix the 4 genome VCFs:
for v in dbsnp known_indels mills germline_resource; do
  pixi run bgzip -f tests/data/spike/genome/$v.vcf && pixi run tabix -p vcf -f tests/data/spike/genome/$v.vcf.gz
done
```

## Run

```bash
# MS-half only (~5 min) — searches the hand-made synth_db.fasta
pixi run nextflow run . -profile pixi,test_spike --outdir results_spike

# db_construct half only (~3.5 min) — builds the cryptic FASTA from synthetic RNAseq
pixi run nextflow run . -profile pixi,test_spike --step db_construct --outdir results_db
```

### Full two-phase e2e (the real known-answer chain)

Phase 2 searches db_construct's **own** cryptic FASTA, concatenated with a background proteome
(a cryptic-only DB is too small for Mokapot's FDR — real runs concatenate with UniProt). The
background here is `synth_db.fasta` **minus** `CRYPTICALLY`, so a hit proves the whole chain.

```bash
SP=tests/data/spike
pixi run nextflow run . -profile pixi,test_spike --step db_construct --outdir $SP/results_e2e/db
awk '/^>/{k=($0!~/CRYPTIC_SPIKE/)} k' $SP/synth_db.fasta > $SP/results_e2e/background.fasta
cat $SP/results_e2e/db/squish/SPIKE_cryptic.fasta $SP/results_e2e/background.fasta > $SP/results_e2e/search.fasta
pixi run nextflow run . -profile pixi,test_spike --step ms_search \
    --search_fasta $SP/results_e2e/search.fasta --outdir $SP/results_e2e/ms
```

## Expected outputs (validated 2026-06-03; Comet 2026.01 · Sage 0.14.7 · MSFragger 4.2)

- **db_construct:** `results_db/squish/SPIKE_cryptic.fasta` contains `CRYPTICALLY` (within a
  squish-concatenated ORF, e.g. `>TCONS_..._f3p1_1 ... WWWWWCRYPTICALLYWWWWW...`).
- **ms_search / e2e:** `…/integrate/integrated_peptides.tsv` contains `CRYPTICALLY`,
  `engine=['comet','sage']`, `peptide_qval ≤ 0.01` (q≈0.0033), chimeric files empty.

Assertion for both: `CRYPTICALLY` present at `peptide_qval ≤ 0.01`. MSFragger is omitted by
default (academic-license JAR); add with `--ms_engines comet,sage,msfragger --msfragger_jar …`.
