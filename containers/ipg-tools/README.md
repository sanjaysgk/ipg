# ipg-tools container

Reproducible build of the five C tools from the `immunopeptidogenomics` repository
(originally written by Scull et al. for their 2020 *Mol Cell Proteomics* paper on
cryptic peptide discovery from RNA-seq).

| Tool | Used in pipeline step |
|---|---|
| `curate_vcf` | 23 — VCF curation before `FastaAlternateReferenceMaker` |
| `revert_headers` | 26 — restore reference chromosome headers on alt FASTA |
| `alt_liftover` | 28 — lift GTF coordinates to alternate reference |
| `triple_translate` | 30 — 3-frame translation with transcript tracking |
| `squish` | 31 — multi-FASTA deduplication into final cryptic peptide database |

## Source pin

Built from `https://github.com/sanjaysgk/immunopeptidogenomics.git` at commit
`fef8e68ab86dac49070542fe40b477ec76058aaa` — a fork of
[kescull/immunopeptidogenomics](https://github.com/kescull/immunopeptidogenomics)
under `sanjaysgk` control so the source cannot disappear from under us.

The pin lives in two places:

- `containers/ipg-tools/Dockerfile` — `ARG IPG_REPO_URL` and `ARG IPG_REPO_SHA`
- `.github/workflows/build-ipg-tools.yml` — `env.IPG_REPO_SHA`

Both must be updated together when bumping the version.

## Building locally

```bash
docker build \
    -t ipg-tools:local \
    -f containers/ipg-tools/Dockerfile \
    containers/ipg-tools
```

To override the source repo or commit at build time (e.g. for testing a fix):

```bash
docker build \
    --build-arg IPG_REPO_URL=https://github.com/<you>/immunopeptidogenomics.git \
    --build-arg IPG_REPO_SHA=<sha> \
    -t ipg-tools:dev \
    -f containers/ipg-tools/Dockerfile \
    containers/ipg-tools
```

## Published image

CI publishes to `ghcr.io/sanjaysgk/ipg-tools:<tag>` on every push to `main`
and on every git tag matching `ipg-tools-v*`. The Nextflow modules under
`modules/local/{curate_vcf,revert_headers,alt_liftover,triple_translate,squish}`
reference this image.

## License

The C sources are (C) their respective authors under the MIT license (see the
upstream `immunopeptidogenomics` repository). The Dockerfile and build workflow
in this directory are part of the `sanjaysgk/ipg` pipeline and covered by that
pipeline's license.
