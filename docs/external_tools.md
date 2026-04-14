# External tools

`sanjaysgk/ipg` does **not** bundle academic-licensed binaries. Users download
them separately and point the pipeline at the paths through CLI flags or a
params YAML.

## What bioconda installs for you (no user action required)

| Tool        | Used by step        | Source                          |
| ----------- | ------------------- | ------------------------------- |
| Comet       | `ms_search`         | `bioconda::comet-ms`            |
| Sage        | `ms_search`         | `bioconda::sage-proteomics`     |
| mokapot     | `ms_search`         | `bioconda::mokapot`             |
| MS2Rescore  | `ms_search`         | `bioconda::ms2rescore`          |
| pymzml      | `ms_search`         | `bioconda::pymzml`              |
| FlashLFQ    | `immunoinformatics` | `bioconda::flashlfq`            |
| BLAST+      | `immunoinformatics` | `bioconda::blast`               |

These are fetched by the singularity / conda / pixi engine profile. Nothing
to install manually.

## What you must supply yourself

| Tool            | Needed when                        | Param                  | Licence    |
| --------------- | ---------------------------------- | ---------------------- | ---------- |
| MSFragger JAR   | `msfragger` in `--ms_engines`      | `--msfragger_jar`      | Academic   |
| netMHCpan-4.1   | `--run_netmhcpan`                  | `--netmhcpan_path`     | Academic   |
| netMHCIIpan-4.3 | `--run_netmhciipan`                | `--netmhciipan_path`   | Academic   |
| GibbsCluster-2.0| `--run_gibbscluster`               | `--gibbscluster_path`  | Academic   |
| BLAST DB        | `--run_blastp_host`                | `--blast_db` (prefix)  | Public     |

## Recommended layout on Monash M3

Place everything under `/fs04/scratch2/xy86/<user>/external_tools/` so the
paths survive home-directory quota resets:

```
external_tools/
├── MSFragger-4.1/
│   └── MSFragger-4.1.jar
├── netMHCpan-4.1/
│   └── netMHCpan                      # the launcher script, not the tar
├── netMHCIIpan-4.3/
│   └── netMHCIIpan
├── gibbscluster-2.0/
│   └── GibbsCluster-2.0e_SA.pl
└── blastdb/
    ├── human.pin
    ├── human.psq
    ├── human.phr
    └── ...                            # --blast_db /path/to/blastdb/human
```

The Li/Purcell lab reference copy already lives at
`/fs04/scratch2/xy86/sanjay/ATLANTIS/RNAseq/Analysis/Cryptic/immunopeptidomics/external_tools/`
and can be pointed at directly — no need to clone.

## Example params YAML

```yaml title="params_full_run.yaml"
# MS search
ms_engines:     msfragger,comet,sage
msfragger_jar:  /fs04/scratch2/xy86/sanjay/ATLANTIS/RNAseq/Analysis/Cryptic/immunopeptidomics/external_tools/MSFragger-4.1/MSFragger-4.1.jar

# Immunoinformatics
run_netmhcpan:       true
run_netmhciipan:     true
run_gibbscluster:    true
run_flashlfq:        true
run_blastp_host:     true
hla:                 HLA-A01:01,HLA-A02:01,HLA-B07:02
netmhcpan_path:      /fs04/scratch2/xy86/.../external_tools/netMHCpan-4.1/netMHCpan
netmhciipan_path:    /fs04/scratch2/xy86/.../external_tools/netMHCIIpan-4.3/netMHCIIpan
gibbscluster_path:   /fs04/scratch2/xy86/.../external_tools/gibbscluster-2.0/GibbsCluster-2.0e_SA.pl
blast_db:            /fs04/scratch2/xy86/.../external_tools/blastdb/human
host_species:        HUMAN
```

Invoke with `-params-file params_full_run.yaml`.

## Pre-flight check

Run `bin/check_external_tools.sh <params.yaml>` before kicking off a
long job — the script resolves every external-tool path referenced in
the params file and fails fast if anything is missing, unreadable, or
not executable. Cheaper than discovering a typo 2 hours into the run.
