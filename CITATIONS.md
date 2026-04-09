# sanjaysgk/ipg: Citations

## Pipeline reference

If you use sanjaysgk/ipg for your analysis, please cite the original
methodological paper:

> Scull KE, Pandey K, Ramarathinam SH, Purcell AW.
> **Immunopeptidogenomics: harnessing RNA-seq to illuminate the dark
> immunopeptidome.** _Mol Cell Proteomics._ 2021;20:100143.
> doi: [10.1016/j.mcpro.2021.100143](https://doi.org/10.1016/j.mcpro.2021.100143)

## Pipeline framework

- [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C.
  > Nextflow enables reproducible computational workflows.
  > Nat Biotechnol. 2017 Apr 11;35(4):316-319.
  > doi: [10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820).

- [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/)

  > Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU,
  > Di Tommaso P, Nahnsen S. The nf-core framework for community-curated
  > bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278.
  > doi: [10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x).

## Alignment, transcript assembly and QC

- [STAR](https://pubmed.ncbi.nlm.nih.gov/23104886/)

  > Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P,
  > Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner.
  > Bioinformatics. 2013 Jan 1;29(1):15-21.
  > doi: [10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635).

- [SAMtools](https://pubmed.ncbi.nlm.nih.gov/19505943/)

  > Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G,
  > Durbin R, 1000 Genome Project Data Processing Subgroup. The Sequence
  > Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9.
  > doi: [10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352).

- [StringTie](https://pubmed.ncbi.nlm.nih.gov/25690850/)

  > Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT, Salzberg SL.
  > StringTie enables improved reconstruction of a transcriptome from RNA-seq
  > reads. Nat Biotechnol. 2015 Mar;33(3):290-5.
  > doi: [10.1038/nbt.3122](https://doi.org/10.1038/nbt.3122).

- [gffcompare](https://pubmed.ncbi.nlm.nih.gov/32489650/)

  > Pertea G, Pertea M. GFF Utilities: GffRead and GffCompare. F1000Res. 2020
  > Apr 28;9:ISCB Comm J-304.
  > doi: [10.12688/f1000research.23297.2](https://doi.org/10.12688/f1000research.23297.2).

- [gffread](https://pubmed.ncbi.nlm.nih.gov/32489650/)

  > Pertea G, Pertea M. GFF Utilities: GffRead and GffCompare. F1000Res. 2020
  > Apr 28;9:ISCB Comm J-304.
  > doi: [10.12688/f1000research.23297.2](https://doi.org/10.12688/f1000research.23297.2).

- [RSeQC](https://pubmed.ncbi.nlm.nih.gov/22743226/)

  > Wang L, Wang S, Li W. RSeQC: quality control of RNA-seq experiments.
  > Bioinformatics. 2012 Aug 15;28(16):2184-5.
  > doi: [10.1093/bioinformatics/bts356](https://doi.org/10.1093/bioinformatics/bts356).

- [gff3sort](https://pubmed.ncbi.nlm.nih.gov/30169744/)

  > Zhu T, Niu D-K. Frequency of intron loss correlates with processed
  > pseudogene abundance: a novel strategy to test the reverse transcriptase
  > model of intron loss. BMC Biol. 2018 Sep 1;16(1):84. (used as the GFF3
  > sorting tool described at <https://github.com/billzt/gff3sort>)

## Variant calling

- [GATK4](https://pubmed.ncbi.nlm.nih.gov/20644199/)

  > McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A,
  > Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. The Genome
  > Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA
  > sequencing data. Genome Res. 2010 Sep;20(9):1297-303.
  > doi: [10.1101/gr.107524.110](https://doi.org/10.1101/gr.107524.110).

- [Mutect2 / Tumour-only somatic calling](https://www.biorxiv.org/content/10.1101/861054v1)

  > Benjamin D, Sato T, Cibulskis K, Getz G, Stewart C, Lichtenstein L.
  > Calling Somatic SNVs and Indels with Mutect2. bioRxiv 861054.
  > doi: [10.1101/861054](https://doi.org/10.1101/861054).

## Cryptic peptide / IPG custom tools

The five custom C tools (`curate_vcf`, `revert_headers`, `alt_liftover`,
`triple_translate`, `squish`) are taken from the
[immunopeptidogenomics](https://github.com/kescull/immunopeptidogenomics)
repository, originally developed by Kate Scull et al. (Purcell Lab,
Monash University) for the cryptic peptide discovery method described in:

> Scull KE, Pandey K, Ramarathinam SH, Purcell AW. Immunopeptidogenomics:
> harnessing RNA-seq to illuminate the dark immunopeptidome. _Mol Cell
> Proteomics._ 2021;20:100143.
> doi: [10.1016/j.mcpro.2021.100143](https://doi.org/10.1016/j.mcpro.2021.100143)

The pinned source for the binaries used by this pipeline lives at
<https://github.com/sanjaysgk/immunopeptidogenomics> (a fork of the
upstream `kescull/immunopeptidogenomics`) and is built into a
reproducible OCI image at `ghcr.io/sanjaysgk/ipg-tools` via
`containers/ipg-tools/Dockerfile` in this repository.

## QC reporting

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

  > Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput
  > Sequence Data [Online].

- [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/)

  > Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis
  > results for multiple tools and samples in a single report. Bioinformatics.
  > 2016 Oct 1;32(19):3047-8.
  > doi: [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354).

## Software packaging / containerisation tools

- [Pixi](https://pixi.sh) — reproducible per-project conda environment used
  for the developer toolchain (`nextflow`, `nf-core`, `nf-test`, linters)
  and locally for non-container pipeline runs.

- [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)

  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH,
  > Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and
  > comprehensive software distribution for the life sciences. Nat Methods.
  > 2018 Jul;15(7):475-476.
  > doi: [10.1038/s41592-018-0046-7](https://doi.org/10.1038/s41592-018-0046-7).

- [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/)

  > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J,
  > Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC,
  > Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI,
  > Perez-Riverol Y. BioContainers: an open-source and community-driven
  > framework for software standardization. Bioinformatics. 2017 Aug
  > 15;33(16):2580-2582.
  > doi: [10.1093/bioinformatics/btx192](https://doi.org/10.1093/bioinformatics/btx192).

- [Singularity / Apptainer](https://pubmed.ncbi.nlm.nih.gov/28494014/)

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for
  > mobility of compute. PLoS One. 2017 May 11;12(5):e0177459.
  > doi: [10.1371/journal.pone.0177459](https://doi.org/10.1371/journal.pone.0177459).
