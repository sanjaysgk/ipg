/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_ipg_pipeline'

//
// Cryptic peptide discovery subworkflows (steps 1-31 of the legacy bash pipeline)
//
include { SAMPLESHEET_TO_CHANNEL } from '../subworkflows/local/samplesheet_to_channel/main'
include { PREPARE_GENOME      } from '../subworkflows/local/prepare_genome/main'
include { ALIGN_QC            } from '../subworkflows/local/align_qc/main'
include { TRANSCRIPT_ASSEMBLY } from '../subworkflows/local/transcript_assembly/main'
include { BAM_PREP            } from '../subworkflows/local/bam_prep/main'
include { BQSR                } from '../subworkflows/local/bqsr/main'
include { BAM_VARIANT_CALLING_MUTECT2 } from '../subworkflows/local/bam_variant_calling_mutect2/main'
include { VCF_CURATE          } from '../subworkflows/local/vcf_curate/main'
include { VCF_ANNOTATE_ALL    } from '../subworkflows/local/vcf_annotate_all/main'
include { DB_CONSTRUCT        } from '../subworkflows/local/db_construct/main'
include { POST_MS_ANALYSIS    } from '../subworkflows/local/post_ms_analysis/main'
include { MS_SEARCH           } from '../subworkflows/local/ms_search/main'
include { DENOVO              } from '../subworkflows/local/denovo/main'
include { DOWNLOAD_UNIPROT_PROTEOME } from '../modules/local/download_uniprot_proteome/main'
include { VALIDATE_CRYPTIC    } from '../subworkflows/local/validate_cryptic/main'
include { ANNOTATE_ORIGIN     } from '../modules/local/annotate_origin/main'
include { ORIGINS as ORIGINS_ENSEMBL } from '../modules/local/origins/main'
include { MERGE_ORIGIN        } from '../modules/local/merge_origin/main'
include { CRYPTIC_REPORT      } from '../modules/local/cryptic_report/main'
include { IMMUNOINFORMATICS   } from '../subworkflows/local/immunoinformatics/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow IPG {

    take:
    ch_samplesheet      // channel: [ val(meta), [r1, r2] ]   samplesheet rows

    main:

    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

    // Placeholder for emit — may or may not be produced depending on step
    ch_cryptic_fasta = channel.empty()

    //
    // Single samplesheet parser: validates --step / --*_input combo and
    // emits typed channels (rnaseq / ms / post_ms) — only the one matching
    // params.step is populated.
    //
    SAMPLESHEET_TO_CHANNEL(
        params.step,
        params.input,
        params.ms_input,
        params.post_ms_input,
        "${projectDir}/assets/schema_input.json",
        "${projectDir}/assets/schema_ms_input.json",
        "${projectDir}/assets/schema_post_ms_input.json"
    )
    ch_rnaseq_input  = SAMPLESHEET_TO_CHANNEL.out.rnaseq
    ch_ms_input      = SAMPLESHEET_TO_CHANNEL.out.ms
    ch_post_ms_input = SAMPLESHEET_TO_CHANNEL.out.post_ms

    if (params.step == 'post_ms') {

        //
        // POST-MS ANALYSIS ENTRY POINT (deprecated)
        // Skips all upstream processing; only runs db_compare + origins.
        // PEAKS is now a first-class engine: prefer `--step ms_search
        // --ms_engines peaks` so PEAKS results flow through MS2Rescore and the
        // class-separated FDR like every other engine. db_compare/post_ms is
        // kept for one release to reproduce the legacy PEAKS-only path.
        //
        log.warn(
            "`--step post_ms` (db_compare + origins) is deprecated and will be removed in a future release. " +
            "Run PEAKS as a search engine instead: `--step ms_search --ms_engines peaks --peaks_psm_csv <db.psms.csv>`."
        )

        ch_post_ms = ch_post_ms_input

        ch_tracking            = channel.fromPath(params.prefix_tracking, checkIfExists: true).collect()
        ch_uniprot_fasta       = channel.fromPath(params.uniprot_fasta, checkIfExists: true).collect()
        ch_transcriptome_fasta = channel.fromPath(params.transcriptome_fasta, checkIfExists: true).collect()

        POST_MS_ANALYSIS(
            ch_post_ms,
            ch_tracking,
            ch_uniprot_fasta,
            ch_transcriptome_fasta
        )
        ch_versions = ch_versions.mix(POST_MS_ANALYSIS.out.versions)

    } else if (params.step == 'ms_search') {

        //
        // MS SEARCH ENTRY POINT
        // Runs search engines on MS data against a FASTA database.
        //
        // Per-sample DB routing. Each MS file declares its database via the `db`
        // column (e.g. results/db_construct/<sample>/<sample>_cryptic.fasta), or
        // falls back to a single --search_fasta for every sample. A sample's MS
        // files (fractions + replicates) are pooled into ONE search against its
        // ONE db; per-rep stays recoverable via the `run` column downstream.
        def fallback_db = params.search_fasta ? file(params.search_fasta, checkIfExists: true) : null

        ch_ms_db = ch_ms_input
            .map { meta, ms_file ->
                def raw = meta.db ?: fallback_db
                if (!raw) {
                    error("ms_search: MS sample '${meta.id}' has no 'db' column and --search_fasta is not set; " +
                        "cannot route this file to a database. Add a 'db' column per row, or pass --search_fasta.")
                }
                [ meta.id, meta, ms_file, file(raw, checkIfExists: true) ]
            }
            .groupTuple(by: 0)
            .map { id, metas, files, dbs ->
                def udb = dbs.collect { it.toString() }.unique()
                if (udb.size() > 1) {
                    error("ms_search: sample '${id}' rows name different databases ${udb}; all MS files of one " +
                        "sample must share a single db (search is keyed on the 'sample' column).")
                }
                def conds = metas.collect { it.condition }.findAll { it }.unique()
                def base  = metas[0] + [ id: id, db: dbs[0].toString(), num_ms_files: files.size() ]
                if (conds.size() == 1) { base = base + [ condition: conds[0] ] }
                [ base, files, dbs[0] ]
            }

        // Log the MS->DB plan once, before any search.
        ch_ms_db
            .map { meta, files, db -> "    ${meta.id}  <- ${files.size()} MS file(s)  vs  ${db.name}" }
            .collect()
            .subscribe { plan ->
                log.info("ms_search DB mapping — ${plan.size()} sample(s) -> ${plan.size()} search(es):\n${plan.sort().join('\n')}")
            }

        // Canonical proteome. Unused in 'cryptic_only'. A local --canonical_fasta
        // wins; otherwise DOWNLOAD_UNIPROT_PROTEOME fetches canonical_proteome_url
        // (default UP000005640) and caches it via storeDir.
        if (params.db_search_mode == 'cryptic_only') {
            ch_canonical_fasta = channel.value([])
        } else if (params.canonical_fasta) {
            ch_canonical_fasta = channel.fromPath(params.canonical_fasta, checkIfExists: true).collect()
        } else {
            DOWNLOAD_UNIPROT_PROTEOME(channel.value(params.canonical_proteome_url))
            ch_versions        = ch_versions.mix(DOWNLOAD_UNIPROT_PROTEOME.out.versions)
            ch_canonical_fasta = DOWNLOAD_UNIPROT_PROTEOME.out.fasta.collect()
        }

        // When --msfragger_jar isn't provided, pass the NO_FILE sentinel
        // so the MSFRAGGER process still schedules (an empty channel
        // would starve the input and the task never runs). The module
        // detects size==0 and falls back to the bioconda wrapper.
        ch_msfragger_jar = params.msfragger_jar
            ? channel.fromPath(params.msfragger_jar, checkIfExists: true).collect()
            : channel.value(file("${projectDir}/assets/NO_FILE"))

        def engine_list = params.ms_engines.tokenize(',')

        MS_SEARCH(
            ch_ms_db,
            ch_canonical_fasta,
            channel.value(engine_list),
            ch_msfragger_jar
        )
        ch_versions = ch_versions.mix(MS_SEARCH.out.versions)

        // De novo discovery lane (Phase 1a) — parallel to the DB-search engines, same
        // spectra + per-sample cryptic DB. One de novo task per spectra file (unique
        // meta.id) so predictions join 1:1 to their DB for classification.
        if (params.run_denovo) {
            ch_denovo_in = ch_ms_db
                .transpose()
                .map { meta, ms_file, db -> [meta + [id: "${meta.id}_${ms_file.baseName}"], ms_file, db] }
            DENOVO(ch_denovo_in, ch_canonical_fasta)
            ch_versions = ch_versions.mix(DENOVO.out.versions)
        }

        //
        // Optional PepQuery2 spectrum-level validation of the cryptic subset.
        // Annotates pepquery_status onto the integrated peptides before any
        // downstream analysis consumes them. Each sample validates against its
        // own search db (joined by meta), matching the per-sample search.
        //
        ch_ms_peptides = MS_SEARCH.out.peptides
        if (params.run_pepquery) {
            VALIDATE_CRYPTIC(MS_SEARCH.out.peptides, MS_SEARCH.out.mgf, MS_SEARCH.out.search_db)
            ch_versions    = ch_versions.mix(VALIDATE_CRYPTIC.out.versions)
            ch_ms_peptides = VALIDATE_CRYPTIC.out.peptides
        }

        //
        // Optional cryptic-origin annotation: backtrack each cryptic peptide to
        // its source ORF/transcript/gene via the 3-frame translation + gffcompare
        // .tracking (both from DB_CONSTRUCT). No network — see ANNOTATE_ORIGIN.
        //
        if (params.run_annotate_origin) {
            if (!params.origin_orf_fasta) error("--run_annotate_origin requires --origin_orf_fasta (3-frame ORF FASTA from db_construct)")
            if (!params.origin_gtf)       error("--run_annotate_origin requires --origin_gtf (gffcompare .combined.gtf from db_construct)")
            if (!params.origin_tmap)      error("--run_annotate_origin requires --origin_tmap (gffcompare .tmap from db_construct)")
            ANNOTATE_ORIGIN(
                ch_ms_peptides,
                file(params.origin_orf_fasta, checkIfExists: true),
                file(params.origin_gtf,       checkIfExists: true),
                file(params.origin_tmap,      checkIfExists: true)
            )
            ch_versions    = ch_versions.mix(ANNOTATE_ORIGIN.out.versions)
            ch_ms_peptides = ANNOTATE_ORIGIN.out.peptides

            //
            // Rich region/frame/biotype origin calls (5'/3'UTR, in/out-of-frame,
            // biotype, ENSG, Ensembl coords) via the origins Ensembl pass — the
            // default. Set --annotate_with_ensembl false for a network-free run.
            //
            if (params.annotate_with_ensembl) {
                def missing = []
                if (!params.origin_transcriptome) missing << '--origin_transcriptome (reference transcriptome FASTA)'
                if (!params.origin_tracking)      missing << '--origin_tracking (gffcompare .tracking)'
                if (!params.canonical_fasta)      missing << '--canonical_fasta (UniProt proteome)'
                if (missing) error("--annotate_with_ensembl (default true) needs ${missing.join(', ')} for the rich origins region calls. Provide them, or pass --annotate_with_ensembl false for local-only origin annotation.")
                ORIGINS_ENSEMBL(
                    ANNOTATE_ORIGIN.out.cryptic_list.filter { _meta, f -> f.size() > 0 },
                    file(params.origin_tracking,      checkIfExists: true),
                    file(params.canonical_fasta,      checkIfExists: true),
                    file(params.origin_transcriptome, checkIfExists: true)
                )
                ch_versions = ch_versions.mix(ORIGINS_ENSEMBL.out.versions)

                //
                // Graft the rich region/frame/biotype/ENSG calls onto the peptide
                // table. origins emits prot + rna CSVs under origins_csv — pick the
                // rna one, join to the table, and merge by sequence.
                //
                def ch_origins_rna = ORIGINS_ENSEMBL.out.origins_csv
                    .map { meta, f -> [meta, (f instanceof List ? f : [f]).find { it.name.endsWith('_rna.csv') }] }
                    .filter { _meta, f -> f != null }
                MERGE_ORIGIN(ANNOTATE_ORIGIN.out.peptides.join(ch_origins_rna))
                ch_versions = ch_versions.mix(MERGE_ORIGIN.out.versions)

                // Samples with no cryptic peptides never ran origins; keep their
                // annotate_origin table (remainder join, as in VALIDATE_CRYPTIC).
                def ch_no_region = ANNOTATE_ORIGIN.out.peptides
                    .join(MERGE_ORIGIN.out.peptides, remainder: true)
                    .filter { it -> it.size() > 2 && it[2] == null }
                    .map { meta, orig, _merged -> [meta, orig] }
                ch_ms_peptides = MERGE_ORIGIN.out.peptides.mix(ch_no_region)
            }
        }

        //
        // Optional cryptic-discovery HTML report: class breakdown, cross-engine
        // corroboration, PepQuery/spectral confidence, novelty class + origin,
        // before/after-rescore shift by class.
        //
        if (params.run_cryptic_report) {
            ch_report_in = ch_ms_peptides.join(
                MS_SEARCH.out.rescored_per_eng
                    .map { meta, _eng, f -> [meta, f] }
                    .groupTuple(by: 0)
            )   // [meta, peptides_tsv, [rescored tsvs]]
            CRYPTIC_REPORT(ch_report_in)
            ch_versions = ch_versions.mix(CRYPTIC_REPORT.out.versions)
        }

        //
        // Optional immunoinformatics analysis on integrated peptides.
        // Only triggers any downstream work if at least one --run_* flag is set.
        //
        if (params.run_netmhcpan || params.run_netmhciipan
                || params.run_gibbscluster || params.run_flashlfq
                || params.run_blastp_host) {
            IMMUNOINFORMATICS(ch_ms_peptides, ch_ms_data)
            ch_versions = ch_versions.mix(IMMUNOINFORMATICS.out.versions)
        }

    } else {

        //
        // DB_CONSTRUCT ENTRY POINT (default)
        // Full pipeline from FASTQ → cryptic peptide database.
        //

        //
        // SUBWORKFLOW: Resolve, decompress, and index reference files.
        //
        PREPARE_GENOME()
        ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

        ch_fasta                 = PREPARE_GENOME.out.fasta
        ch_fai                   = PREPARE_GENOME.out.fasta_fai
        ch_dict                  = PREPARE_GENOME.out.fasta_dict
        ch_fasta_fai             = PREPARE_GENOME.out.fasta_fai_combo
        ch_star_index            = PREPARE_GENOME.out.star_index
        ch_gtf                   = PREPARE_GENOME.out.gtf
        ch_rseqc_bed             = PREPARE_GENOME.out.rseqc_bed
        ch_known_sites           = PREPARE_GENOME.out.known_sites
        ch_known_sites_tbi       = PREPARE_GENOME.out.known_sites_tbi
        ch_germline_resource     = PREPARE_GENOME.out.germline_resource
        ch_germline_resource_tbi = PREPARE_GENOME.out.germline_resource_tbi

        //
        // MODULE: FastQC — per-sample raw read QC
        //
        FASTQC(ch_rnaseq_input)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { row -> row[1] })
        ch_versions      = ch_versions.mix(FASTQC.out.versions.first())

        ch_reads = ch_rnaseq_input

        //
        // SUBWORKFLOW: align_qc — steps 1-3
        //
        ALIGN_QC(
            ch_reads,
            ch_star_index,
            ch_gtf,
            ch_fasta_fai,
            ch_rseqc_bed
        )
        ch_versions = ch_versions.mix(ALIGN_QC.out.versions)

        //
        // SUBWORKFLOW: bam_prep — steps 6-12. STAR-aligned reps are merged at
        // MarkDuplicates; the merged BAM feeds both transcript assembly and
        // variant calling.
        //
        BAM_PREP(
            ALIGN_QC.out.star_bam_unsorted,
            ch_reads,
            ch_fasta,
            ch_fai,
            ch_dict
        )
        ch_versions = ch_versions.mix(BAM_PREP.out.versions)

        //
        // Resolve strandedness from RSeQC infer_experiment when the samplesheet
        // value is 'auto' (explicit forward/reverse/unstranded are honoured, with
        // a warning on disagreement). Reps share a library, so reduce to one call
        // per sample, then re-stamp BOTH merged BAMs — the one feeding StringTie
        // and the one feeding the variant branch — so they carry identical meta:
        // db_construct keys its cross-branch joins on the full meta.
        //
        def strand_thr = (params.strandedness_threshold ?: 0.8) as double
        ch_sample_strand = ALIGN_QC.out.strandedness_txt
            .map { meta, txt ->
                def fr = (txt.text =~ /Fraction of reads explained by "[^"]*":\s*([0-9.]+)/).collect { (it[1] as double) }
                def fwd = fr ? fr[0] : 0.0
                def rev = fr.size() > 1 ? fr[1] : 0.0
                def inferred = fwd >= strand_thr ? 'forward' : (rev >= strand_thr ? 'reverse' : 'unstranded')
                def declared = meta.strandedness
                if (declared == 'auto')
                    log.info("strandedness[${meta.id}]: auto -> ${inferred} (forward=${fwd}, reverse=${rev}, threshold=${strand_thr})")
                else if (declared != inferred)
                    log.warn("strandedness[${meta.id}]: declared '${declared}' disagrees with inferred '${inferred}' (forward=${fwd}, reverse=${rev})")
                [ (meta.sample ?: meta.id), (declared == 'auto' ? inferred : declared) ]
            }
            .groupTuple()
            .map { sample, vals -> [ sample, vals[0] ] }

        ch_markdup_strand = BAM_PREP.out.markdup_bam
            .map { meta, bam -> [ meta.id, meta, bam ] }
            .join(ch_sample_strand)
            .map { _id, meta, bam, s -> [ meta + [strandedness: s], bam ] }

        ch_splitncigar_strand = BAM_PREP.out.splitncigar_bam
            .map { meta, bam -> [ meta.id, meta, bam ] }
            .join(ch_sample_strand)
            .map { _id, meta, bam, s -> [ meta + [strandedness: s], bam ] }

        //
        // SUBWORKFLOW: transcript_assembly — steps 4-5, on the merged
        // (coord-sorted, dup-marked) MarkDuplicates BAM.
        //
        TRANSCRIPT_ASSEMBLY(
            ch_markdup_strand,
            ch_gtf,
            ch_fasta_fai
        )
        ch_versions = ch_versions.mix(TRANSCRIPT_ASSEMBLY.out.versions)

        //
        // SUBWORKFLOW: bqsr — steps 13-16
        //
        BQSR(
            ch_splitncigar_strand,
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_known_sites,
            ch_known_sites_tbi
        )
        ch_versions = ch_versions.mix(BQSR.out.versions)

        //
        // SUBWORKFLOW: bam_variant_calling_mutect2 — steps 17-22
        //
        BAM_VARIANT_CALLING_MUTECT2(
            BQSR.out.recal_bam_bai,
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_germline_resource,
            ch_germline_resource_tbi
        )
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING_MUTECT2.out.versions)

        //
        // SUBWORKFLOW: vcf_curate — step 23 (curate_vcf -> unmasked + indel VCFs)
        //
        VCF_CURATE(BAM_VARIANT_CALLING_MUTECT2.out.vcf)
        ch_versions = ch_versions.mix(VCF_CURATE.out.versions)

        //
        // SUBWORKFLOW: vcf_annotate_all — leaf branch, published artifacts only
        //
        if (params.tools) {
            VCF_ANNOTATE_ALL(
                BAM_VARIANT_CALLING_MUTECT2.out.vcf.join(BAM_VARIANT_CALLING_MUTECT2.out.tbi),
                ch_fasta,
                ch_fai,
                params.vep_cache         ? file(params.vep_cache)         : [],
                params.vep_genome,
                params.vep_species,
                params.vep_cache_version,
                params.snpeff_cache      ? file(params.snpeff_cache)      : [],
                params.snpeff_db,
                params.alphamissense_tsv ? file(params.alphamissense_tsv) : [],
                params.tools
            )
            ch_versions = ch_versions.mix(VCF_ANNOTATE_ALL.out.versions)
        }

        //
        // SUBWORKFLOW: db_construct — steps 24-31
        //
        DB_CONSTRUCT(
            VCF_CURATE.out.unmasked_vcf,
            VCF_CURATE.out.indel_vcf,
            TRANSCRIPT_ASSEMBLY.out.combined_gtf,
            TRANSCRIPT_ASSEMBLY.out.tracking,
            ch_fasta,
            ch_fai,
            ch_dict
        )
        ch_versions = ch_versions.mix(DB_CONSTRUCT.out.versions)
        ch_cryptic_fasta = DB_CONSTRUCT.out.cryptic_fasta

    } // end step branching

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'sanjaysgk_ipg_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config
        ? channel.fromPath(params.multiqc_config, checkIfExists: true)
        : channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : channel.empty()

    summary_params      = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList()
    cryptic_fasta  = ch_cryptic_fasta
    versions       = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
