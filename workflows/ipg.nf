/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { samplesheetToList     } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_ipg_pipeline'

//
// Cryptic peptide discovery subworkflows (steps 1-31 of the legacy bash pipeline)
//
include { PREPARE_GENOME      } from '../subworkflows/local/prepare_genome/main'
include { ALIGN_QC            } from '../subworkflows/local/align_qc/main'
include { TRANSCRIPT_ASSEMBLY } from '../subworkflows/local/transcript_assembly/main'
include { BAM_PREP            } from '../subworkflows/local/bam_prep/main'
include { BQSR                } from '../subworkflows/local/bqsr/main'
include { MUTECT_CALLING      } from '../subworkflows/local/mutect_calling/main'
include { DB_CONSTRUCT        } from '../subworkflows/local/db_construct/main'
include { POST_MS_ANALYSIS    } from '../subworkflows/local/post_ms_analysis/main'
include { MS_SEARCH           } from '../subworkflows/local/ms_search/main'
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

    if (params.step == 'post_ms') {

        //
        // POST-MS ANALYSIS ENTRY POINT
        // Skips all upstream processing; only runs db_compare + origins.
        //
        Channel
            .fromList(samplesheetToList(params.post_ms_input, "${projectDir}/assets/schema_post_ms_input.json"))
            .map { meta, cryptic_psm, uniprot_psm, cryptic_decoy_score, uniprot_decoy_score ->
                def new_meta = meta + [
                    cryptic_decoy_score: cryptic_decoy_score,
                    uniprot_decoy_score: uniprot_decoy_score
                ]
                [ new_meta, file(cryptic_psm), file(uniprot_psm) ]
            }
            .set { ch_post_ms }

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
        Channel
            .fromList(samplesheetToList(params.ms_input, "${projectDir}/assets/schema_ms_input.json"))
            .map { meta, ms_file ->
                [ meta, file(ms_file) ]
            }
            .set { ch_ms_data }

        ch_search_fasta = channel.fromPath(params.search_fasta, checkIfExists: true).collect()

        ch_msfragger_jar = params.msfragger_jar
            ? channel.fromPath(params.msfragger_jar, checkIfExists: true).collect()
            : channel.empty()

        def engine_list = params.ms_engines.tokenize(',')

        MS_SEARCH(
            ch_ms_data,
            ch_search_fasta,
            channel.value(engine_list),
            ch_msfragger_jar
        )
        ch_versions = ch_versions.mix(MS_SEARCH.out.versions)

        //
        // Optional immunoinformatics analysis on integrated peptides.
        // Only triggers any downstream work if at least one --run_* flag is set.
        //
        if (params.run_netmhcpan || params.run_netmhciipan
                || params.run_gibbscluster || params.run_flashlfq
                || params.run_blastp_host) {
            IMMUNOINFORMATICS(MS_SEARCH.out.peptides, ch_ms_data)
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
        FASTQC(ch_samplesheet)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { row -> row[1] })
        ch_versions      = ch_versions.mix(FASTQC.out.versions.first())

        ch_reads = ch_samplesheet

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
        // SUBWORKFLOW: transcript_assembly — steps 4-5
        //
        TRANSCRIPT_ASSEMBLY(
            ALIGN_QC.out.sorted_bam,
            ch_gtf,
            ch_fasta_fai
        )
        ch_versions = ch_versions.mix(TRANSCRIPT_ASSEMBLY.out.versions)

        //
        // SUBWORKFLOW: bam_prep — steps 6-12
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
        // SUBWORKFLOW: bqsr — steps 13-16
        //
        BQSR(
            BAM_PREP.out.splitncigar_bam,
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_known_sites,
            ch_known_sites_tbi
        )
        ch_versions = ch_versions.mix(BQSR.out.versions)

        //
        // SUBWORKFLOW: mutect_calling — steps 17-23
        //
        MUTECT_CALLING(
            BQSR.out.recal_bam_bai,
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_germline_resource,
            ch_germline_resource_tbi
        )
        ch_versions = ch_versions.mix(MUTECT_CALLING.out.versions)

        //
        // SUBWORKFLOW: db_construct — steps 24-31
        //
        DB_CONSTRUCT(
            MUTECT_CALLING.out.unmasked_vcf,
            MUTECT_CALLING.out.indel_vcf,
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
