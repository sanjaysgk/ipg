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
include { PREPARE_GENOME      } from '../subworkflows/local/prepare_genome/main'
include { ALIGN_QC            } from '../subworkflows/local/align_qc/main'
include { TRANSCRIPT_ASSEMBLY } from '../subworkflows/local/transcript_assembly/main'
include { BAM_PREP            } from '../subworkflows/local/bam_prep/main'
include { BQSR                } from '../subworkflows/local/bqsr/main'
include { MUTECT_CALLING      } from '../subworkflows/local/mutect_calling/main'
include { DB_CONSTRUCT        } from '../subworkflows/local/db_construct/main'

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

    //
    // SUBWORKFLOW: Resolve, decompress, and index reference files.
    // Builds fai/dict/STAR index on-the-fly if not provided via params.
    // When all paths are explicit (e.g. from a params YAML), this is a
    // zero-cost pass-through.
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

    //
    // Build the reads-only channel for the aligners: [ meta, [r1, r2] ]
    //
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
    // StringTie needs meta.strandedness set. Merge the rseqc result back
    // into the sorted-bam meta if the user did not hard-code it in the
    // samplesheet. For Day 2 we trust the samplesheet value — strand
    // auto-detection from infer_experiment.txt is a Day 3+ enhancement.
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
    // Use the gffcompare COMBINED GTF (not the StringTie raw GTF) so that
    // downstream gff3sort -> alt_liftover -> gffread -> triple_translate
    // see TCONS_NNNN consensus transcript IDs and the associated
    // ENSG/ENST annotations from the .tracking file. This matches the
    // legacy bash pipeline at step 27 (gff3sort prefix.combined.gtf) and
    // produces FASTA headers like
    //     >TCONS_00092914_f1p15_1 ENSG00000138722.10|ENST00000264790.7:=
    // Feeding the StringTie GTF instead would produce STRG.* IDs with
    // no gene/transcript annotation in the cryptic peptide DB headers.
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
    cryptic_fasta  = DB_CONSTRUCT.out.cryptic_fasta
    versions       = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
