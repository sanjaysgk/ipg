//
// BQSR: Base Quality Score Recalibration (two passes + covariates plot)
//
// Implements steps 13–16 of the legacy cryptic peptide pipeline:
//   13. GATK4 BaseRecalibrator      first-pass recalibration table
//   14. GATK4 ApplyBQSR             apply the table to produce corrected BAM
//   15. GATK4 BaseRecalibrator      second-pass table on the corrected BAM
//   16. GATK4 AnalyzeCovariates     plot before/after for QC
//
// All three known-sites VCFs from the legacy script (dbSNP, known_indels,
// Mills) are folded into a single `ch_known_sites` list-channel; the GATK4
// BaseRecalibrator module accepts a list of known-sites paths.
//

include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_INPUT       } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_RECAL       } from '../../../modules/nf-core/samtools/index/main'
include { GATK4_BASERECALIBRATOR as BASERECALIBRATOR_FIRST  } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_BASERECALIBRATOR as BASERECALIBRATOR_SECOND } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR          } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { GATK4_ANALYZECOVARIATES  } from '../../../modules/nf-core/gatk4/analyzecovariates/main'

workflow BQSR {

    take:
    ch_input_bam        // channel: [ val(meta), path(bam) ]                    from BAM_PREP.splitncigar_bam
    ch_fasta            // channel: [ val(meta), path(fasta) ]
    ch_fai              // channel: [ val(meta), path(fai) ]
    ch_dict             // channel: [ val(meta), path(dict) ]
    ch_known_sites      // channel: [ val(meta), [path(vcf1), path(vcf2), ...] ]
    ch_known_sites_tbi  // channel: [ val(meta), [path(tbi1), path(tbi2), ...] ]

    main:

    // All modules in this subworkflow use topic-channel versions —
    // auto-collected. See align_qc for rationale.
    ch_versions = channel.empty()

    //
    // Index the SplitNCigar BAM for BaseRecalibrator
    //
    SAMTOOLS_INDEX_INPUT(ch_input_bam)

    ch_bam_bai_in = ch_input_bam
        .join(SAMTOOLS_INDEX_INPUT.out.index)
        .map { meta, bam, bai -> [meta, bam, bai, []] }   // no intervals

    //
    // Step 13: first-pass BaseRecalibrator
    //
    BASERECALIBRATOR_FIRST(
        ch_bam_bai_in,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_known_sites,
        ch_known_sites_tbi
    )

    //
    // Step 14: ApplyBQSR — module takes bare path fasta / fai / dict
    //
    ch_applybqsr_input = ch_input_bam
        .join(SAMTOOLS_INDEX_INPUT.out.index)
        .join(BASERECALIBRATOR_FIRST.out.table)
        .map { meta, bam, bai, table -> [meta, bam, bai, table, []] }

    GATK4_APPLYBQSR(
        ch_applybqsr_input,
        ch_fasta.map { _meta, fasta -> fasta },
        ch_fai.map   { _meta, fai   -> fai   },
        ch_dict.map  { _meta, dict  -> dict  }
    )

    //
    // Index the recalibrated BAM so BaseRecalibrator (second pass) can read it
    //
    SAMTOOLS_INDEX_RECAL(GATK4_APPLYBQSR.out.bam)

    //
    // Step 15: second-pass BaseRecalibrator on the corrected BAM
    //
    ch_bam_bai_recal = GATK4_APPLYBQSR.out.bam
        .join(SAMTOOLS_INDEX_RECAL.out.index)
        .map { meta, bam, bai -> [meta, bam, bai, []] }

    BASERECALIBRATOR_SECOND(
        ch_bam_bai_recal,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_known_sites,
        ch_known_sites_tbi
    )

    //
    // Step 16: AnalyzeCovariates — plot before vs after tables
    //
    // Module takes: tuple(meta, before_table, after_table, additional_table)
    // We don't have a third BQSR table, so we pass an empty list.
    //
    ch_covariates_input = BASERECALIBRATOR_FIRST.out.table
        .join(BASERECALIBRATOR_SECOND.out.table)
        .map { meta, before, after -> [meta, before, after, []] }

    GATK4_ANALYZECOVARIATES(ch_covariates_input)

    //
    // Join the recalibrated BAM with its BAI for downstream consumers
    //
    ch_recal_bam_bai = GATK4_APPLYBQSR.out.bam.join(SAMTOOLS_INDEX_RECAL.out.index)

    emit:
    recal_bam       = GATK4_APPLYBQSR.out.bam                 // [meta, bam]
    recal_bai       = SAMTOOLS_INDEX_RECAL.out.index          // [meta, bai]
    recal_bam_bai   = ch_recal_bam_bai                        // [meta, bam, bai] — input to MUTECT_CALLING
    before_table    = BASERECALIBRATOR_FIRST.out.table
    after_table     = BASERECALIBRATOR_SECOND.out.table
    covariates_pdf  = GATK4_ANALYZECOVARIATES.out.plots
    versions        = ch_versions
}
