//
// MUTECT_CALLING: Mutect2 tumour-only somatic calling + contamination + curation
//
// Implements steps 17–23 of the legacy cryptic peptide pipeline:
//   17. GATK4 Mutect2                 tumour-only (no matched normal)
//   18. GATK4 LearnReadOrientationModel
//   19. GATK4 GetPileupSummaries
//   20. GATK4 CalculateContamination
//   21. GATK4 FilterMutectCalls
//   22. GATK4 SelectVariants           --exclude-filtered --sites-only-vcf-output
//   23. curate_vcf                     IPG local module — produces *_unmasked.vcf
//                                                         and *_indel.vcf for the
//                                                         FastaAlternateReferenceMaker
//                                                         step in DB_CONSTRUCT
//

include { GATK4_MUTECT2                    } from '../../../modules/nf-core/gatk4/mutect2/main'
include { GATK4_LEARNREADORIENTATIONMODEL  } from '../../../modules/nf-core/gatk4/learnreadorientationmodel/main'
include { GATK4_GETPILEUPSUMMARIES         } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION     } from '../../../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS          } from '../../../modules/nf-core/gatk4/filtermutectcalls/main'
include { GATK4_SELECTVARIANTS             } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { CURATE_VCF                       } from '../../../modules/local/curate_vcf/main'

workflow MUTECT_CALLING {

    take:
    ch_recal_bam_bai     // channel: [ val(meta), path(bam), path(bai) ]    from BQSR.recal_bam_bai
    ch_fasta             // channel: [ val(meta), path(fasta) ]
    ch_fai               // channel: [ val(meta), path(fai) ]
    ch_dict              // channel: [ val(meta), path(dict) ]
    ch_germline_resource // channel: path(vcf.gz)  (e.g. small_exac_common_3.hg38.vcf.gz)
    ch_germline_tbi      // channel: path(vcf.gz.tbi)

    main:

    ch_versions = Channel.empty()

    //
    // Step 17: Mutect2 tumour-only
    //
    // Mutect2 input: tuple(meta, input, input_index, intervals) + fasta/fai/dict tuples
    //                + germline_resource (plain path) + germline_tbi (plain path)
    //                + panel_of_normals (plain path, can be []) + pon_tbi
    //
    ch_mutect_input = ch_recal_bam_bai
        .map { meta, bam, bai -> [meta, bam, bai, []] }  // no intervals

    GATK4_MUTECT2(
        ch_mutect_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_germline_resource,
        ch_germline_tbi,
        [],   // panel_of_normals
        []    // panel_of_normals_tbi
    )
    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions)

    //
    // Step 18: LearnReadOrientationModel from Mutect2 f1r2
    //
    GATK4_LEARNREADORIENTATIONMODEL(GATK4_MUTECT2.out.f1r2)
    ch_versions = ch_versions.mix(GATK4_LEARNREADORIENTATIONMODEL.out.versions_gatk4)

    //
    // Step 19: GetPileupSummaries at common exac sites
    //
    // Uses the same recalibrated BAM, -V and -L both pointing at the
    // germline exac VCF per the legacy script.
    //
    ch_pileup_input = ch_recal_bam_bai
        .map { meta, bam, bai -> [meta, bam, bai, []] }  // no intervals

    GATK4_GETPILEUPSUMMARIES(
        ch_pileup_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_germline_resource,
        ch_germline_tbi
    )
    ch_versions = ch_versions.mix(GATK4_GETPILEUPSUMMARIES.out.versions_gatk4)

    //
    // Step 20: CalculateContamination — tumour-only (no matched normal)
    //
    ch_contam_input = GATK4_GETPILEUPSUMMARIES.out.table
        .map { meta, pileup -> [meta, pileup, []] }     // no matched-normal pileup

    GATK4_CALCULATECONTAMINATION(ch_contam_input)
    ch_versions = ch_versions.mix(GATK4_CALCULATECONTAMINATION.out.versions_gatk4)

    //
    // Step 21: FilterMutectCalls
    //
    // Module input tuple:
    //   (meta, vcf, vcf_tbi, stats, orientationbias, segmentation, table, estimate)
    //
    ch_filter_input = GATK4_MUTECT2.out.vcf
        .join(GATK4_MUTECT2.out.tbi)
        .join(GATK4_MUTECT2.out.stats)
        .join(GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior)
        .join(GATK4_CALCULATECONTAMINATION.out.segmentation)
        .join(GATK4_CALCULATECONTAMINATION.out.contamination)
        .map { meta, vcf, tbi, stats, bias, seg, table -> [meta, vcf, tbi, stats, bias, seg, table, []] }

    GATK4_FILTERMUTECTCALLS(
        ch_filter_input,
        ch_fasta,
        ch_fai,
        ch_dict
    )
    ch_versions = ch_versions.mix(GATK4_FILTERMUTECTCALLS.out.versions_gatk4)

    //
    // Step 22: SelectVariants — keep PASS-only, emit sites-only VCF
    //
    // Configure via ext.args in conf/modules.config:
    //     '--exclude-filtered true --sites-only-vcf-output true'
    //
    ch_select_input = GATK4_FILTERMUTECTCALLS.out.vcf
        .join(GATK4_FILTERMUTECTCALLS.out.tbi)
        .map { meta, vcf, tbi -> [meta, vcf, tbi, []] }

    GATK4_SELECTVARIANTS(ch_select_input)
    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions_gatk4)

    //
    // Step 23: curate_vcf (local) — produces unmasked + indel variants
    //
    // CURATE_VCF decompresses .vcf.gz internally if needed (added in the
    // module script). We pass the selectvariants VCF directly.
    //
    CURATE_VCF(GATK4_SELECTVARIANTS.out.vcf)
    ch_versions = ch_versions.mix(CURATE_VCF.out.versions)

    emit:
    unmasked_vcf  = CURATE_VCF.out.unmasked   // [meta, vcf]   — input to DB_CONSTRUCT (step 25)
    indel_vcf     = CURATE_VCF.out.indel      // [meta, vcf]   — input to DB_CONSTRUCT (step 25)
    filtered_vcf  = GATK4_FILTERMUTECTCALLS.out.vcf
    contamination = GATK4_CALCULATECONTAMINATION.out.contamination
    versions      = ch_versions
}
