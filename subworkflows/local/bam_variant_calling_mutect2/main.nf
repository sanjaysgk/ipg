//
// BAM_VARIANT_CALLING_MUTECT2: Mutect2 tumour-only somatic calling + contamination
//
// Steps 17–22 of the legacy cryptic peptide pipeline:
//   17. GATK4 Mutect2                  tumour-only (no matched normal)
//   18. GATK4 LearnReadOrientationModel
//   19. GATK4 GetPileupSummaries
//   20. GATK4 CalculateContamination
//   21. GATK4 FilterMutectCalls
//   22. GATK4 SelectVariants           --exclude-filtered --sites-only-vcf-output
// Curation (step 23) is now a separate subworkflow: VCF_CURATE.
//

include { GATK4_MUTECT2                   } from '../../../modules/nf-core/gatk4/mutect2/main'
include { GATK4_LEARNREADORIENTATIONMODEL } from '../../../modules/nf-core/gatk4/learnreadorientationmodel/main'
include { GATK4_GETPILEUPSUMMARIES        } from '../../../modules/nf-core/gatk4/getpileupsummaries/main'
include { GATK4_CALCULATECONTAMINATION    } from '../../../modules/nf-core/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS         } from '../../../modules/nf-core/gatk4/filtermutectcalls/main'
include { GATK4_SELECTVARIANTS            } from '../../../modules/nf-core/gatk4/selectvariants/main'

workflow BAM_VARIANT_CALLING_MUTECT2 {

    take:
    ch_recal_bam_bai     // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta             // channel: [ val(meta), path(fasta) ]
    ch_fai               // channel: [ val(meta), path(fai) ]
    ch_dict              // channel: [ val(meta), path(dict) ]
    ch_germline_resource // channel: path(vcf.gz)
    ch_germline_tbi      // channel: path(vcf.gz.tbi)

    main:

    ch_versions = Channel.empty()

    // Step 17: Mutect2 tumour-only
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

    // Step 18: LearnReadOrientationModel from Mutect2 f1r2
    GATK4_LEARNREADORIENTATIONMODEL(GATK4_MUTECT2.out.f1r2)

    // Step 19: GetPileupSummaries at common exac sites
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

    // Step 20: CalculateContamination — tumour-only
    ch_contam_input = GATK4_GETPILEUPSUMMARIES.out.table
        .map { meta, pileup -> [meta, pileup, []] }     // no matched-normal pileup

    GATK4_CALCULATECONTAMINATION(ch_contam_input)

    // Step 21: FilterMutectCalls
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

    // Step 22: SelectVariants — PASS-only, sites-only VCF
    ch_select_input = GATK4_FILTERMUTECTCALLS.out.vcf
        .join(GATK4_FILTERMUTECTCALLS.out.tbi)
        .map { meta, vcf, tbi -> [meta, vcf, tbi, []] }

    GATK4_SELECTVARIANTS(ch_select_input)

    emit:
    vcf           = GATK4_SELECTVARIANTS.out.vcf          // [meta, vcf]
    tbi           = GATK4_SELECTVARIANTS.out.tbi          // [meta, tbi]
    filtered_vcf  = GATK4_FILTERMUTECTCALLS.out.vcf
    contamination = GATK4_CALCULATECONTAMINATION.out.contamination
    versions      = ch_versions
}
