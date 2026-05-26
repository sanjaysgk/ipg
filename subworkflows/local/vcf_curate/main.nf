//
// VCF_CURATE: kescull curate_vcf — produces unmasked + indel variant tracks
//
// Step 23 of the legacy pipeline. Runs the curate_vcf C tool twice:
//   - without -d  -> *_indel.vcf     (deletions preserved; indel-aware track)
//   - with    -d  -> *_unmasked.vcf  (deletions deprioritised; substitution track)
// Both feed DB_CONSTRUCT's FastaAlternateReferenceMaker step.
//

include { CURATE_VCF } from '../../../modules/local/curate_vcf/main'

workflow VCF_CURATE {

    take:
    ch_vcf      // channel: [ val(meta), path(vcf) ]   from BAM_VARIANT_CALLING_MUTECT2.out.vcf

    main:

    ch_versions = Channel.empty()

    CURATE_VCF(ch_vcf)
    ch_versions = ch_versions.mix(CURATE_VCF.out.versions)

    emit:
    unmasked_vcf = CURATE_VCF.out.unmasked   // [meta, *_unmasked.vcf] -> DB_CONSTRUCT SNV track
    indel_vcf    = CURATE_VCF.out.indel      // [meta, *_indel.vcf]    -> DB_CONSTRUCT indel track
    curate_log   = CURATE_VCF.out.log
    versions     = ch_versions
}
