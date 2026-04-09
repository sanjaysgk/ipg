//
// DB_CONSTRUCT: build the final cryptic peptide FASTA database
//
// Implements steps 24–31 of the legacy cryptic peptide pipeline. The
// reference branch (gff3sort -> gffread -> triple_translate) ALWAYS runs.
// The variant-derived branches (FastaAlternateReferenceMaker -> revert_headers
// -> alt_liftover -> gffread -> triple_translate, fanned out for both the
// `unmasked` and `indel` curate_vcf outputs) only run when the global
// pipeline parameter `params.include_variant_peptides` is true.
//
// This is intentionally a pipeline-level flag, NOT a per-sample column.
// Mixing samples with different cryptic-DB intents in one run is confusing
// and error-prone — run the pipeline twice if you need both modes for a
// mixed cohort.
//
// SCIENTIFIC NOTE: variant calling itself is ALWAYS performed in tumour-only
// Mutect2 mode against a gnomAD-style germline allele-frequency database.
// This pipeline does NOT support matched tumour-normal calling. The
// include_variant_peptides flag controls only whether the discovered
// variants get folded into the final cryptic peptide DB. The legacy Scull
// et al. 2021 / D122_Lung pipeline always used reference-only (verified
// empirically against the legacy squish.log) — enabling include_variant_peptides
// is an EXTENSION suitable for samples expected to harbour biologically
// meaningful somatic variants (e.g. tumour tissue, hypermutated samples).
//
// Step mapping:
//   24. GATK4 IndexFeatureFile             index each curated VCF       (variant branch only)
//   25. GATK4 FastaAlternateReferenceMaker apply variants to reference  (variant branch only)
//   26. revert_headers                     restore chromosome names     (variant branch only)
//   27. gff3sort                           sort assembly GTF             (always)
//   28. alt_liftover                       lift GTF to alt-reference    (variant branch only)
//   29. gffread                            extract transcriptome FASTA  (always for ref; variant branch for unmasked/indel)
//   30. triple_translate                   3-frame peptide translation  (always for ref; variant branch for unmasked/indel)
//   31. squish                             dedupe into final DB         (always)
//

include { GATK4_INDEXFEATUREFILE as INDEX_UNMASKED            } from '../../../modules/nf-core/gatk4/indexfeaturefile/main'
include { GATK4_INDEXFEATUREFILE as INDEX_INDEL               } from '../../../modules/nf-core/gatk4/indexfeaturefile/main'
include { GATK4_FASTAALTERNATEREFERENCEMAKER as FARM_UNMASKED } from '../../../modules/local/gatk4_fastaalternatereferencemaker/main'
include { GATK4_FASTAALTERNATEREFERENCEMAKER as FARM_INDEL    } from '../../../modules/local/gatk4_fastaalternatereferencemaker/main'
include { REVERT_HEADERS as REVERT_UNMASKED                   } from '../../../modules/local/revert_headers/main'
include { REVERT_HEADERS as REVERT_INDEL                      } from '../../../modules/local/revert_headers/main'
include { GFF3SORT                                            } from '../../../modules/local/gff3sort/main'
include { ALT_LIFTOVER as LIFTOVER_UNMASKED                   } from '../../../modules/local/alt_liftover/main'
include { ALT_LIFTOVER as LIFTOVER_INDEL                      } from '../../../modules/local/alt_liftover/main'
include { GFFREAD as GFFREAD_UNMASKED                         } from '../../../modules/nf-core/gffread/main'
include { GFFREAD as GFFREAD_INDEL                            } from '../../../modules/nf-core/gffread/main'
include { GFFREAD as GFFREAD_REFERENCE                        } from '../../../modules/nf-core/gffread/main'
include { TRIPLE_TRANSLATE as TT_UNMASKED                     } from '../../../modules/local/triple_translate/main'
include { TRIPLE_TRANSLATE as TT_INDEL                        } from '../../../modules/local/triple_translate/main'
include { TRIPLE_TRANSLATE as TT_REFERENCE                    } from '../../../modules/local/triple_translate/main'
include { SQUISH                                              } from '../../../modules/local/squish/main'

workflow DB_CONSTRUCT {

    take:
    ch_unmasked_vcf       // channel: [ val(meta), path(vcf) ]      from MUTECT_CALLING.unmasked_vcf
    ch_indel_vcf          // channel: [ val(meta), path(vcf) ]      from MUTECT_CALLING.indel_vcf
    ch_assembly_gtf       // channel: [ val(meta), path(gtf) ]      from TRANSCRIPT_ASSEMBLY.combined_gtf
    ch_tracking           // channel: [ val(meta), path(tracking) ] from TRANSCRIPT_ASSEMBLY.tracking
    ch_fasta              // channel: [ val(meta), path(fasta) ]
    ch_fai                // channel: [ val(meta), path(fai) ]
    ch_dict               // channel: [ val(meta), path(dict) ]

    main:

    ch_versions = channel.empty()

    //
    // Step 27 (always): gff3sort the assembly GTF
    //
    GFF3SORT(ch_assembly_gtf)
    ch_versions = ch_versions.mix(GFF3SORT.out.versions)

    //
    // Reference branch (always runs for every sample) ----------------------
    //
    // Step 29 (reference): gffread extract transcriptome FASTA from sorted
    //                      assembly GTF against the original reference FASTA
    //
    GFFREAD_REFERENCE(
        GFF3SORT.out.gtf,
        ch_fasta.map { _meta, fa -> fa }
    )

    // Step 30 (reference): triple_translate the reference transcriptome
    //                      using the gffcompare .tracking file for class
    //                      codes and ENSG/ENST annotations
    ch_tt_reference_in = GFFREAD_REFERENCE.out.gffread_fasta.join(ch_tracking)
    TT_REFERENCE(ch_tt_reference_in)
    ch_versions = ch_versions.mix(TT_REFERENCE.out.versions)

    //
    // Variant branches (only when params.include_variant_peptides == true) ---
    //
    // Branch the curated VCFs by include_variant_peptides. Samples that
    // do NOT request variant peptides bypass the entire FARM/revert/liftover/
    // gffread/tt-unmasked-indel chain. Their cryptic.fasta is built from
    // the reference branch only — matching legacy Scull et al. 2021.
    //
    ch_unmasked_vcf
        .branch { meta, _vcf ->
            with_variants: params.include_variant_peptides == true
            reference_only: true
        }
        .set { ch_unmasked_branch }

    ch_indel_vcf
        .branch { meta, _vcf ->
            with_variants: params.include_variant_peptides == true
            reference_only: true
        }
        .set { ch_indel_branch }

    //
    // Step 24: index curated VCFs (variant branch only)
    //
    INDEX_UNMASKED(ch_unmasked_branch.with_variants)
    INDEX_INDEL(ch_indel_branch.with_variants)

    //
    // Step 25: FastaAlternateReferenceMaker (variant branch only)
    //
    ch_farm_unmasked_in = ch_unmasked_branch.with_variants.join(INDEX_UNMASKED.out.index)
    ch_farm_indel_in    = ch_indel_branch.with_variants.join(INDEX_INDEL.out.index)

    FARM_UNMASKED(ch_farm_unmasked_in, ch_fasta, ch_fai, ch_dict)
    FARM_INDEL(ch_farm_indel_in,       ch_fasta, ch_fai, ch_dict)

    //
    // Step 26: revert_headers — restore chromosome names on alt genomes
    //
    REVERT_UNMASKED(FARM_UNMASKED.out.fasta, ch_fasta.map { _meta, fa -> fa })
    REVERT_INDEL(FARM_INDEL.out.fasta,       ch_fasta.map { _meta, fa -> fa })
    ch_versions = ch_versions.mix(REVERT_UNMASKED.out.versions)
    ch_versions = ch_versions.mix(REVERT_INDEL.out.versions)

    //
    // Step 28: alt_liftover — lift sorted GTF into each alt-reference
    //
    ch_liftover_unmasked_in = ch_unmasked_branch.with_variants.join(GFF3SORT.out.gtf)
    ch_liftover_indel_in    = ch_indel_branch.with_variants.join(GFF3SORT.out.gtf)

    LIFTOVER_UNMASKED(ch_liftover_unmasked_in, /* suffix = */ '_unmasked')
    LIFTOVER_INDEL(ch_liftover_indel_in,       /* suffix = */ '_indel')
    ch_versions = ch_versions.mix(LIFTOVER_UNMASKED.out.versions)
    ch_versions = ch_versions.mix(LIFTOVER_INDEL.out.versions)

    //
    // Step 29 (variant branches): gffread extract transcriptome FASTAs
    // from each lifted GTF against the corresponding alt-reference FASTA
    //
    GFFREAD_UNMASKED(
        LIFTOVER_UNMASKED.out.gtf,
        REVERT_UNMASKED.out.fasta.map { _meta, fa -> fa }
    )
    GFFREAD_INDEL(
        LIFTOVER_INDEL.out.gtf,
        REVERT_INDEL.out.fasta.map { _meta, fa -> fa }
    )

    //
    // Step 30 (variant branches): triple_translate each variant transcriptome
    //
    ch_tt_unmasked_in = GFFREAD_UNMASKED.out.gffread_fasta.join(ch_tracking)
    ch_tt_indel_in    = GFFREAD_INDEL.out.gffread_fasta.join(ch_tracking)

    TT_UNMASKED(ch_tt_unmasked_in)
    TT_INDEL(ch_tt_indel_in)
    ch_versions = ch_versions.mix(TT_UNMASKED.out.versions)
    ch_versions = ch_versions.mix(TT_INDEL.out.versions)

    //
    // Step 31: squish — dedupe peptide DBs into the final cryptic FASTA
    //
    // When params.include_variant_peptides is true, gather all three
    // translation outputs (reference + unmasked + indel) keyed by meta.id.
    // When false, the unmasked/indel branches never ran, so the join below
    // is a no-op and squish receives only the reference translation.
    //
    // Implementation: build the squish input as `reference fasta + optional
    // variant fastas`. Use groupTuple keyed by meta.id to collect everything
    // for one sample into a single SQUISH invocation.
    //
    ch_ref_only = TT_REFERENCE.out.fasta
        .map { meta, fa -> [meta, [fa]] }

    ch_variant_combined = TT_UNMASKED.out.fasta
        .join(TT_INDEL.out.fasta)
        .map { meta, unmasked, indel -> [meta, [unmasked, indel]] }

    // Outer-join reference (always) with variant pair (only for variant
    // samples). For reference-only samples, the variant side will be null
    // and we use just the reference fasta.
    ch_squish_input = ch_ref_only
        .join(ch_variant_combined, remainder: true)
        .map { meta, ref_list, variant_list ->
            def fastas = variant_list ? ref_list + variant_list : ref_list
            [meta, fastas]
        }

    SQUISH(ch_squish_input)
    ch_versions = ch_versions.mix(SQUISH.out.versions)

    emit:
    cryptic_fasta       = SQUISH.out.fasta                      // [meta, fasta]  — THE DELIVERABLE
    alt_fasta_unmasked  = REVERT_UNMASKED.out.fasta
    alt_fasta_indel     = REVERT_INDEL.out.fasta
    sorted_gtf          = GFF3SORT.out.gtf
    lifted_gtf_unmasked = LIFTOVER_UNMASKED.out.gtf
    lifted_gtf_indel    = LIFTOVER_INDEL.out.gtf
    versions            = ch_versions
}
