//
// DB_CONSTRUCT: build the final cryptic peptide FASTA database
//
// Implements steps 24–31 of the legacy cryptic peptide pipeline. This is the
// IPG-specific "cryptic peptide discovery" heart of the pipeline: two parallel
// paths fan out from the curated VCFs (unmasked + indel), each producing an
// alternate-reference FASTA → lifted GTF → transcriptome → three-frame
// peptide translation, and all three translations are then squished into one
// deduplicated cryptic peptide database.
//
//   24. GATK4 IndexFeatureFile             index each curated VCF
//   25. GATK4 FastaAlternateReferenceMaker apply variants to reference
//   26. revert_headers                     restore chromosome names
//   27. gff3sort                           sort assembly GTF for tabix
//   28. alt_liftover                       lift GTF to alt-reference coords
//   29. gffread                            extract transcriptome FASTA
//   30. triple_translate                   3-frame peptide translation
//   31. squish                             dedupe all three translations
//

include { GATK4_INDEXFEATUREFILE as INDEX_UNMASKED           } from '../../../modules/nf-core/gatk4/indexfeaturefile/main'
include { GATK4_INDEXFEATUREFILE as INDEX_INDEL              } from '../../../modules/nf-core/gatk4/indexfeaturefile/main'
include { GATK4_FASTAALTERNATEREFERENCEMAKER as FARM_UNMASKED } from '../../../modules/local/gatk4_fastaalternatereferencemaker/main'
include { GATK4_FASTAALTERNATEREFERENCEMAKER as FARM_INDEL    } from '../../../modules/local/gatk4_fastaalternatereferencemaker/main'
include { REVERT_HEADERS as REVERT_UNMASKED                  } from '../../../modules/local/revert_headers/main'
include { REVERT_HEADERS as REVERT_INDEL                     } from '../../../modules/local/revert_headers/main'
include { GFF3SORT                                           } from '../../../modules/local/gff3sort/main'
include { ALT_LIFTOVER as LIFTOVER_UNMASKED                  } from '../../../modules/local/alt_liftover/main'
include { ALT_LIFTOVER as LIFTOVER_INDEL                     } from '../../../modules/local/alt_liftover/main'
include { GFFREAD as GFFREAD_UNMASKED                        } from '../../../modules/nf-core/gffread/main'
include { GFFREAD as GFFREAD_INDEL                           } from '../../../modules/nf-core/gffread/main'
include { GFFREAD as GFFREAD_REFERENCE                       } from '../../../modules/nf-core/gffread/main'
include { TRIPLE_TRANSLATE as TT_UNMASKED                    } from '../../../modules/local/triple_translate/main'
include { TRIPLE_TRANSLATE as TT_INDEL                       } from '../../../modules/local/triple_translate/main'
include { TRIPLE_TRANSLATE as TT_REFERENCE                   } from '../../../modules/local/triple_translate/main'
include { SQUISH                                             } from '../../../modules/local/squish/main'

workflow DB_CONSTRUCT {

    take:
    ch_unmasked_vcf       // channel: [ val(meta), path(vcf) ]    from MUTECT_CALLING.unmasked_vcf
    ch_indel_vcf          // channel: [ val(meta), path(vcf) ]    from MUTECT_CALLING.indel_vcf
    ch_assembly_gtf       // channel: [ val(meta), path(gtf) ]    from TRANSCRIPT_ASSEMBLY.assembly_gtf
    ch_tracking           // channel: [ val(meta), path(tracking) ] from TRANSCRIPT_ASSEMBLY.tracking
    ch_fasta              // channel: [ val(meta), path(fasta) ]
    ch_fai                // channel: [ val(meta), path(fai) ]
    ch_dict               // channel: [ val(meta), path(dict) ]

    main:

    ch_versions = Channel.empty()

    //
    // Step 24: index both curated VCFs
    //
    INDEX_UNMASKED(ch_unmasked_vcf)
    INDEX_INDEL(ch_indel_vcf)

    //
    // Step 25: FastaAlternateReferenceMaker for each curated VCF
    //
    ch_farm_unmasked_in = ch_unmasked_vcf.join(INDEX_UNMASKED.out.index)
    ch_farm_indel_in    = ch_indel_vcf.join(INDEX_INDEL.out.index)

    FARM_UNMASKED(
        ch_farm_unmasked_in,
        ch_fasta,
        ch_fai,
        ch_dict
    )
    FARM_INDEL(
        ch_farm_indel_in,
        ch_fasta,
        ch_fai,
        ch_dict
    )

    //
    // Step 26: revert_headers — restore chromosome names on both alt genomes
    //
    REVERT_UNMASKED(
        FARM_UNMASKED.out.fasta,
        ch_fasta.map { _meta, fa -> fa }
    )
    REVERT_INDEL(
        FARM_INDEL.out.fasta,
        ch_fasta.map { _meta, fa -> fa }
    )
    ch_versions = ch_versions.mix(REVERT_UNMASKED.out.versions)
    ch_versions = ch_versions.mix(REVERT_INDEL.out.versions)

    //
    // Step 27: gff3sort the assembly GTF once (shared by both liftovers + ref)
    //
    GFF3SORT(ch_assembly_gtf)
    ch_versions = ch_versions.mix(GFF3SORT.out.versions)

    //
    // Step 28: alt_liftover — lift the sorted GTF into each alt-reference
    //
    // alt_liftover takes (meta, vcf, gtf) + a suffix value.
    //
    ch_liftover_unmasked_in = ch_unmasked_vcf.join(GFF3SORT.out.gtf)
    ch_liftover_indel_in    = ch_indel_vcf.join(GFF3SORT.out.gtf)

    LIFTOVER_UNMASKED(
        ch_liftover_unmasked_in,
        /* suffix = */ '_unmasked'
    )
    LIFTOVER_INDEL(
        ch_liftover_indel_in,
        /* suffix = */ '_indel'
    )
    ch_versions = ch_versions.mix(LIFTOVER_UNMASKED.out.versions)
    ch_versions = ch_versions.mix(LIFTOVER_INDEL.out.versions)

    //
    // Step 29: gffread extract transcriptome FASTA from each lifted GTF
    //
    GFFREAD_UNMASKED(
        LIFTOVER_UNMASKED.out.gtf,
        REVERT_UNMASKED.out.fasta.map { _meta, fa -> fa }
    )
    GFFREAD_INDEL(
        LIFTOVER_INDEL.out.gtf,
        REVERT_INDEL.out.fasta.map { _meta, fa -> fa }
    )
    // Reference transcriptome (no variants) — uses the sorted assembly GTF
    // against the original reference FASTA so the final squish has three
    // inputs matching the legacy step 29 third gffread invocation.
    GFFREAD_REFERENCE(
        GFF3SORT.out.gtf,
        ch_fasta.map { _meta, fa -> fa }
    )

    //
    // Step 30: triple_translate each transcriptome
    //
    // triple_translate takes a (meta, transcriptome_fasta, tracking) tuple.
    //
    ch_tt_unmasked_in  = GFFREAD_UNMASKED.out.gffread_fasta.join(ch_tracking)
    ch_tt_indel_in     = GFFREAD_INDEL.out.gffread_fasta.join(ch_tracking)
    ch_tt_reference_in = GFFREAD_REFERENCE.out.gffread_fasta.join(ch_tracking)

    TT_UNMASKED(ch_tt_unmasked_in)
    TT_INDEL(ch_tt_indel_in)
    TT_REFERENCE(ch_tt_reference_in)
    ch_versions = ch_versions.mix(TT_UNMASKED.out.versions)
    ch_versions = ch_versions.mix(TT_INDEL.out.versions)
    ch_versions = ch_versions.mix(TT_REFERENCE.out.versions)

    //
    // Step 31: squish — dedupe the three translations into the final DB
    //
    // Gather the three translated FASTAs into a single [meta, [fa1, fa2, fa3]]
    // channel. Join by meta so each sample squishes its own triple.
    //
    ch_squish_input = TT_REFERENCE.out.fasta
        .join(TT_UNMASKED.out.fasta)
        .join(TT_INDEL.out.fasta)
        .map { meta, ref, unmasked, indel -> [meta, [ref, unmasked, indel]] }

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
