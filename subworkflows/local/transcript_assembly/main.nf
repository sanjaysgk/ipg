//
// TRANSCRIPT_ASSEMBLY: de-novo transcript assembly and comparison to reference
//
// Implements steps 4–5 of the legacy cryptic peptide pipeline:
//   04. StringTie   assembles novel transcripts from the coordinate-sorted
//                   BAM (strand-aware via meta.strandedness)
//   05. gffcompare  compares the assembly to the reference annotation,
//                   producing the .tracking file consumed downstream by
//                   triple_translate and the .combined.gtf used for liftover
//

include { STRINGTIE_STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie/main'
include { GFFCOMPARE          } from '../../../modules/nf-core/gffcompare/main'

workflow TRANSCRIPT_ASSEMBLY {

    take:
    ch_sorted_bam       // channel: [ val(meta), path(bam) ]    meta.strandedness must be set
    ch_reference_gtf    // channel: [ val(meta), path(gtf) ]
    ch_fasta_fai        // channel: [ val(meta), path(fasta), path(fai) ]

    main:

    ch_versions = Channel.empty()

    //
    // Step 04: StringTie assembly
    //
    STRINGTIE_STRINGTIE(
        ch_sorted_bam,
        ch_reference_gtf.map { _meta, gtf -> gtf }    // StringTie takes a bare path
    )
    ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions_stringtie)

    //
    // Step 05: GFFCompare against the reference annotation
    //
    GFFCOMPARE(
        STRINGTIE_STRINGTIE.out.transcript_gtf,
        ch_fasta_fai,
        ch_reference_gtf
    )
    ch_versions = ch_versions.mix(GFFCOMPARE.out.versions_gffcompare)

    emit:
    assembly_gtf = STRINGTIE_STRINGTIE.out.transcript_gtf   // [meta, gtf]  — input to gff3sort in DB_CONSTRUCT
    combined_gtf = GFFCOMPARE.out.combined_gtf              // [meta, gtf]
    tracking     = GFFCOMPARE.out.tracking                  // [meta, tracking] — input to triple_translate
    stats        = GFFCOMPARE.out.stats                     // [meta, stats]
    versions     = ch_versions
}
