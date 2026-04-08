//
// BAM_PREP: prepare a variant-calling-ready BAM from the STAR-aligned reads
//
// Implements steps 6–12 of the legacy cryptic peptide pipeline. Step 06 of
// the legacy script re-runs STAR on the same input FASTQs in a separate
// directory — we skip that duplication and consume the STAR output from
// ALIGN_QC instead.
//
//   07. samtools sort (queryname)    prepare aligned BAM for MergeBamAlignment
//   08. GATK4 FastqToSam             convert raw FASTQ to an unmapped BAM
//                                    containing sample / read group metadata
//   09. GATK4 MergeBamAlignment      merge the aligned BAM with the unmapped
//                                    BAM so alignments carry RG tags
//   10. GATK4 MarkDuplicates         mark duplicate reads
//   11. GATK4 SplitNCigarReads       split reads with N in CIGAR (RNA-seq)
//   12. GATK4 ValidateSamFile        SUMMARY validation of the final BAM
//

include { SAMTOOLS_SORT              } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX             } from '../../../modules/nf-core/samtools/index/main'
include { GATK4_FASTQTOSAM           } from '../../../modules/nf-core/gatk4/fastqtosam/main'
include { GATK4_MERGEBAMALIGNMENT    } from '../../../modules/nf-core/gatk4/mergebamalignment/main'
include { GATK4_MARKDUPLICATES       } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_SPLITNCIGARREADS     } from '../../../modules/nf-core/gatk4/splitncigarreads/main'
include { GATK4_VALIDATESAMFILE      } from '../../../modules/local/gatk4_validatesamfile/main'

workflow BAM_PREP {

    take:
    ch_star_bam         // channel: [ val(meta), path(bam) ]     from ALIGN_QC (coordinate-sorted is fine — GATK will re-sort)
    ch_reads            // channel: [ val(meta), [path(r1), path(r2)] ]
    ch_fasta            // channel: [ val(meta), path(fasta) ]
    ch_fai              // channel: [ val(meta), path(fai) ]
    ch_dict             // channel: [ val(meta), path(dict) ]

    main:

    ch_versions = Channel.empty()

    //
    // Step 07: sort aligned BAM by queryname for MergeBamAlignment
    //
    // We override ext.args in conf/modules.config to pass '-n' and ext.prefix
    // to add a '.qnsort' suffix.
    SAMTOOLS_SORT(
        ch_star_bam,
        ch_fasta.combine(ch_fai.map { _meta, fai -> fai }).map { meta, fasta, fai -> [meta, fasta, fai] },
        /* index_format = */ ''
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions_samtools)

    //
    // Step 08: FastqToSam -> unmapped BAM with read-group + sample metadata
    //
    GATK4_FASTQTOSAM(ch_reads)
    ch_versions = ch_versions.mix(GATK4_FASTQTOSAM.out.versions_gatk4)

    //
    // Step 09: merge aligned + unmapped BAMs
    //
    ch_merge_input = SAMTOOLS_SORT.out.bam.join(GATK4_FASTQTOSAM.out.bam)
    GATK4_MERGEBAMALIGNMENT(
        ch_merge_input,
        ch_fasta,
        ch_dict
    )
    ch_versions = ch_versions.mix(GATK4_MERGEBAMALIGNMENT.out.versions_gatk4)

    //
    // Step 10: mark duplicates
    // GATK4_MARKDUPLICATES takes bare path fasta / path fai (not tuples)
    //
    GATK4_MARKDUPLICATES(
        GATK4_MERGEBAMALIGNMENT.out.bam,
        ch_fasta.map { _meta, fasta -> fasta },
        ch_fai.map  { _meta, fai   -> fai   }
    )
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions_gatk4)

    //
    // Index the MarkDuplicates output BAM so SplitNCigarReads has an index
    //
    SAMTOOLS_INDEX(GATK4_MARKDUPLICATES.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions_samtools)

    //
    // Step 11: SplitNCigarReads — required for RNA-seq variant calling
    //
    ch_split_input = GATK4_MARKDUPLICATES.out.bam
        .join(SAMTOOLS_INDEX.out.index)
        .map { meta, bam, bai -> [meta, bam, bai, []] }   // no intervals

    GATK4_SPLITNCIGARREADS(
        ch_split_input,
        ch_fasta,
        ch_fai,
        ch_dict
    )
    ch_versions = ch_versions.mix(GATK4_SPLITNCIGARREADS.out.versions_gatk4)

    //
    // Step 12: ValidateSamFile (summary) — audit only, no channel output used downstream
    //
    GATK4_VALIDATESAMFILE(
        GATK4_SPLITNCIGARREADS.out.bam,
        ch_fasta,
        ch_fai,
        ch_dict
    )

    emit:
    splitncigar_bam = GATK4_SPLITNCIGARREADS.out.bam             // [meta, bam]  — input to BQSR
    markdup_metrics = GATK4_MARKDUPLICATES.out.metrics           // [meta, metrics.txt]
    validate_report = GATK4_VALIDATESAMFILE.out.summary          // [meta, txt]
    versions        = ch_versions
}
