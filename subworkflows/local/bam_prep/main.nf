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

    // All modules in this subworkflow use the nf-core topic-channel
    // version pattern — auto-collected downstream. See align_qc for
    // rationale on why we do not mix them into ch_versions.
    ch_versions = channel.empty()

    // References are singletons — value channels so they broadcast to every
    // rep (queue channels would throttle the per-rep processes to one run).
    ch_fasta = ch_fasta.first()
    ch_fai   = ch_fai.first()
    ch_dict  = ch_dict.first()

    //
    // Step 07: sort aligned BAM by queryname for MergeBamAlignment
    //
    // We override ext.args in conf/modules.config to pass '-n' and ext.prefix
    // to add a '.qnsort' suffix.
    SAMTOOLS_SORT(
        ch_star_bam,
        // .first() re-broadcasts: combine() yields a queue channel even from
        // value inputs, which would throttle this to one rep.
        ch_fasta.combine(ch_fai.map { _meta, fai -> fai }).map { meta, fasta, fai -> [meta, fasta, fai] }.first(),
        /* index_format = */ ''
    )

    //
    // Step 08: FastqToSam -> unmapped BAM with read-group + sample metadata
    //
    GATK4_FASTQTOSAM(ch_reads)

    //
    // Step 09: merge aligned + unmapped BAMs
    //
    ch_merge_input = SAMTOOLS_SORT.out.bam.join(GATK4_FASTQTOSAM.out.bam)
    GATK4_MERGEBAMALIGNMENT(
        ch_merge_input,
        ch_fasta,
        ch_dict
    )

    //
    // Step 10: mark duplicates, merging a sample's technical reps here (rather
    // than carrying them separately through to BaseRecalibrator). Read groups
    // were assigned at Step 8 (FastqToSam), so MarkDuplicates can take one
    // --INPUT per rep BAM and emit a single merged BAM per sample. Group per-rep
    // BAMs by meta.sample; single-rep samples form a one-element list (the
    // original single-BAM path). The merged BAM then feeds BOTH transcript
    // assembly (StringTie) and the variant-calling branch. fasta/fai are bare
    // paths, not tuples.
    //
    ch_markdup_input = GATK4_MERGEBAMALIGNMENT.out.bam
        .map { meta, bam -> [ meta.sample ?: meta.id, meta, bam ] }
        .groupTuple()
        .map { sample, metas, bams ->
            def keep = metas[0].keySet() - ['rep', 'read_group', 'num_reps', 'sample']
            [ metas[0].subMap(keep) + [ id: sample ], bams ]
        }

    GATK4_MARKDUPLICATES(
        ch_markdup_input,
        ch_fasta.map { _meta, fasta -> fasta },
        ch_fai.map  { _meta, fai   -> fai   }
    )

    //
    // Index the MarkDuplicates output BAM so SplitNCigarReads has an index
    //
    SAMTOOLS_INDEX(GATK4_MARKDUPLICATES.out.bam)

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
    markdup_bam     = GATK4_MARKDUPLICATES.out.bam               // [meta, bam]  — merged, coord-sorted; feeds transcript assembly
    markdup_metrics = GATK4_MARKDUPLICATES.out.metrics           // [meta, metrics.txt]
    validate_report = GATK4_VALIDATESAMFILE.out.summary          // [meta, txt]
    versions        = ch_versions
}
