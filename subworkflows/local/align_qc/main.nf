//
// ALIGN_QC: RNA-seq alignment and strandedness inference
//
// Implements steps 1–3 of the legacy cryptic peptide pipeline:
//   01. STAR two-pass alignment against the GRCh38 STAR index
//   02. RSeQC infer_experiment.py for strand inference
//   03. samtools sort (coordinate) + samtools index
//

include { STAR_ALIGN             } from '../../../modules/nf-core/star/align/main'
include { SAMTOOLS_SORT          } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX         } from '../../../modules/nf-core/samtools/index/main'
include { RSEQC_INFEREXPERIMENT  } from '../../../modules/nf-core/rseqc/inferexperiment/main'

workflow ALIGN_QC {

    take:
    ch_reads        // channel: [ val(meta), [ reads ] ]
    ch_star_index   // channel: [ val(meta2), path(star_index_dir) ]
    ch_gtf          // channel: [ val(meta3), path(gtf) ]
    ch_fasta        // channel: [ val(meta4), path(fasta), path(fai) ] — for samtools sort optional ref
    ch_rseqc_bed    // channel: path(bed)                              — single value, shared across samples

    main:

    // All modules in this subworkflow use the nf-core topic-channel
    // version pattern (`eval(...) , topic: versions`). Those are
    // auto-collected by Nextflow into the 'versions' topic and merged
    // into the global versions table by `softwareVersionsToYAML`
    // downstream — we do NOT mix them into a ch_versions file channel
    // because doing so breaks the YAML loader in utils_nfcore_pipeline.
    ch_versions = channel.empty()

    //
    // Step 01: STAR two-pass alignment
    //
    STAR_ALIGN(
        ch_reads,
        ch_star_index,
        ch_gtf,
        /* star_ignore_sjdbgtf = */ false
    )

    //
    // Step 03: sort by coordinate and index
    // (out of order with the legacy script, which does infer_experiment
    //  on the unsorted SAM — RSeQC actually wants a coordinate-sorted,
    //  indexed BAM, so we sort first.)
    //
    SAMTOOLS_SORT(
        STAR_ALIGN.out.bam,
        ch_fasta,
        /* index_format = */ ''
    )

    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    ch_sorted_bam_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.index)

    //
    // Step 02: infer strandedness from a sorted+indexed BAM
    //
    RSEQC_INFEREXPERIMENT(
        ch_sorted_bam_bai,
        ch_rseqc_bed
    )

    emit:
    star_bam_unsorted = STAR_ALIGN.out.bam              // [meta, bam]  — consumed by BAM_PREP (step 06+)
    sorted_bam        = SAMTOOLS_SORT.out.bam           // [meta, bam]  — consumed by TRANSCRIPT_ASSEMBLY
    sorted_bai        = SAMTOOLS_INDEX.out.index        // [meta, bai]
    sorted_bam_bai    = ch_sorted_bam_bai               // [meta, bam, bai] — for downstream use
    strandedness_txt  = RSEQC_INFEREXPERIMENT.out.txt   // [meta, txt]
    star_log          = STAR_ALIGN.out.log_final        // [meta, log]
    versions          = ch_versions
}
