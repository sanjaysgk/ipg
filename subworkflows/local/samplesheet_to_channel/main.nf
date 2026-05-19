/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SAMPLESHEET_TO_CHANNEL — parse + validate + meta-tag input rows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Single entry point for samplesheet parsing across the three pipeline
    entry steps (db_construct, ms_search, post_ms). Validates that the
    required --*_input param is set for the chosen --step, runs the
    nf-schema parser, normalises meta tags, and emits one channel per
    entry shape (plus an empty channel for the unused shapes).

    Patterned on nf-core/sarek subworkflows/local/samplesheet_to_channel/main.nf
    — single subworkflow that catches user errors before any compute starts.
----------------------------------------------------------------------------------------
*/

include { samplesheetToList } from 'plugin/nf-schema'

workflow SAMPLESHEET_TO_CHANNEL {

    take:
    step                  // String : params.step
    input                 // Path   : params.input              (RNAseq samplesheet)
    ms_input              // Path   : params.ms_input           (MS samplesheet)
    post_ms_input         // Path   : params.post_ms_input      (post-MS samplesheet)
    schema_input          // Path   : projectDir-relative schema for input
    schema_ms_input       // Path   : projectDir-relative schema for ms_input
    schema_post_ms_input  // Path   : projectDir-relative schema for post_ms_input

    main:

    // ---- Cross-cutting validation ----------------------------------------------
    // Refuse to start if the active --step has no matching input.
    if ( step == 'db_construct' && !input ) {
        error "--step db_construct requires --input <samplesheet.csv> (RNAseq FASTQ rows)."
    }
    if ( step == 'ms_search' && !ms_input ) {
        error "--step ms_search requires --ms_input <samplesheet.csv> (MS file rows)."
    }
    if ( step == 'post_ms' && !post_ms_input ) {
        error "--step post_ms requires --post_ms_input <samplesheet.csv> (cryptic/uniprot PSM rows)."
    }

    // ---- RNAseq samplesheet (params.input) -------------------------------------
    if ( step == 'db_construct' ) {
        ch_rnaseq = Channel
            .fromList(samplesheetToList(input, schema_input))
            .map { meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end: true  ], [ file(fastq_1) ] ]
                } else {
                    return [ meta.id, meta + [ single_end: false ], [ file(fastq_1), file(fastq_2) ] ]
                }
            }
            .groupTuple()
            .map { _id, metas, fastqs -> [ metas[0], fastqs.flatten() ] }
    } else {
        ch_rnaseq = Channel.empty()
    }

    // ---- MS samplesheet (params.ms_input) --------------------------------------
    if ( step == 'ms_search' ) {
        ch_ms = Channel
            .fromList(samplesheetToList(ms_input, schema_ms_input))
            .map { row ->
                // schema_ms_input emits: meta, ms_file [, condition, fraction, replicate]
                // We keep meta as-is; downstream subworkflows can read meta.fraction etc.
                def meta    = row[0]
                def ms_file = row[1]
                [ meta, file(ms_file) ]
            }
    } else {
        ch_ms = Channel.empty()
    }

    // ---- Post-MS samplesheet (params.post_ms_input) ----------------------------
    if ( step == 'post_ms' ) {
        ch_post_ms = Channel
            .fromList(samplesheetToList(post_ms_input, schema_post_ms_input))
            .map { meta, cryptic_psm, uniprot_psm, cryptic_decoy_score, uniprot_decoy_score ->
                def new_meta = meta + [
                    cryptic_decoy_score: cryptic_decoy_score,
                    uniprot_decoy_score: uniprot_decoy_score,
                ]
                [ new_meta, file(cryptic_psm), file(uniprot_psm) ]
            }
    } else {
        ch_post_ms = Channel.empty()
    }

    emit:
    rnaseq  = ch_rnaseq    // channel: [ meta, [r1, r2] ]
    ms      = ch_ms        // channel: [ meta, ms_file ]
    post_ms = ch_post_ms   // channel: [ meta, cryptic_psm, uniprot_psm ]
}
