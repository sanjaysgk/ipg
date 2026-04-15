/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    POST_MS_ANALYSIS — two-phase post-MS search analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Phase 1: db_compare → origins (simple mode)
    Phase 2: db_compare (with discard/unconventional) → origins (full Ensembl mode)

    Follows the workflow clarified by Kate Scull (2025-07-15).
----------------------------------------------------------------------------------------
*/

include { DB_COMPARE as DB_COMPARE_PHASE1 } from '../../../modules/local/db_compare/main'
include { DB_COMPARE as DB_COMPARE_PHASE2 } from '../../../modules/local/db_compare/main'
include { ORIGINS    as ORIGINS_SIMPLE     } from '../../../modules/local/origins/main'
include { ORIGINS    as ORIGINS_FULL       } from '../../../modules/local/origins/main'

workflow POST_MS_ANALYSIS {

    take:
    ch_post_ms      // channel: [ val(meta), path(cryptic_psm), path(uniprot_psm) ]
                    //   meta must include: id, cryptic_decoy_score, uniprot_decoy_score
    ch_tracking             // path: gffcompare .tracking file
    ch_uniprot_fasta        // path: UniProt FASTA
    ch_transcriptome_fasta  // path: transcriptome FASTA from DB construction

    main:

    ch_versions = Channel.empty()

    //
    // Phase 1: db_compare (no discard/unconventional lists — use NO_FILE
    // sentinel so the module's `discard_list.name != 'NO_FILE'` gate
    // evaluates false and the -j / -u flags are omitted. Passing bare `[]`
    // here makes Nextflow stage a file whose .name is not 'NO_FILE', which
    // produced `db_compare_v2.R -j  -u  ` and failed optparse.
    //
    DB_COMPARE_PHASE1(
        ch_post_ms,
        file("${projectDir}/assets/NO_FILE"),
        file("${projectDir}/assets/NO_FILE")
    )
    ch_versions = ch_versions.mix(DB_COMPARE_PHASE1.out.versions)

    //
    // Phase 1: origins in simple mode (-s flag set via ext.args in post_ms.config)
    //
    ORIGINS_SIMPLE(
        DB_COMPARE_PHASE1.out.cryptic_only,
        ch_tracking,
        ch_uniprot_fasta,
        ch_transcriptome_fasta
    )
    ch_versions = ch_versions.mix(ORIGINS_SIMPLE.out.versions)

    //
    // Phase 2: db_compare with discard + unconventional from Phase 1 origins
    //
    // Join the original PSM channels back with the origins output by meta.id
    // so each sample gets its matching discard/unconventional lists.
    //
    ch_phase2_input = ch_post_ms
        .map { meta, c, u -> [ meta.id, meta, c, u ] }
        .join(ORIGINS_SIMPLE.out.discard.map { meta, f -> [ meta.id, f ] })
        .join(ORIGINS_SIMPLE.out.unconventional.map { meta, f -> [ meta.id, f ] })
        .map { id, meta, c, u, d, uc -> [ meta, c, u, d, uc ] }
        // shape: [ meta, cryptic_psm, uniprot_psm, discard, unconventional ]

    DB_COMPARE_PHASE2(
        ch_phase2_input.map { meta, c, u, d, uc -> [ meta, c, u ] },
        ch_phase2_input.map { meta, c, u, d, uc -> d },
        ch_phase2_input.map { meta, c, u, d, uc -> uc }
    )
    ch_versions = ch_versions.mix(DB_COMPARE_PHASE2.out.versions)

    //
    // Phase 2: origins in full Ensembl mode (no -s flag)
    //
    ORIGINS_FULL(
        DB_COMPARE_PHASE2.out.unambiguous_unconventional,
        ch_tracking,
        ch_uniprot_fasta,
        ch_transcriptome_fasta
    )
    ch_versions = ch_versions.mix(ORIGINS_FULL.out.versions)

    emit:
    phase1_cryptic_only = DB_COMPARE_PHASE1.out.cryptic_only
    phase2_unconventional = DB_COMPARE_PHASE2.out.unambiguous_unconventional
    origins_full = ORIGINS_FULL.out.all_txt
    versions = ch_versions
}
