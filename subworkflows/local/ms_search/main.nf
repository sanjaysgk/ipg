/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MS_SEARCH — open-source MS database search subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs search engines (MSFragger, Comet, Sage) on MS data against a
    target-decoy FASTA, then applies mokapot for FDR control.

    Sprint 1: MSFragger → mokapot (single engine)
    Sprint 2+: adds Comet, Sage, MS2Rescore, INTEGRATE_ENGINES
----------------------------------------------------------------------------------------
*/

include { PREPARE_FASTA } from '../../../modules/local/prepare_fasta/main'
include { MSFRAGGER     } from '../../../modules/local/msfragger/main'
include { MOKAPOT       } from '../../../modules/local/mokapot/main'

workflow MS_SEARCH {

    take:
    ch_ms_data       // channel: [ val(meta), path(ms_files) ]
    ch_fasta         // path: search FASTA database
    ch_engines       // val: list of engine names e.g. ['msfragger']
    ch_msfragger_jar // path: MSFragger JAR (may be empty channel)

    main:

    ch_versions = Channel.empty()

    //
    // STEP 1: Prepare target-decoy FASTA
    //
    PREPARE_FASTA(ch_fasta)
    ch_versions = ch_versions.mix(PREPARE_FASTA.out.versions)

    ch_tda_fasta = PREPARE_FASTA.out.fasta

    //
    // STEP 2: Run search engines
    // Currently: MSFragger only (Sprint 1)
    // Sprint 2+: Comet, Sage run in parallel on calibrated mzMLs
    //

    // MSFragger branch
    ch_msfragger_pin = Channel.empty()
    ch_calibrated_mzml = Channel.empty()

    if (params.engines.tokenize(',').contains('msfragger')) {

        // MSFragger params file — shipped in assets/ms_search_params/
        // For Sprint 1, use ext.args to pass params file via config
        ch_params = Channel.fromPath(
            "${projectDir}/assets/ms_search_params/${params.instrument ?: 'orbitrap'}_${params.mod_type}.params",
            checkIfExists: true
        )

        MSFRAGGER(
            ch_ms_data,
            ch_tda_fasta,
            ch_msfragger_jar,
            ch_params.collect()
        )
        ch_versions = ch_versions.mix(MSFRAGGER.out.versions)

        ch_msfragger_pin = MSFRAGGER.out.pin
        ch_calibrated_mzml = MSFRAGGER.out.mzml

        //
        // STEP 3: Mokapot FDR control on MSFragger results
        //
        MOKAPOT(
            ch_msfragger_pin,
            'msfragger'
        )
        ch_versions = ch_versions.mix(MOKAPOT.out.versions)
    }

    emit:
    psms             = params.engines.tokenize(',').contains('msfragger') ? MOKAPOT.out.psms : Channel.empty()
    peptides         = params.engines.tokenize(',').contains('msfragger') ? MOKAPOT.out.peptides : Channel.empty()
    calibrated_mzml  = ch_calibrated_mzml
    versions         = ch_versions
}
