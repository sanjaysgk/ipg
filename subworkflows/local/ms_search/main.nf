/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MS_SEARCH — open-source MS database search subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    1. PREPARE_FASTA           target/decoy database
    2. MSFRAGGER               → calibrated mzMLs + PIN
    3. COMET + SAGE            parallel on calibrated mzMLs → PIN
    4. CONVERT_MZML            mzML → MGF + scans.pkl + index2scan.pkl
    5. CONVERT_PEAKS           PEAKS db.psms.csv → PIN (uses index2scan)
    6. MS2RESCORE (per engine) rescore PINs with spectrum predictions; runs
                               mokapot internally (no standalone mokapot pass)
    7. INTEGRATE_ENGINES       merge across engines at 1% peptide-level FDR
----------------------------------------------------------------------------------------
*/

include { PREPARE_FASTA                } from '../../../modules/local/prepare_fasta/main'
include { COMBINE_FASTA                } from '../../../modules/local/combine_fasta/main'
include { MSFRAGGER                    } from '../../../modules/local/msfragger/main'
include { COMET                        } from '../../../modules/local/comet/main'
include { SAGE                         } from '../../../modules/local/sage/main'
include { CONVERT_MZML                 } from '../../../modules/local/convert_mzml/main'
include { CONVERT_PEAKS                } from '../../../modules/local/convert_peaks/main'
include { MS2RESCORE as MS2RESCORE_MSFRAGGER } from '../../../modules/local/ms2rescore/main'
include { MS2RESCORE as MS2RESCORE_COMET     } from '../../../modules/local/ms2rescore/main'
include { MS2RESCORE as MS2RESCORE_SAGE      } from '../../../modules/local/ms2rescore/main'
include { MS2RESCORE as MS2RESCORE_PEAKS     } from '../../../modules/local/ms2rescore/main'
include { INTEGRATE_ENGINES            } from '../../../modules/local/integrate_engines/main'

workflow MS_SEARCH {

    take:
    ch_ms_data         // channel: [ val(meta), path(ms_files) ]
    ch_fasta           // path: cryptic FASTA database (from db_construct)
    ch_canonical_fasta // path: canonical proteome FASTA (unused in 'cryptic_only' mode)
    ch_engines         // val: list of engine names e.g. ['msfragger','comet','sage']
    ch_msfragger_jar   // path: MSFragger JAR (may be empty channel)

    main:

    ch_versions = Channel.empty()
    def engines    = params.ms_engines.tokenize(',').collect { it.trim() }
    def instrument = params.instrument ?: 'orbitrap'
    def mod_type   = params.mod_type

    //
    // STEP 1: Assemble the search DB per --db_search_mode, then build target-decoy.
    //   cryptic_only : cryptic DB alone (original behaviour)
    //   appended     : canonical + cryptic in one DB → class-separated FDR in INTEGRATE
    //   separate     : (follow-up) independent canonical/cryptic searches + reconcile
    //
    def db_mode = params.db_search_mode ?: 'appended'
    if (db_mode == 'cryptic_only') {
        ch_search_db = ch_fasta
    } else if (db_mode == 'appended') {
        COMBINE_FASTA(
            ch_canonical_fasta.combine(ch_fasta).map { canon, cryptic -> [canon, cryptic].flatten() }
        )
        ch_search_db = COMBINE_FASTA.out.fasta
        ch_versions  = ch_versions.mix(COMBINE_FASTA.out.versions)
    } else {
        error("--db_search_mode '${db_mode}' not implemented yet (separate-search " +
              "mode is a follow-up; use 'appended' or 'cryptic_only'). " +
              "See .claude/specs/canonical-db-search.md")
    }

    PREPARE_FASTA(ch_search_db)
    ch_versions  = ch_versions.mix(PREPARE_FASTA.out.versions)
    ch_tda_fasta = PREPARE_FASTA.out.fasta

    //
    // STEP 2a: MSFragger → calibrated mzMLs + PIN
    //
    ch_calibrated_mzml = Channel.empty()
    ch_msfragger_log   = Channel.empty()
    ch_pin_msfragger   = Channel.empty()

    if (engines.contains('msfragger')) {
        ch_msfragger_params = Channel.fromPath(
            "${projectDir}/assets/ms_search_params/${instrument}_${mod_type}.params",
            checkIfExists: true
        ).collect()
        MSFRAGGER(ch_ms_data, ch_tda_fasta, ch_msfragger_jar, ch_msfragger_params)
        ch_versions        = ch_versions.mix(MSFRAGGER.out.versions)
        ch_calibrated_mzml = MSFRAGGER.out.mzml
        ch_msfragger_log   = MSFRAGGER.out.log

        ch_pin_msfragger = MSFRAGGER.out.pin
    }

    // Always feed raw input to Comet/Sage — they handle their own calibration.
    // MSFragger's calibrated mzML output is reserved for CONVERT_MZML downstream.
    ch_engine_input = ch_ms_data

    //
    // STEP 2b: Comet
    //
    ch_pin_comet = Channel.empty()
    if (engines.contains('comet')) {
        ch_comet_params = Channel.fromPath(
            "${projectDir}/assets/ms_search_params/comet/${instrument}_${mod_type}.params",
            checkIfExists: true
        ).collect()
        COMET(ch_engine_input, ch_tda_fasta, ch_comet_params)
        ch_versions = ch_versions.mix(COMET.out.versions)

        ch_pin_comet = COMET.out.pin
    }

    //
    // STEP 2c: Sage
    //
    ch_pin_sage = Channel.empty()
    if (engines.contains('sage')) {
        ch_sage_params = Channel.fromPath(
            "${projectDir}/assets/ms_search_params/sage/${instrument}_${mod_type}.json",
            checkIfExists: true
        ).collect()
        ch_sage_log = engines.contains('msfragger')
            ? ch_msfragger_log.collect()
            : Channel.value(file("${projectDir}/assets/NO_FILE", checkIfExists: false))
        SAGE(ch_engine_input, ch_tda_fasta, ch_sage_params, ch_sage_log)
        ch_versions = ch_versions.mix(SAGE.out.versions)

        ch_pin_sage = SAGE.out.pin
    }

    //
    // STEP 3: CONVERT_MZML — fan out over mzMLs, emit MGF + scans pkl.
    //         Only needed for MS2Rescore and PEAKS branch; skip otherwise.
    //
    ch_scans_per_sample    = Channel.empty()
    ch_mgf_per_sample      = Channel.empty()
    ch_idx2scan_per_sample = Channel.empty()

    if (!params.skip_ms2rescore) {
        ch_one_mzml = ch_engine_input
            .flatMap { meta, files ->
                files instanceof List ? files.collect { [meta, it] } : [[meta, files]]
            }
        CONVERT_MZML(ch_one_mzml)
        ch_versions = ch_versions.mix(CONVERT_MZML.out.versions.first())

        ch_scans_per_sample = CONVERT_MZML.out.scans
            .map { meta, f -> [meta, f] }
            .groupTuple(by: 0)
        ch_mgf_per_sample = CONVERT_MZML.out.mgf
            .map { meta, f -> [meta, f] }
            .groupTuple(by: 0)
        ch_idx2scan_per_sample = CONVERT_MZML.out.index2scan
            .map { meta, f -> [meta, f] }
            .groupTuple(by: 0)
    }

    //
    // STEP 3b: PEAKS branch (optional).
    // User-supplied db.psms.csv is converted to PIN using index2scan.pkl,
    // then run through its own MOKAPOT instance.
    //
    ch_pin_peaks = Channel.empty()
    if (engines.contains('peaks')) {
        ch_peaks_csv = ch_ms_data
            .map { meta, _files -> [meta, file(params.peaks_psm_csv, checkIfExists: true)] }
        ch_peaks_in = ch_peaks_csv.join(ch_idx2scan_per_sample)
        CONVERT_PEAKS(
            ch_peaks_in.map { meta, csv, _pkls -> [meta, csv] },
            ch_peaks_in.map { meta, _csv, pkls -> pkls }.flatMap { it }.collect()
        )
        ch_versions = ch_versions.mix(CONVERT_PEAKS.out.versions)

        ch_pin_peaks = CONVERT_PEAKS.out.pin
    }

    ch_ms2rescore_cfg = Channel.fromPath(
        "${projectDir}/assets/ms_search_params/ms2rescore/${instrument}_${mod_type}.json",
        checkIfExists: true
    ).collect()

    //
    // STEP 4: MS2Rescore per engine.
    //
    ch_rescored = Channel.empty()

    // Rescoring needs scans covering every PSM in the mokapot output.
    // Rather than trust per-sample meta.join (which is flaky when
    // groupTuple + stageAs interact), broadcast the ENTIRE set of
    // per-sample scans + MGFs to every MS2Rescore invocation. The prep
    // script merges all scans into one dict and looks up by run_scan
    // key, so a superset is safe.
    ch_all_scans = ch_scans_per_sample.map { meta, f -> f }.flatMap { it }.collect()
    ch_all_mgfs  = ch_mgf_per_sample.map   { meta, f -> f }.flatMap { it }.collect()

    if (!params.skip_ms2rescore) {
        if (engines.contains('msfragger')) {
            MS2RESCORE_MSFRAGGER(ch_pin_msfragger, 'msfragger', ch_all_scans, ch_all_mgfs, ch_ms2rescore_cfg)
            ch_versions = ch_versions.mix(MS2RESCORE_MSFRAGGER.out.versions)
            ch_rescored = ch_rescored.mix(MS2RESCORE_MSFRAGGER.out.psms.map { meta, f -> [meta, 'msfragger', f] })
        }
        if (engines.contains('comet')) {
            MS2RESCORE_COMET(ch_pin_comet, 'comet', ch_all_scans, ch_all_mgfs, ch_ms2rescore_cfg)
            ch_versions = ch_versions.mix(MS2RESCORE_COMET.out.versions)
            ch_rescored = ch_rescored.mix(MS2RESCORE_COMET.out.psms.map { meta, f -> [meta, 'comet', f] })
        }
        if (engines.contains('sage')) {
            MS2RESCORE_SAGE(ch_pin_sage, 'sage', ch_all_scans, ch_all_mgfs, ch_ms2rescore_cfg)
            ch_versions = ch_versions.mix(MS2RESCORE_SAGE.out.versions)
            ch_rescored = ch_rescored.mix(MS2RESCORE_SAGE.out.psms.map { meta, f -> [meta, 'sage', f] })
        }
        if (engines.contains('peaks')) {
            MS2RESCORE_PEAKS(ch_pin_peaks, 'peaks', ch_all_scans, ch_all_mgfs, ch_ms2rescore_cfg)
            ch_versions = ch_versions.mix(MS2RESCORE_PEAKS.out.versions)
            ch_rescored = ch_rescored.mix(MS2RESCORE_PEAKS.out.psms.map { meta, f -> [meta, 'peaks', f] })
        }
    }

    //
    // STEP 5: INTEGRATE_ENGINES — merge per-engine MS2Rescore TSVs into one.
    // Skipped when --skip_ms2rescore is set (no rescored TSVs to merge).
    //
    ch_integrated_psms     = Channel.empty()
    ch_integrated_peptides = Channel.empty()
    ch_integrated_chimeric = Channel.empty()

    if (!params.skip_ms2rescore && !params.skip_integrate_engines) {
        ch_integrate_in = ch_rescored
            .groupTuple(by: 0)   // [meta, [engines...], [tsvs...]]
            .map { meta, eng_list, tsv_list -> [meta, tsv_list, eng_list] }

        INTEGRATE_ENGINES(
            ch_integrate_in.map { meta, tsvs, engs -> [meta, tsvs] },
            ch_integrate_in.map { meta, tsvs, engs -> engs },
            ch_tda_fasta
        )
        ch_versions            = ch_versions.mix(INTEGRATE_ENGINES.out.versions)
        ch_integrated_psms     = INTEGRATE_ENGINES.out.psms
        ch_integrated_peptides = INTEGRATE_ENGINES.out.peptides
        ch_integrated_chimeric = INTEGRATE_ENGINES.out.chimeric
    }

    emit:
    psms              = ch_integrated_psms
    peptides          = ch_integrated_peptides
    chimeric          = ch_integrated_chimeric
    calibrated_mzml   = ch_calibrated_mzml
    mgf               = ch_mgf_per_sample
    rescored_per_eng  = ch_rescored
    versions          = ch_versions
}
