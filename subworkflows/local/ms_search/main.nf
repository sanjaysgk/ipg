/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MS_SEARCH — open-source MS database search subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    1. PREPARE_FASTA           target/decoy database (once per distinct db)
    2. MSFRAGGER               → calibrated mzMLs + PIN
    3. COMET + SAGE            parallel on calibrated mzMLs → PIN
    4. CONVERT_MZML            mzML → MGF + scans.pkl + index2scan.pkl
    5. CONVERT_PEAKS           PEAKS db.psms.csv → PIN (uses index2scan)
    6. MS2RESCORE (per engine) rescore PINs with spectrum predictions; runs
        mokapot internally (no standalone mokapot pass)
    7. INTEGRATE_ENGINES       merge across engines at the configured peptide-level FDR (default 1%)

    Per-sample DB routing: each input carries [meta, [ms_files], db]. All of a
    sample's MS files (fractions + replicates) are pooled into ONE search against
    that sample's ONE database — pooling maximises mokapot's target/decoy training
    data, and every PSM keeps its source `run` so per-replicate results stay
    recoverable downstream. The db is prepared (target-decoy, optional canonical
    append) ONCE per distinct database, then joined back to each sample by db path.
----------------------------------------------------------------------------------------
*/

include { THERMORAWFILEPARSER          } from '../../../modules/nf-core/thermorawfileparser/main'
include { MSCONVERT                    } from '../../../modules/local/msconvert/main'
include { PREPARE_FASTA                } from '../../../modules/local/prepare_fasta/main'
include { COMBINE_FASTA                } from '../../../modules/local/combine_fasta/main'
include { MSFRAGGER                    } from '../../../modules/local/msfragger/main'
include { FRAGGER_SPLIT_SEARCH         } from '../fragger_split_search/main'
include { COMET                        } from '../../../modules/local/comet/main'
include { SAGE                         } from '../../../modules/local/sage/main'
include { CONVERT_MZML                 } from '../../../modules/local/convert_mzml/main'
include { CONVERT_PEAKS                } from '../../../modules/local/convert_peaks/main'
include { MS2RESCORE as MS2RESCORE_MSFRAGGER } from '../../../modules/local/ms2rescore/main'
include { MS2RESCORE as MS2RESCORE_COMET     } from '../../../modules/local/ms2rescore/main'
include { MS2RESCORE as MS2RESCORE_SAGE      } from '../../../modules/local/ms2rescore/main'
include { MS2RESCORE as MS2RESCORE_PEAKS     } from '../../../modules/local/ms2rescore/main'
include { INTEGRATE_ENGINES            } from '../../../modules/local/integrate_engines/main'
include { PICKED_GROUP_FDR             } from '../../../modules/local/picked_group_fdr/main'

workflow MS_SEARCH {

    take:
    ch_ms_db           // channel: [ val(meta), [ms_files], path(db) ]  meta.db = resolved db path
    ch_canonical_fasta // path: canonical proteome FASTA (unused in 'cryptic_only' mode)
    ch_engines         // val: list of engine names e.g. ['msfragger','comet','sage']
    ch_msfragger_jar   // path: MSFragger JAR (may be empty channel)

    main:

    ch_versions = Channel.empty()
    def engines         = params.ms_engines.tokenize(',').collect { it.trim() }
    def instrument      = params.instrument ?: 'orbitrap'
    def mod_type        = params.mod_type
    def msfragger_split = (params.msfragger_split ?: 1) as int

    //
    // STEP 0: normalise spectra to mzML. Every engine (MSFragger/Comet/Sage) reads
    // mzML; Comet/Sage cannot read vendor formats. Detect format by extension per
    // file and convert the vendor ones:
    //   .raw (Thermo)  -> THERMORAWFILEPARSER -> mzML
    //   .d   (Bruker)  -> MSCONVERT           -> mzML
    //   .mzML/.mgf     -> passthrough
    // Convert per-file (unique id = basename, so outputs don't collide), carrying
    // the sample meta + db, then regroup by sample so the rest of the subworkflow
    // sees [meta, [mzML...], db] exactly as before.
    //
    ch_one = ch_ms_db.flatMap { meta, files, db ->
        (files instanceof List ? files : [files]).collect { f -> [ meta, db, f ] }
    }
    ch_fmt = ch_one.branch { _meta, _db, f ->
        raw:    f.name.toLowerCase().endsWith('.raw')
        bruker: f.name.toLowerCase().endsWith('.d')
        ready:  true
    }
    THERMORAWFILEPARSER( ch_fmt.raw.map    { meta, db, f -> [ [id: f.baseName, orig: meta, db: db], f ] } )
    MSCONVERT(           ch_fmt.bruker.map { meta, db, f -> [ [id: f.baseName, orig: meta, db: db], f ] } )
    ch_versions = ch_versions.mix(MSCONVERT.out.versions)
    ch_converted = THERMORAWFILEPARSER.out.spectra.map { m, mz -> [ m.orig, m.db, mz ] }
        .mix( MSCONVERT.out.spectra.map { m, mz -> [ m.orig, m.db, mz ] } )
    ch_ms_db = ch_converted
        .mix( ch_fmt.ready.map { meta, db, f -> [ meta, db, f ] } )
        .map { meta, db, f -> [ meta.id, meta, db, f ] }
        .groupTuple( by: 0 )
        .map { _id, metas, dbs, fs -> [ metas[0], fs, dbs[0] ] }   // [meta, [mzML files], db]

    // [meta, files] view for spectrum-side processes (engines, convert, peaks).
    ch_ms_data = ch_ms_db.map { meta, files, _db -> [ meta, files ] }

    //
    // STEP 1: Build a target-decoy search DB per DISTINCT database, then pair it
    // back to each sample. Key on the db path so samples sharing a db prep once.
    //   cryptic_only : cryptic DB alone
    //   appended     : canonical + cryptic in one DB → class-separated FDR
    //   separate     : (follow-up) independent canonical/cryptic searches + reconcile
    //
    def db_mode = params.db_search_mode ?: 'appended'

    ch_db_distinct = ch_ms_db
        .map { _meta, _files, db -> [ db.toString(), db ] }
        .unique { it[0] }

    if (db_mode == 'cryptic_only') {
        ch_prep_in = ch_db_distinct.map { key, db -> [ [id: db.baseName, key: key], db ] }
    } else if (db_mode == 'appended') {
        COMBINE_FASTA(
            ch_db_distinct.combine(ch_canonical_fasta).map { key, db, canon ->
                [ [id: db.baseName, key: key], [canon, db].flatten() ]
            }
        )
        ch_versions = ch_versions.mix(COMBINE_FASTA.out.versions)
        ch_prep_in  = COMBINE_FASTA.out.fasta
    } else {
        error("--db_search_mode '${db_mode}' not implemented yet (separate-search " +
            "mode is a follow-up; use 'appended' or 'cryptic_only').")
    }

    PREPARE_FASTA(ch_prep_in)
    ch_versions  = ch_versions.mix(PREPARE_FASTA.out.versions)
    ch_tda_keyed = PREPARE_FASTA.out.fasta.map { meta, tda -> [ meta.key, tda ] }   // [dbPathKey, tda]

    // Pair every sample group to its prepared target-decoy db (N samples : 1 db).
    ch_search = ch_ms_db
        .map { meta, files, db -> [ db.toString(), meta, files ] }
        .combine( ch_tda_keyed, by: 0 )
        .map { _key, meta, files, tda -> [ meta, files, tda ] }

    // [meta.id, tda] lookup for per-sample consumers downstream (INTEGRATE).
    ch_tda_by_id = ch_search.map { meta, _files, tda -> [ meta.id, tda ] }

    //
    // STEP 2a: MSFragger → calibrated mzMLs + PIN. Each sample's pooled files are
    // searched against its own target-decoy db (fork the joined channel by field).
    //
    ch_calibrated_mzml = Channel.empty()
    ch_msfragger_log   = Channel.empty()
    ch_pin_msfragger   = Channel.empty()

    if (engines.contains('msfragger')) {
        ch_msfragger_params = Channel.fromPath(
            "${projectDir}/assets/ms_search_params/msfragger/${instrument}_${mod_type}.params",
            checkIfExists: true
        ).collect()

        if (msfragger_split > 1) {
            // Split-database search: large no-enzyme cryptic DBs (~200K seqs) OOM
            // MSFragger's in-memory fragment index. Split into N chunks, search each
            // --partial, then merge with full-DB e-value recompute (statistically
            // identical to an unsplit search). Needs the REAL MSFragger (--partial +
            // --generate_expect_functions): a user JAR, or bioconda with the
            // MSFRAGGER_LICENSE secret (the modules' `secret` directive errors clearly
            // if it is unset).
            FRAGGER_SPLIT_SEARCH(ch_search, ch_msfragger_jar, ch_msfragger_params, msfragger_split)
            ch_versions      = ch_versions.mix(FRAGGER_SPLIT_SEARCH.out.versions)
            ch_pin_msfragger = FRAGGER_SPLIT_SEARCH.out.pin
        } else {
            MSFRAGGER(
                ch_search.map { meta, files, _tda -> [ meta, files ] },
                ch_search.map { _meta, _files, tda -> tda },
                ch_msfragger_jar,
                ch_msfragger_params
            )
            ch_versions        = ch_versions.mix(MSFRAGGER.out.versions)
            ch_calibrated_mzml = MSFRAGGER.out.mzml
            ch_msfragger_log   = MSFRAGGER.out.log
            ch_pin_msfragger   = MSFRAGGER.out.pin
        }
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
        COMET(
            ch_search.map { meta, files, _tda -> [ meta, files ] },
            ch_search.map { _meta, _files, tda -> tda },
            ch_comet_params
        )
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
        // SAGE optionally reuses THIS sample's MSFragger calibration (fragment tol /
        // topN) when MSFragger ran un-split. Join the log per sample by meta.id — never
        // collect() across samples: every MSFragger emits a log named search_log.txt, so
        // collecting >1 sample's logs into one SAGE input collides on filename.
        def no_file = file("${projectDir}/assets/NO_FILE", checkIfExists: false)
        ch_sage_in = (engines.contains('msfragger') && msfragger_split == 1)
            ? ch_search.map { meta, files, tda -> [ meta.id, meta, files, tda ] }
                .join( ch_msfragger_log.map { meta, log -> [ meta.id, log ] } )
                .map { _id, meta, files, tda, log -> [ meta, files, tda, log ] }
            : ch_search.map { meta, files, tda -> [ meta, files, tda, no_file ] }
        SAGE(
            ch_sage_in.map { meta, files, _tda, _log -> [ meta, files ] },
            ch_sage_in.map { _meta, _files, tda, _log -> tda },
            ch_sage_params,
            ch_sage_in.map { _meta, _files, _tda, log -> log }
        )
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
    // [meta, engine, mokapot_psms.txt, mokapot_decoy.txt] for tier-2 protein-group FDR.
    ch_mokapot  = Channel.empty()

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
            ch_mokapot  = ch_mokapot.mix(MS2RESCORE_MSFRAGGER.out.mokapot_psms.join(MS2RESCORE_MSFRAGGER.out.mokapot_decoy).map { meta, p, d -> [meta, 'msfragger', p, d] })
        }
        if (engines.contains('comet')) {
            MS2RESCORE_COMET(ch_pin_comet, 'comet', ch_all_scans, ch_all_mgfs, ch_ms2rescore_cfg)
            ch_versions = ch_versions.mix(MS2RESCORE_COMET.out.versions)
            ch_rescored = ch_rescored.mix(MS2RESCORE_COMET.out.psms.map { meta, f -> [meta, 'comet', f] })
            ch_mokapot  = ch_mokapot.mix(MS2RESCORE_COMET.out.mokapot_psms.join(MS2RESCORE_COMET.out.mokapot_decoy).map { meta, p, d -> [meta, 'comet', p, d] })
        }
        if (engines.contains('sage')) {
            MS2RESCORE_SAGE(ch_pin_sage, 'sage', ch_all_scans, ch_all_mgfs, ch_ms2rescore_cfg)
            ch_versions = ch_versions.mix(MS2RESCORE_SAGE.out.versions)
            ch_rescored = ch_rescored.mix(MS2RESCORE_SAGE.out.psms.map { meta, f -> [meta, 'sage', f] })
            ch_mokapot  = ch_mokapot.mix(MS2RESCORE_SAGE.out.mokapot_psms.join(MS2RESCORE_SAGE.out.mokapot_decoy).map { meta, p, d -> [meta, 'sage', p, d] })
        }
        if (engines.contains('peaks')) {
            MS2RESCORE_PEAKS(ch_pin_peaks, 'peaks', ch_all_scans, ch_all_mgfs, ch_ms2rescore_cfg)
            ch_versions = ch_versions.mix(MS2RESCORE_PEAKS.out.versions)
            ch_rescored = ch_rescored.mix(MS2RESCORE_PEAKS.out.psms.map { meta, f -> [meta, 'peaks', f] })
            ch_mokapot  = ch_mokapot.mix(MS2RESCORE_PEAKS.out.mokapot_psms.join(MS2RESCORE_PEAKS.out.mokapot_decoy).map { meta, p, d -> [meta, 'peaks', p, d] })
        }
    }

    //
    // STEP 5b: optional tier-2 protein/ORF-group FDR (picked) on the cryptic class.
    // Per engine, restrict the mokapot target+decoy PSMs to cryptic and run
    // PickedGroupFDR against the per-sample cryptic ORF FASTA (= the search db).
    //
    ch_protein_groups = Channel.empty()
    if (params.run_protein_group_fdr && !params.skip_ms2rescore) {
        ch_db_by_id = ch_ms_db.map { meta, _files, db -> [ meta.id, db ] }
        ch_pgfdr_in = ch_mokapot
            .map { meta, eng, p, d -> [ meta.id, meta, eng, p, d ] }
            .combine(ch_db_by_id, by: 0)
            .map { _id, meta, eng, p, d, db -> [ meta, eng, p, d, db ] }
        PICKED_GROUP_FDR(ch_pgfdr_in)
        ch_versions       = ch_versions.mix(PICKED_GROUP_FDR.out.versions)
        ch_protein_groups = PICKED_GROUP_FDR.out.protein_groups
    }

    //
    // STEP 5: INTEGRATE_ENGINES — merge per-engine MS2Rescore TSVs into one.
    // Skipped when --skip_ms2rescore is set (no rescored TSVs to merge). Each
    // sample integrates against its own target-decoy db (joined by meta.id).
    //
    ch_integrated_psms     = Channel.empty()
    ch_integrated_peptides = Channel.empty()
    ch_integrated_chimeric = Channel.empty()

    if (!params.skip_ms2rescore && !params.skip_integrate_engines) {
        ch_integrate_in = ch_rescored
            .groupTuple(by: 0)   // [meta, [engines...], [tsvs...]]
            .map { meta, eng_list, tsv_list -> [ meta.id, meta, tsv_list, eng_list ] }
            .join( ch_tda_by_id )   // [id, meta, tsvs, engines, tda]

        INTEGRATE_ENGINES(
            ch_integrate_in.map { _id, meta, tsvs, _engs, _tda -> [meta, tsvs] },
            ch_integrate_in.map { _id, _meta, _tsvs, engs, _tda -> engs },
            ch_integrate_in.map { _id, _meta, _tsvs, _engs, tda -> tda }
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
    protein_groups    = ch_protein_groups   // [meta, engine, *_protein_groups.tsv] — cryptic picked group FDR
    search_db         = ch_search.map { meta, _files, tda -> [ meta, tda ] }   // [meta, tda] per sample
    versions          = ch_versions
}
