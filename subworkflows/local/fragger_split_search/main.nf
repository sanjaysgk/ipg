/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FRAGGER_SPLIT_SEARCH — MSFragger split-database search (bioconda, no FragPipe)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Large no-enzyme cryptic DBs (~200K seqs) OOM MSFragger's in-memory fragment
    index. Split the DB into N chunks, search each chunk in --partial mode, then
    MERGE with full-DB e-value recompute -> statistically identical to an unsplit
    search (no FDR inflation). Drop-in for the single MSFragger arm of ms_search:
    emits [meta, [merged PINs]] per sample, which feeds the existing MS2Rescore.

        per DISTINCT db : SPLIT_FASTA(db, N) -> N chunks
        per (run x chunk): MSFRAGGER_PARTIAL --partial -> chunk_<id>/ (raw scores + histogram)
        per run         : MERGE_SPLIT_SEARCH (sum histograms -> full-DB e-value) -> merged PIN
        per sample      : collect run PINs -> [meta, [pins]]
----------------------------------------------------------------------------------------
*/

include { SPLIT_FASTA        } from '../../../modules/local/split_fasta/main'
include { MSFRAGGER_PARTIAL  } from '../../../modules/local/msfragger_partial/main'
include { MERGE_SPLIT_SEARCH } from '../../../modules/local/merge_split_search/main'

workflow FRAGGER_SPLIT_SEARCH {

    take:
    ch_ms_db_tda        // channel: [ val(meta), [mzml_files], path(tda) ]  per sample
    ch_msfragger_jar    // path:    MSFragger JAR or file('NO_FILE') (bioconda + --key)
    ch_msfragger_params // path:    MSFragger native params file
    num_slices          // val:     number of database chunks (plain int, > 1)

    main:
    ch_versions = Channel.empty()

    // 1. Split each DISTINCT target-decoy db into N chunks (once per db). We split
    //    the target-decoy db so every chunk carries its own targets+decoys; the
    //    merge's full-DB histogram is what restores correct statistics.
    ch_db_distinct = ch_ms_db_tda
        .map { _meta, _files, tda -> [ tda.toString(), tda ] }
        .unique { it[0] }

    SPLIT_FASTA( ch_db_distinct.map { key, tda -> [ [id: tda.baseName, key: key], tda ] }, num_slices )
    ch_versions = ch_versions.mix(SPLIT_FASTA.out.versions)

    // chunk fastas keyed by db path: [dbKey, chunk_id, chunk_fasta]
    ch_chunks = SPLIT_FASTA.out.fasta_chunks
        .flatMap { meta, fastas ->
            (fastas instanceof List ? fastas : [fastas]).collect { fa ->
                [ meta.key, fa.parent.name as int, fa ]   // split_db/<chunk_idx>/<db>.fasta
            }
        }

    // 2. Flatten samples to per-run (one spectrum file each), keyed by db, then
    //    cross each run with its db's chunks -> one --partial search per (run x chunk).
    ch_runs = ch_ms_db_tda
        .flatMap { meta, files, tda ->
            (files instanceof List ? files : [files]).collect { f -> [ tda.toString(), meta, f ] }
        }

    ch_partial_in = ch_runs
        .combine( ch_chunks, by: 0 )                  // [dbKey, meta, run_file, chunk_id, chunk_fasta]
        .map { _dbKey, meta, run_file, chunk_id, chunk_fasta ->
            def run   = run_file.baseName
            def pmeta = meta + [ original_id: meta.id, run: run, chunk_id: chunk_id,
                id: "${meta.id}__${run}__c${chunk_id}" ]
            [ pmeta, run_file, chunk_fasta ]
        }

    MSFRAGGER_PARTIAL(
        ch_partial_in.map { m, f, _fa  -> [ m, f ] },
        ch_partial_in.map { _m, _f, fa -> fa },
        ch_msfragger_jar,
        ch_msfragger_params
    )
    ch_versions = ch_versions.mix(MSFRAGGER_PARTIAL.out.versions)

    // 3. Merge per run: gather the N chunk dirs for each (sample, run) and recompute
    //    full-DB e-values. sample_name = run basename so the script finds <run>.* inside.
    ch_merge_in = MSFRAGGER_PARTIAL.out.chunk_dir
        .map { m, dir -> [ "${m.original_id}__${m.run}", m, dir ] }
        .groupTuple( by: 0, size: num_slices )
        .map { _runKey, metas, dirs ->
            def rmeta = metas[0].findAll { k, _v -> k != 'chunk_id' } + [ id: metas[0].original_id ]
            [ rmeta, dirs ]
        }

    MERGE_SPLIT_SEARCH( ch_merge_in, num_slices, ch_msfragger_jar, ch_msfragger_params )
    ch_versions = ch_versions.mix(MERGE_SPLIT_SEARCH.out.versions)

    // 4. Collect per-run merged PINs back into one [meta, [pins]] per sample —
    //    same shape the non-split MSFragger arm hands to MS2Rescore.
    ch_pin = MERGE_SPLIT_SEARCH.out.pin
        .map { m, pin -> [ m.original_id, m, pin ] }
        .groupTuple( by: 0 )
        .map { _sid, metas, pins ->
            def smeta = metas[0].findAll { k, _v -> !(k in ['run', 'original_id']) }
            [ smeta, pins ]
        }

    emit:
    pin      = ch_pin       // channel: [ val(meta), [merged PINs] ] per sample
    versions = ch_versions
}
