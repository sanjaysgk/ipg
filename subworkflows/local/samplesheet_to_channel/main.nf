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
    // A sample may appear on multiple rows (lanes / technical reps). Each row
    // is one rep, aligned separately with its own read group (rg1, rg2... by
    // fastq name so STAR and FastqToSam agree), then merged at MarkDuplicates.
    // A single-row sample keeps meta.id = sample, unchanged from before.
    if ( step == 'db_construct' ) {
        ch_rnaseq = Channel
            .fromList(samplesheetToList(input, schema_input))
            .map { meta, fastq_1, fastq_2 ->
                def reads = fastq_2 ? [ file(fastq_1), file(fastq_2) ] : [ file(fastq_1) ]
                [ meta.id, meta + [ single_end: !fastq_2 ], reads ]
            }
            .groupTuple()
            .flatMap { id, metas, reads_list ->
                def n = metas.size()

                // Condition guard: every rep under one sample id must declare the
                // same condition. Merge is keyed on `sample`, so two conditions
                // sharing a sample name would otherwise be pooled into one DB
                // silently — refuse it here instead.
                def conds = metas.collect { it.condition }.findAll { it != null }.unique()
                if ( conds.size() > 1 ) {
                    error "Sample '${id}' mixes conditions ${conds} across its reps. A sample id is a single condition — give each condition a distinct sample name (merge is keyed on the 'sample' column)."
                }

                // Read-group numbering: use the explicit `rep` column when every
                // rep of the sample provides it (deterministic rg<rep>); else fall
                // back to fastq-name order. rep must be unique within a sample.
                def haveRep = n > 0 && metas.every { it.rep != null }
                if ( haveRep && metas.collect { it.rep }.unique().size() != n ) {
                    error "Sample '${id}' has duplicate rep numbers ${metas.collect { it.rep }.sort()}; each rep within a sample must be unique."
                }
                def reps = (0..<n).collect { [ m: metas[it], reads: reads_list[it] ] }
                reps = haveRep ? reps.sort { it.m.rep } : reps.sort { it.reads[0].name }

                reps.withIndex().collect { rep, i ->
                    def k = haveRep ? rep.m.rep : (i + 1)
                    // Single-rep sample: leave meta untouched (id = sample), so
                    // existing samplesheets produce byte-identical metas. Only
                    // multi-rep samples gain the per-rep tags + unique id.
                    def m = n == 1
                        ? rep.m
                        : rep.m + [ sample: id, rep: k, num_reps: n, read_group: "rg${k}", id: "${id}_rg${k}" ]
                    [ m, rep.reads ]
                }
            }

        // Merge plan, logged once before any compute. Reps are grouped strictly
        // by the `sample` column: a shared name is pooled into one DB, distinct
        // names yield independent DBs. Printing the grouping makes a naming
        // mistake — e.g. two conditions sharing one sample id, silently merged
        // into a single DB — visible before a long run, complementing the
        // db_construct completeness guard.
        def lines    = file(input).readLines().findAll { it?.trim() }
        def header   = lines[0].split(',').collect { it.trim() }
        def iSample  = header.indexOf('sample')
        def iFq1     = header.indexOf('fastq_1')
        def iCond    = header.indexOf('condition')
        def dataRows = lines.drop(1)
        def repsBySample = [:]
        def condBySample = [:]
        dataRows.each { row ->
            def cols = row.split(',')
            def s = cols[iSample].trim()
            repsBySample.get(s, []) << (iFq1 >= 0 && cols.size() > iFq1 ? file(cols[iFq1].trim()).name : '')
            if ( iCond >= 0 && cols.size() > iCond ) { condBySample[s] = cols[iCond].trim() }
        }
        def planLines = repsBySample.collect { s, fqs ->
            def c = condBySample[s] ? " [condition=${condBySample[s]}]" : ''
            "    ${s}${c}  <- ${fqs.size()} rep(s)  ${fqs}"
        }.join('\n')
        log.info(
            "db_construct merge plan — ${dataRows.size()} rep-row(s) across ${repsBySample.size()} " +
            "sample(s) -> ${repsBySample.size()} cryptic DB(s):\n${planLines}\n" +
            "  reps are grouped by the 'sample' column; distinct names => independent DBs, " +
            "a shared name => reps pooled into one DB (read groups from the 'rep' column)."
        )
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
