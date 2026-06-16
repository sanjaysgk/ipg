//
// DENOVO: de novo peptide sequencing discovery lane (Phase 1a — InstaNovo + classify).
//
// De novo engines predict peptides directly from spectra (no database). Each
// prediction is then classified against the sample's cryptic ORF DB + the canonical
// proteome (canonical/cryptic/novel). Confidence is the raw de novo log-probability;
// calibrated Winnow FDR is a separate, later step. Runs in parallel to the
// database-search engines and never feeds them — the two lanes meet only at the report.
//

include { INSTANOVO_PREDICT } from '../../../modules/local/instanovo_predict/main'
include { DENOVO_CLASSIFY    } from '../../../modules/local/denovo_classify/main'

workflow DENOVO {

    take:
    ch_spectra_db       // channel: [ val(meta), path(spectra), path(cryptic_fasta) ]
    ch_canonical_fasta  // path:    canonical proteome FASTA (or [] to skip)

    main:
    ch_versions = Channel.empty()

    // 1. De novo sequencing (currently InstaNovo; Casanovo joins in Phase 2).
    INSTANOVO_PREDICT(ch_spectra_db.map { meta, spectra, _db -> [meta, spectra] })
    ch_versions = ch_versions.mix(INSTANOVO_PREDICT.out.versions)

    // 2. Classify each sample's predictions against its own cryptic DB (join by meta).
    ch_classify_in = INSTANOVO_PREDICT.out.predictions
        .join(ch_spectra_db.map { meta, _spectra, db -> [meta, db] }, by: 0)

    DENOVO_CLASSIFY(ch_classify_in, ch_canonical_fasta)
    ch_versions = ch_versions.mix(DENOVO_CLASSIFY.out.versions)

    emit:
    classified  = DENOVO_CLASSIFY.out.classified     // channel: [ val(meta), path(tsv) ]
    predictions = INSTANOVO_PREDICT.out.predictions   // channel: [ val(meta), path(csv) ]
    versions    = ch_versions                         // channel: path(versions.yml)
}
