/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE_CRYPTIC — orthogonal spectrum-level validation of cryptic peptides
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    1. PEPQUERY_PREP      extract cryptic peptide list + canonical-only reference
    2. PEPQUERY           re-match each cryptic candidate against its spectrum,
        competing with the canonical proteome (± modifications)
    3. PEPQUERY_ANNOTATE  write pepquery_status onto the integrated peptide table

    Cryptic calls that do not win this competitive re-match are flagged
    'not_confident'; canonical peptides are not tested. Samples with no cryptic
    peptides pass through unchanged.
----------------------------------------------------------------------------------------
*/

include { PEPQUERY_PREP     } from '../../../modules/local/pepquery_prep/main'
include { PEPQUERY          } from '../../../modules/local/pepquery/main'
include { PEPQUERY_ANNOTATE } from '../../../modules/local/pepquery_annotate/main'

workflow VALIDATE_CRYPTIC {

    take:
    ch_peptides      // channel: [ val(meta), path(integrated_peptides.tsv) ]
    ch_mgf           // channel: [ val(meta), path(mgf_files) ] from CONVERT_MZML
    ch_search_db     // channel: [ val(meta), path(target search FASTA) ] per sample

    main:

    ch_versions = Channel.empty()

    //
    // STEP 1: cryptic peptide list + canonical-only reference DB. Each sample's
    // peptides are prepped against ITS own search db (joined by meta).
    //
    ch_prep_in = ch_peptides.join(ch_search_db)   // [meta, peptides, fasta]
    PEPQUERY_PREP(
        ch_prep_in.map { meta, pep, _fa -> [meta, pep] },
        ch_prep_in.map { _meta, _pep, fa -> fa }
    )
    ch_versions = ch_versions.mix(PEPQUERY_PREP.out.versions)

    // Only run PepQuery where cryptic peptides actually exist.
    ch_have_cryptic = PEPQUERY_PREP.out.peptides.filter { _meta, f -> f.size() > 0 }

    //
    // STEP 2: PepQuery2 competitive re-match (peptides + spectra + ref, per sample).
    //
    ch_pepquery_in = ch_have_cryptic
        .join(PEPQUERY_PREP.out.ref_fasta)   // [meta, peptides, ref]
        .join(ch_mgf)                         // [meta, peptides, ref, mgf]
        .map { meta, peps, ref, mgf -> [meta, peps, mgf, ref] }

    PEPQUERY(ch_pepquery_in)
    ch_versions = ch_versions.mix(PEPQUERY.out.versions)

    //
    // STEP 3: write pepquery_status onto the integrated peptide table.
    //
    ch_annotate_in = ch_peptides.join(PEPQUERY.out.psm_rank)   // [meta, peptides_tsv, psm_rank]
    PEPQUERY_ANNOTATE(ch_annotate_in)
    ch_versions = ch_versions.mix(PEPQUERY_ANNOTATE.out.versions)

    // Samples that produced no cryptic peptides (so PepQuery never ran) keep
    // their original, unannotated table.
    ch_no_cryptic = ch_peptides
        .join(PEPQUERY_ANNOTATE.out.peptides, remainder: true)
        .filter { it -> it.size() > 2 && it[2] == null }
        .map { meta, orig, _validated -> [meta, orig] }

    ch_peptides_out = PEPQUERY_ANNOTATE.out.peptides.mix(ch_no_cryptic)

    emit:
    peptides = ch_peptides_out             // [meta, integrated_peptides(_validated).tsv]
    psm_rank = PEPQUERY.out.psm_rank       // [meta, psm_rank.txt]
    versions = ch_versions
}
