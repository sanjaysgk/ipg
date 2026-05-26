//
// VCF_ANNOTATE_ALL: annotate Mutect2 calls with VEP + SnpEff (+ optional AlphaMissense plugin).
//
// Leaf branch — outputs are published artifacts only and do not feed the
// cryptic-peptide pipeline. Engine selection via `tools` (comma list).
//

include { ENSEMBLVEP_VEP      } from '../../../modules/nf-core/ensemblvep/vep/main'
include { ENSEMBLVEP_DOWNLOAD } from '../../../modules/nf-core/ensemblvep/download/main'
include { SNPEFF_SNPEFF       } from '../../../modules/nf-core/snpeff/snpeff/main'
include { SNPEFF_DOWNLOAD     } from '../../../modules/nf-core/snpeff/download/main'
include { TABIX_BGZIPTABIX    } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow VCF_ANNOTATE_ALL {

    take:
    ch_vcf                // channel: [ val(meta), path(vcf), path(tbi) ]
    ch_fasta              // value:   [ val(meta), path(fasta) ]
    ch_fai                // value:   [ val(meta), path(fai) ]
    ch_vep_cache          // value:   path(cache) | []
    val_vep_genome        // val:     'GRCh38'
    val_vep_species       // val:     'homo_sapiens'
    val_vep_cache_version // val:     113
    ch_snpeff_cache       // value:   path(cache) | []
    val_snpeff_db         // val:     'GRCh38.105'
    ch_alphamissense_tsv  // value:   path(tsv.gz) [+ .tbi] | []
    val_tools             // val:     'vep,snpeff,alphamissense'

    main:

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()
    def tool_list = val_tools ? val_tools.split(',').collect { it.trim() } : []

    //
    // VEP (with optional AlphaMissense plugin)
    //
    ch_vep_vcf = Channel.empty()
    if (tool_list.contains('vep')) {
        def use_am = tool_list.contains('alphamissense') && ch_alphamissense_tsv

        // Resolve cache to a bare-path channel: use provided cache, else download once.
        // ENSEMBLVEP_VEP takes a bare `path cache`; ENSEMBLVEP_DOWNLOAD emits [meta, path].
        def ch_vep_cache_resolved
        if (ch_vep_cache) {
            ch_vep_cache_resolved = Channel.value(ch_vep_cache)
        } else {
            ENSEMBLVEP_DOWNLOAD(
                Channel.value([[id:'vep_cache'], val_vep_genome, val_vep_species, val_vep_cache_version])
            )
            ch_vep_cache_resolved = ENSEMBLVEP_DOWNLOAD.out.cache.map { _meta, cache -> cache }
            ch_versions = ch_versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)
        }

        ENSEMBLVEP_VEP(
            ch_vcf.map { meta, vcf, _tbi -> [meta, vcf, []] },
            val_vep_genome,
            val_vep_species,
            val_vep_cache_version,
            ch_vep_cache_resolved,
            ch_fasta,
            use_am ? ch_alphamissense_tsv : []
        )
        ch_vep_vcf  = ENSEMBLVEP_VEP.out.vcf
        ch_reports  = ch_reports.mix(ENSEMBLVEP_VEP.out.report)
        ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)
    }

    //
    // SnpEff
    //
    ch_snpeff_vcf = Channel.empty()
    if (tool_list.contains('snpeff')) {
        // SNPEFF_SNPEFF takes a [meta, cache] tuple; SNPEFF_DOWNLOAD emits the same.
        def ch_snpeff_cache_resolved
        if (ch_snpeff_cache) {
            ch_snpeff_cache_resolved = Channel.value([[id:'snpeff_cache'], ch_snpeff_cache])
        } else {
            SNPEFF_DOWNLOAD(
                Channel.value([[id:'snpeff_cache'], val_snpeff_db])
            )
            ch_snpeff_cache_resolved = SNPEFF_DOWNLOAD.out.cache
            ch_versions = ch_versions.mix(SNPEFF_DOWNLOAD.out.versions)
        }

        SNPEFF_SNPEFF(
            ch_vcf.map { meta, vcf, _tbi -> [meta, vcf] },
            val_snpeff_db,
            ch_snpeff_cache_resolved
        )
        ch_snpeff_vcf = SNPEFF_SNPEFF.out.vcf
        ch_reports    = ch_reports.mix(SNPEFF_SNPEFF.out.report)
        ch_versions   = ch_versions.mix(SNPEFF_SNPEFF.out.versions)
    }

    //
    // bgzip + tabix every annotated VCF
    //
    TABIX_BGZIPTABIX(ch_vep_vcf.mix(ch_snpeff_vcf))
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    emit:
    vcf_vep     = ch_vep_vcf
    vcf_snpeff  = ch_snpeff_vcf
    vcf_indexed = TABIX_BGZIPTABIX.out.gz_tbi
    reports     = ch_reports
    versions    = ch_versions
}
