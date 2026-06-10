/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE_GENOME — resolve, decompress, and index reference files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Conditionally builds missing derivable indices (fai, dict, STAR index)
    from the primary reference files (fasta, gtf). Pass-through for files
    that are already provided via params or resolved from --genome config.

    Three usage modes:
        1. --genome GRCh38              resolves all paths from genomes.config
        2. --fasta /path --gtf /path    explicit paths, builds missing indices
        3. --genome GRCh38 --star_index /my/star   mix of config + overrides
----------------------------------------------------------------------------------------
*/

include { GUNZIP as GUNZIP_FASTA } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF   } from '../../../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX          } from '../../../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { STAR_GENOMEGENERATE     } from '../../../modules/nf-core/star/genomegenerate/main'
include { TABIX_TABIX as TABIX_DBSNP            } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_KNOWN_INDELS     } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_MILLS            } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GERMLINE_RESOURCE } from '../../../modules/nf-core/tabix/tabix/main'

/*
 * Resolve a reference parameter: explicit CLI param wins, then --genome
 * config lookup, then null. This replaces the getGenomeAttribute() approach
 * that relied on mutable params (broken in Nextflow 25.x).
 */
def resolve(String attr) {
    if ( params[attr] ) return params[attr]
    if ( params.genome && params.genomes && params.genomes.containsKey(params.genome) ) {
        def g = params.genomes[ params.genome ]
        if ( g.containsKey(attr) ) return g[ attr ]
        if ( attr == 'star_index' && g.containsKey('star') ) return g[ 'star' ]
    }
    return null
}

workflow PREPARE_GENOME {

    main:

    ch_versions = Channel.empty()
    def ref_meta = [ id: 'reference' ]

    // ---- Resolve all references (explicit params > genome config > null) ---------
    def r_fasta                 = resolve('fasta')
    def r_fasta_fai             = resolve('fasta_fai')
    def r_fasta_dict            = resolve('fasta_dict')
    def r_gtf                   = resolve('gtf')
    def r_star_index            = resolve('star_index')
    def r_rseqc_bed             = resolve('rseqc_bed')
    def r_dbsnp                 = resolve('dbsnp')
    def r_dbsnp_tbi             = resolve('dbsnp_tbi')
    def r_known_indels          = resolve('known_indels')
    def r_known_indels_tbi      = resolve('known_indels_tbi')
    def r_mills                 = resolve('mills')
    def r_mills_tbi             = resolve('mills_tbi')
    def r_germline_resource     = resolve('germline_resource')
    def r_germline_resource_tbi = resolve('germline_resource_tbi')

    // ---- Validate primary inputs ------------------------------------------------
    if ( !r_fasta ) {
        error "Reference FASTA is required. Provide --fasta or set --genome (e.g. --genome GRCh38)."
    }
    if ( !r_gtf ) {
        error "Gene annotation GTF is required. Provide --gtf or set --genome (e.g. --genome GRCh38)."
    }

    // ---- FASTA: decompress if .gz -----------------------------------------------
    if ( r_fasta.toString().endsWith('.gz') ) {
        GUNZIP_FASTA( [ [id: 'genome_fasta'], file(r_fasta) ] )
        ch_fasta = GUNZIP_FASTA.out.gunzip.map { meta, fasta -> [ ref_meta, fasta ] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions_gunzip)
    } else {
        ch_fasta = Channel.value( [ ref_meta, file(r_fasta) ] )
    }

    // ---- GTF: decompress if .gz -------------------------------------------------
    if ( r_gtf.toString().endsWith('.gz') ) {
        GUNZIP_GTF( [ [id: 'genome_gtf'], file(r_gtf) ] )
        ch_gtf = GUNZIP_GTF.out.gunzip.map { meta, gtf -> [ ref_meta, gtf ] }
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions_gunzip)
    } else {
        ch_gtf = Channel.value( [ ref_meta, file(r_gtf) ] )
    }

    // ---- FASTA_FAI: build if not provided ---------------------------------------
    if ( r_fasta_fai ) {
        ch_fai = Channel.value( [ ref_meta, file(r_fasta_fai) ] )
    } else {
        SAMTOOLS_FAIDX( ch_fasta.map { meta, fa -> [ meta, fa, [] ] }, false )
        ch_fai = SAMTOOLS_FAIDX.out.fai
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions_samtools)
    }

    // ---- FASTA_DICT: build if not provided --------------------------------------
    if ( r_fasta_dict ) {
        ch_dict = Channel.value( [ ref_meta, file(r_fasta_dict) ] )
    } else {
        GATK4_CREATESEQUENCEDICTIONARY( ch_fasta )
        ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions_gatk4)
    }

    // ---- STAR_INDEX: build if not provided (expensive: ~32 GB RAM, 30-60 min) ---
    if ( r_star_index ) {
        ch_star_index = Channel.value( [ ref_meta, file(r_star_index) ] )
    } else {
        STAR_GENOMEGENERATE( ch_fasta, ch_gtf )
        ch_star_index = STAR_GENOMEGENERATE.out.index
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions_star)
    }

    // ---- RSeQC BED: pass-through (cannot be auto-generated) ---------------------
    if ( !r_rseqc_bed ) {
        error "RSeQC BED file is required. Provide --rseqc_bed or set --genome (e.g. --genome GRCh38)."
    }
    ch_rseqc_bed = Channel.value( file(r_rseqc_bed) )

    // ---- GATK variant calling resources -----------------------------------------
    // Each VCF is required; the matching .tbi is auto-built via TABIX_TABIX if
    // not supplied. Accepts both bundled .tbi (fast) and missing-tbi (one-time
    // build, ~30s per index).
    if ( !r_dbsnp ) {
        error "dbSNP VCF is required. Provide --dbsnp or set --genome."
    }
    if ( !r_known_indels ) {
        error "Known indels VCF is required. Provide --known_indels or set --genome."
    }
    if ( !r_mills ) {
        error "Mills VCF is required. Provide --mills or set --genome."
    }
    if ( !r_germline_resource ) {
        error "Germline resource VCF is required. Provide --germline_resource or set --genome."
    }

    // Helper: emit [meta, vcf] for TABIX input shape
    def tabix_in = { name, path -> [ [id: name], file(path), [], [] ] }

    if ( r_dbsnp_tbi ) {
        ch_dbsnp_tbi = Channel.value( file(r_dbsnp_tbi) )
    } else {
        TABIX_DBSNP( Channel.value( tabix_in.call('dbsnp', r_dbsnp) ) )
        ch_dbsnp_tbi = TABIX_DBSNP.out.index.map { _meta, tbi -> tbi }
    }

    if ( r_known_indels_tbi ) {
        ch_known_indels_tbi = Channel.value( file(r_known_indels_tbi) )
    } else {
        TABIX_KNOWN_INDELS( Channel.value( tabix_in.call('known_indels', r_known_indels) ) )
        ch_known_indels_tbi = TABIX_KNOWN_INDELS.out.index.map { _meta, tbi -> tbi }
    }

    if ( r_mills_tbi ) {
        ch_mills_tbi = Channel.value( file(r_mills_tbi) )
    } else {
        TABIX_MILLS( Channel.value( tabix_in.call('mills', r_mills) ) )
        ch_mills_tbi = TABIX_MILLS.out.index.map { _meta, tbi -> tbi }
    }

    if ( r_germline_resource_tbi ) {
        ch_germline_resource_tbi = Channel.value( file(r_germline_resource_tbi) )
    } else {
        TABIX_GERMLINE_RESOURCE( Channel.value( tabix_in.call('germline_resource', r_germline_resource) ) )
        ch_germline_resource_tbi = TABIX_GERMLINE_RESOURCE.out.index.map { _meta, tbi -> tbi }
    }

    // known_sites is a list-channel [meta, [3 VCFs], [3 TBIs]] used by BaseRecalibrator
    ch_known_sites = Channel.value( [ ref_meta, [
        file(r_dbsnp),
        file(r_known_indels),
        file(r_mills),
    ] ] )
    // Combine .tbi paths (which may be original files or TABIX outputs) into a list channel
    ch_known_sites_tbi = ch_dbsnp_tbi
        .combine(ch_known_indels_tbi)
        .combine(ch_mills_tbi)
        .map { dbsnp_tbi, known_tbi, mills_tbi ->
            [ ref_meta, [ dbsnp_tbi, known_tbi, mills_tbi ] ]
        }
    ch_germline_resource     = Channel.value( file(r_germline_resource) )

    // ---- Composite channel: [meta, fasta, fai] for processes needing both -------
    ch_fasta_fai = ch_fasta.combine( ch_fai ).map { fasta_meta, fasta, fai_meta, fai ->
        [ fasta_meta, fasta, fai ]
    }

    emit:
    fasta                = ch_fasta
    fasta_fai            = ch_fai
    fasta_dict           = ch_dict
    fasta_fai_combo      = ch_fasta_fai
    star_index           = ch_star_index
    gtf                  = ch_gtf
    rseqc_bed            = ch_rseqc_bed
    known_sites          = ch_known_sites
    known_sites_tbi      = ch_known_sites_tbi
    germline_resource     = ch_germline_resource
    germline_resource_tbi = ch_germline_resource_tbi
    versions             = ch_versions
}
