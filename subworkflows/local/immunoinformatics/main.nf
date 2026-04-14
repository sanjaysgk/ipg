/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMMUNOINFORMATICS — downstream analysis of integrated peptides
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Consumes INTEGRATE_ENGINES.out.peptides and runs:
      - NETMHCPAN / NETMHCIIPAN  HLA class I / II binding (academic tools)
      - GIBBSCLUSTER             motif discovery
      - FLASHLFQ                 label-free quantification
      - BLASTP_HOST              non-host contamination check

    Each tool is individually gated by a --run_* param so users can pick
    exactly the subset they are licensed or configured for.
----------------------------------------------------------------------------------------
*/

include { NETMHCPAN                } from '../../../modules/local/netmhcpan/main'
include { NETMHCIIPAN              } from '../../../modules/local/netmhciipan/main'
include { GIBBSCLUSTER             } from '../../../modules/local/gibbscluster/main'
include { FLASHLFQ                 } from '../../../modules/local/flashlfq/main'
include { BLASTP_HOST              } from '../../../modules/local/blastp_host/main'
include { IMMUNOINFORMATICS_REPORT } from '../../../modules/local/immunoinformatics_report/main'

workflow IMMUNOINFORMATICS {

    take:
    ch_peptides    // [meta, integrated_peptides.tsv]
    ch_ms_files    // [meta, [ms files...]]  — for FlashLFQ

    main:

    ch_versions = Channel.empty()

    ch_netmhcpan_best    = Channel.empty()
    ch_netmhciipan_best  = Channel.empty()
    ch_gibbs             = Channel.empty()
    ch_flashlfq          = Channel.empty()
    ch_blastp            = Channel.empty()

    if (params.run_netmhcpan) {
        if (!params.netmhcpan_path) error("--run_netmhcpan requires --netmhcpan_path")
        if (!params.hla)            error("--run_netmhcpan requires --hla")
        NETMHCPAN(
            ch_peptides,
            params.hla,
            file(params.netmhcpan_path, checkIfExists: true)
        )
        ch_versions       = ch_versions.mix(NETMHCPAN.out.versions)
        ch_netmhcpan_best = NETMHCPAN.out.best
    }

    if (params.run_netmhciipan) {
        if (!params.netmhciipan_path) error("--run_netmhciipan requires --netmhciipan_path")
        if (!params.hla)              error("--run_netmhciipan requires --hla")
        NETMHCIIPAN(
            ch_peptides,
            params.hla,
            file(params.netmhciipan_path, checkIfExists: true)
        )
        ch_versions         = ch_versions.mix(NETMHCIIPAN.out.versions)
        ch_netmhciipan_best = NETMHCIIPAN.out.best
    }

    if (params.run_gibbscluster) {
        if (!params.gibbscluster_path) error("--run_gibbscluster requires --gibbscluster_path")
        GIBBSCLUSTER(
            ch_peptides,
            file(params.gibbscluster_path, checkIfExists: true)
        )
        ch_versions = ch_versions.mix(GIBBSCLUSTER.out.versions)
        ch_gibbs    = GIBBSCLUSTER.out.clusters
    }

    if (params.run_flashlfq) {
        // Pair each sample's peptides table with the MS files for that sample.
        ch_flashlfq_in = ch_peptides.join(ch_ms_files)
            .map { meta, pep, ms -> [meta, pep, ms] }
        FLASHLFQ(ch_flashlfq_in)
        ch_versions = ch_versions.mix(FLASHLFQ.out.versions)
        ch_flashlfq = FLASHLFQ.out.peptides
    }

    if (params.run_blastp_host) {
        if (!params.blast_db) error("--run_blastp_host requires --blast_db")
        // Stage every sibling of the db prefix — BLAST DBs are multi-file.
        def db_prefix = file(params.blast_db).name
        ch_blast_files = Channel.fromPath("${params.blast_db}*", checkIfExists: true).collect()
        BLASTP_HOST(
            ch_peptides,
            ch_blast_files,
            db_prefix,
            params.host_species
        )
        ch_versions = ch_versions.mix(BLASTP_HOST.out.versions)
        ch_blastp   = BLASTP_HOST.out.peptides
    }

    //
    // Final HTML report per sample. Pads missing inputs with the assets/NO_FILE
    // sentinel so every branch produces a report even when some tools are off.
    //
    def no_file = file("${projectDir}/assets/NO_FILE")
    ch_report_in = ch_peptides
        .map { meta, pep -> [meta, pep] }
        .join(ch_netmhcpan_best.map   { meta, f -> [meta, f] }.ifEmpty { [null, no_file] }, remainder: true)
        .map   { meta, pep, nmhc1 -> [meta, pep, nmhc1 ?: no_file] }
        .join(ch_netmhciipan_best.map { meta, f -> [meta, f] }.ifEmpty { [null, no_file] }, remainder: true)
        .map   { meta, pep, nmhc1, nmhc2 -> [meta, pep, nmhc1, nmhc2 ?: no_file] }
        .join(ch_gibbs.map            { meta, f -> [meta, f] }.ifEmpty { [null, no_file] }, remainder: true)
        .map   { meta, pep, nmhc1, nmhc2, gbs -> [meta, pep, nmhc1, nmhc2, gbs ?: no_file] }
        .join(ch_flashlfq.map         { meta, f -> [meta, f] }.ifEmpty { [null, no_file] }, remainder: true)
        .map   { meta, pep, nmhc1, nmhc2, gbs, lfq -> [meta, pep, nmhc1, nmhc2, gbs, lfq ?: no_file] }
        .join(ch_blastp.map           { meta, f -> [meta, f] }.ifEmpty { [null, no_file] }, remainder: true)
        .map   { meta, pep, nmhc1, nmhc2, gbs, lfq, bp -> [meta, pep, nmhc1, nmhc2, gbs, lfq, bp ?: no_file] }

    IMMUNOINFORMATICS_REPORT(ch_report_in)
    ch_versions = ch_versions.mix(IMMUNOINFORMATICS_REPORT.out.versions)

    emit:
    netmhcpan_best    = ch_netmhcpan_best
    netmhciipan_best  = ch_netmhciipan_best
    gibbs_clusters    = ch_gibbs
    flashlfq_peptides = ch_flashlfq
    blastp_peptides   = ch_blastp
    report            = IMMUNOINFORMATICS_REPORT.out.html
    versions          = ch_versions
}
