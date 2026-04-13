#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sanjaysgk/ipg
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github  : https://github.com/sanjaysgk/ipg
    Authors : Sanjay SG Krishna, Kate Scull, Chen Li (Monash University)
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { IPG  } from './workflows/ipg'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_ipg_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_ipg_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_ipg_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Resolve reference paths from --genome config lookup (iGenomes + genomes.config).
    The elvis operator (?:) ensures explicit CLI params always win, e.g.
    --genome GRCh38 --star_index /my/star  uses the override for star_index only.
    Params that remain null after resolution are either built by PREPARE_GENOME
    (fasta_fai, fasta_dict, star_index) or cause a validation error (fasta, gtf, etc.).
*/

// Genome attribute resolution is handled in nextflow.config (not here)
// because Nextflow 25.x treats params as immutable after config loading.
// See the params {} block in nextflow.config for the resolution logic.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_IPG {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    IPG (
        samplesheet
    )
    emit:
    multiqc_report = IPG.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_IPG (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_IPG.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
