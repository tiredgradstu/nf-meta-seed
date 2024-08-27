include { FASTP } from './modules.nf'
include { DOWNLOAD_KRAKEN_DB } from './modules.nf'
include { KRAKEN } from './modules.nf'
include { KRAKEN_DL_TAX } from './modules.nf'
include { KRAKEN_DBDL_LIB } from './modules.nf'
include { KRAKEN_DBDL_ADDGENO } from './modules.nf'
include { KRAKEN_DB_BUILD_BUILD } from './modules.nf'


workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    
    fastp_ch = FASTP(reads_ch)

    // Download Kraken db *if* no db is provided to the pipeline
    kraken_db_ch = params.kraken_db 
        // Turn into a value channel with .first() when we create teh channel from a file.
        ? Channel.fromPath(params.kraken_db, checkIfExists: true).first()
        : DOWNLOAD_KRAKEN_DB(params.kraken_db_url)

    KRAKEN(fastp_ch.fastq, kraken_db_ch)

    // Build a custom Kraken Db
    tax_ch = KRAKEN_DL_TAX()
    lib_ch = KRAKEN_DBDL_LIB(params.library)
    db_unbuilt_ch = KRAKEN_DBDL_ADDGENO(params.genome_dir, tax_ch, lib_ch)
    KRAKEN_DB_BUILD_BUILD(db_unbuilt_ch)

    //ASSEMBLY (fastp_ch.fastq)

}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/mockdata_preports\n" : "Oops .. something went wrong" )
}


#!/usr/bin/env nextflow

// Process parameters
krakendb_libs = params.krakendb_libs ? params.krakendb_libs?.split(',') as List : null
//krakendb_libs = params.krakendb_libs?.split(',') as List

// Report
log.info """
    M E T A G E N O M I C S - N F   P I P E L I N E
    ===============================================
    Reads                                           : ${params.reads}
    Outdir                                          : ${params.outdir}
    Path to pre-existing Kraken DB (if any)         : ${params.krakendb_dir}
    URL to pre-existing Kraken DB (if any)          : ${params.krakendb_url}
    Dir with genomes for custom Kraken DB (if any)  : ${params.krakendb_add}
    Library/ies for custom Kraken DB (if any)       : ${params.krakendb_libs}
    ===============================================
    """
    .stripIndent(true)

// Import processes
include { FASTP } from './modules.nf'
include { KRAKEN_DL_DB } from './modules.nf'
include { KRAKEN_DL_TAX } from './modules.nf'
include { KRAKEN_DL_LIB } from './modules.nf'
include { KRAKEN_ADD } from './modules.nf'
include { KRAKEN_BUILD } from './modules.nf'
include { KRAKEN_RUN } from './modules.nf'
include { ASSEMBLY } from './modules.nf'
include { MULTIQC } from './modules.nf'

// Define the workflow
workflow  {
    // Create channels from input files
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    krakendb_ch = params.krakendb_dir
        ? Channel.fromPath(params.krakendb_dir, checkIfExists: true).first()
        : null
    krakendb_add_ch = params.krakendb_add
        ? Channel.fromPath(params.krakendb_add, checkIfExists: true).first()
        : null
    krakendb_liblist_ch = krakendb_libs ? Channel.fromList(krakendb_libs) : null

    // Read preprocessing and QC
    fastp_ch = FASTP(reads_ch)
    
    // Kraken
    if (!krakendb_ch) {
        if (krakendb_add_ch && krakendb_liblist_ch) {
            println "\nNOTE: The pipeline will build and use a custom Kraken DB\n"
            krakendb_tax_ch = KRAKEN_DL_TAX()
            krakendb_lib_ch = KRAKEN_DL_LIB(krakendb_liblist_ch)
            krakendb_unbuilt_ch = KRAKEN_ADD(krakendb_add_ch, krakendb_tax_ch, krakendb_lib_ch)
            krakendb_ch = KRAKEN_BUILD(krakendb_unbuilt_ch)
        } else {
            println "\nNOTE: The pipeline will download and use an online Kraken DB\n"
            krakendb_ch = KRAKEN_DL_DB(params.krakendb_url)
        }
    } else {
        println "\nNOTE: The pipeline will use the provided Kraken DB ($params.krakendb_dir)\n"
    }
    kraken_ch = KRAKEN_RUN(fastp_ch.fastq, krakendb_ch)

    // Assembly and binning
    ASSEMBLY(fastp_ch.fastq)

    // MultiQC
    mqc_in_ch = kraken_ch.report.mix(fastp_ch.report).collect()
    MULTIQC(mqc_in_ch)
}

// Report
workflow.onComplete {
    if (workflow.success) {
        log.info ("\nThe pipeline has finished successfully! Final outputs are in the $params.outdir dir.")
    } else {
        log.info ("\nThe pipeline encountered an error and did not finish successfully")
    }
}
