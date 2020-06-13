#!/usr/bin/env nextflow
/*
========================================================================================
                         slowQC
========================================================================================
 Quality control for RNAseq 
 #### Homepage / Documentation
https://github.com/vpeddu/slowqc
----------------------------------------------------------------------------------------
*/

// Using the Nextflow DSL-2 to account for the logic flow of this workflow
//test command
//nextflow run main.nf --INPUT /Users/gerbix/Documents/vikas/slowqc/test_data/paired/  -with-docker ubuntu:18.04 --OUTDIR test/ --PAIRED --KALLISTO "/Users/gerbix/Documents/vikas/bin/grch38_transcriptome_kallisto_index/transcriptome.idx" -with-trace -resume 

nextflow.preview.dsl=2


def helpMessage() {
    log.info"""
    slowQC 

    Usage:

    An example command for running the pipeline is as follows:

    nextflow run FredHutch/slowqc \\
        --PAIRED_END \\
        --HOST_FILTER_FIRST \\
        --OUTDIR output/
        

        --INPUT Input folder
        --OUTDIR Output folder
        --PAIRED Samples are paired end. IMPORTANT: files must have _1 and _2 in them as delimiter; Default: False
        --KALLISTO <transcriptome_location> Runs Kallisto with the specified transcriptome location and creates 
                                            DEBrowser ready input; Default: False 
                                            TODO put the kallisto transcriptome for the most recent 
                                            release up on the slowqc github page to avoid versioning fuckery
        --KALLISTO_ARGS Arguments that go after "kallisto quant". If you're running a single ended set 
                        -l and -s go here
    
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}
params.PAIRED = false
params.KALLISTO = false
params.KALLISTO_ARGS = ""

/*
 * Import the processes used in this workflow
 */

include trim_files from './Modules.nf'
include kallisto_qc from './Modules.nf'
include kallisto_human from './Modules.nf'
include DEbrowser from './Modules.nf'



// Error handling for input flags
//if CONTROL_FASTQ not set 


// Make sure OUTDIR ends with trailing slash

if (!params.OUTDIR.endsWith("/")){
   params.OUTDIR = "${params.OUTDIR}/"
}






workflow {
        //fml() 
        if(params.PAIRED){ 
        input_read_ch = Channel
            .fromFilePairs("${params.INPUT}*_{1,2}*.gz")
            .ifEmpty { error "Cannot find any FASTQ pairs in ${params.INPUT} ending with .gz" }
            .map { it -> [it[0], it[1][0], it[1][1]]}
        } else { 
        input_read_ch = Channel
            .fromPath("${params.INPUT}*.gz")
            .map { it -> [ file(it)]}
}
    if(params.PAIRED){ 
    trim_files( 
        input_read_ch
    )

    kallisto_qc( 
        trim_files.out[0]
    )

    if(params.KALLISTO) { 
    kallisto_human(
        trim_files.out[0],
        file(params.KALLISTO),
        params.KALLISTO_ARGS
    )
    DEbrowser ( 
        kallisto_human.out[1].collect()

    )
    }
    //end PAIRED
    } 
    publish:
        DEbrowser.out to: "${params.OUTDIR}" , mode: 'copy'
}      

