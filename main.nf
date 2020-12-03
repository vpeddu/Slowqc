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
// single test command  nextflow run main.nf --INPUT /Users/gerbix/Documents/vikas/slowqc/test_data/paired/  --OUTDIR test/  --KALLISTO "/Users/gerbix/Documents/vikas/bin/grch38_transcriptome_kallisto_index/transcriptome.idx" -with-trace -with-docker ubuntu:18.04 --KALLISTO_ARGS "-l 300 -s 20" -resume 

// https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-96/homo_sapiens.tar.gz



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
        --TRANSCRIPTOME <transcriptome_location> Runs Kallisto with the specified transcriptome location and creates 
                                            DEBrowser ready input; Default: False 
                                            TODO put the kallisto transcriptome for the most recent 
                                            release up on the slowqc github page to avoid versioning fuckery
        --KALLISTO_ARGS Arguments that go after "kallisto quant". If you're running a single ended set 
                        -l and -s go here. Must be in double quotes. 
    
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

include trim_files_SE from './Modules.nf'
include trim_files_PE from './Modules.nf'

include kallisto_qc_SE from './Modules.nf'
include kallisto_qc_PE from './Modules.nf'

include kallisto_human_SE from './Modules.nf'
include kallisto_human_PE from './Modules.nf'

include DEbrowser from './Modules.nf'

include Final_SE from './Modules.nf'
include Final_PE from './Modules.nf'


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
            .fromFilePairs("${params.INPUT}*_R{1,2}*.gz")
            .ifEmpty { error "Cannot find any FASTQ pairs in ${params.INPUT} ending with .gz" }
            .map { it -> [it[0], it[1][0], it[1][1]]}
        } else { 
        input_read_ch = Channel
            .fromPath("${params.INPUT}*.gz")
            .map { it -> [ file(it)]}
        }
    if(params.PAIRED){ 
        trim_files_PE( 
            input_read_ch
        )

        kallisto_qc_PE( 
            trim_files_PE.out[0]
        )

        kallisto_human_PE(
            trim_files_PE.out[0],
            file(params.TRANSCRIPTOME),
            params.KALLISTO_ARGS
        )
        DEbrowser ( 
            kallisto_human_PE.out[1].collect()
        )

        Final_PE( 
            trim_files_PE.out.toList(), 
            kallisto_qc_PE.out[0].collect(),
            DEbrowser.out[0],
            kallisto_human_PE.out[1].collect()

        )

        
        // publish:
        //     Final_PE.out to: "${params.OUTDIR}" , mode: 'copy'
        //end PAIRED
    } 
   else{ 

        trim_files_SE( 
                input_read_ch
            )

            kallisto_qc_SE( 
                trim_files_SE.out[0],
                params.KALLISTO_ARGS
            )

            kallisto_human_SE(
                trim_files_SE.out[0],
                file(params.TRANSCRIPTOME),
                params.KALLISTO_ARGS
            )
            DEbrowser ( 
                kallisto_human_SE.out[1].collect()
            )

            Final_SE( 
                trim_files_SE.out.toList(), 
                kallisto_qc_SE.out[0].collect(),
                DEbrowser.out[0],
                kallisto_human_SE.out[1].collect()
            )
            // publish:
            //     Final_SE.out to: "${params.OUTDIR}" , mode: 'copy'

            
   }
}      

