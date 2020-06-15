# Slowqc
RNAseq Pipeline to qc RNAseq FASTQ files with options for trimming, assembly, and export to iDEP/DEBrowser ready input. 

This script accepts untrimmed FASTQ files and performs the following: 
* Quality and adapter trimming using `BBDuk` 
* Transcript abundance quantification using `Kallisto`
* Creats `DEBrowser` ready input using the estimated counts from `Kallisto`
* Alignment and quantification to Human mitochondrial, and ribosomal sequences (5.8s, 18s, 28s) 
* Generate `FASTQC` output 
* Returns a csv with read counts of the original FASTQ file, reads aligning to the ribosome and mitochondrial sequences, and mapped reads to the Human transcriptome. 

## Running the pipeline 

This script is written in `NextFlow` and uses `Docker` containers to simplify workflow execution. To run, `Nextflow` and `Docker` must be installed. At runtime `Docker` must be running. 

### Flags: 
* --paired <True/False> Default: Single
* --trimmomatic "< Any Trimmomatic options to be entered after the 'trimmomatic PE' command>" 
* --kallisto "<Any Kallisto options to be entered following 'kallisto quant'>" 
* --Debrowser <True/False>
* --fastqc <True/False> 
* --spades <True/False> 

### Example command: 
#### Single end: 
```
nextflow run main.nf --INPUT <input_folder>  --OUTDIR test/  --TRANSCRIPTOME <path_to_kallisto_transcriptome> -with-trace -with-docker ubuntu:18.04 --KALLISTO_ARGS "-l <number> -s <number>" -resume 
```
#### Paired end: 
```
nextflow run main.nf --INPUT <input_folder>  --OUTDIR test/  --TRANSCRIPTOME <path_to_kallisto_transcriptome> -with-trace -with-docker ubuntu:18.04 --PAIRED -resume 
```