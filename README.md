# Slowqc
RNAseq Pipeline to qc RNAseq FASTQ files with options for trimming, assembly, and export to iDEP/DEBrowser ready input. 

This script accepts untrimmed FASTQ files and performs the following: 
* Quality and adapter trimming using `Trimmomatic` 
* Transcript abundance quantification using Kallisto 
* Alignment and quantification to Human mitochondrial, and ribosomal sequences (5.8s, 18s, 28s) using  `Bowtie2` 
* Generate FASTQC output 
* Returns a csv with read counts of the original FASTQ file, reads aligning to the ribosome and mitochondrial sequences, and mapped reads to the Human transcriptome. 

## Running the pipeline 

All paired end FASTQ files in the folder with the python script will be run through the pipeline. Wildcard searches are used to distinguish R1 and R2 files, so they must have "R1" and "R2" in their filenames. `Trimmomatic`, `Bowtie2`, `Kallisto`, and `SPADES` must be installed and added to `$PATH`. The python script pulls from NCBI and the bin folder located inside this repo. As a result internet access is required to run this script. 

#### *If all of this is fine all you have to do is place the python script in the folder of FASTQ files* 

The required reference sequences for the `Bowtie2`  alignment and adapters.fa for the trimming will be automatically downloaded and built by the script. The only file that must be provided is for the optional kallisto alignment: [Transcriptome.idx](https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-96/homo_sapiens.tar.gz). 

### Flags: 
* --paired <True/False> Default: Single
* --trimmomatic "< Any Trimmomatic options to be entered after the 'trimmomatic PE' command>" 
* --kallisto "<Any Kallisto options to be entered following 'kallisto quant'>" 
* --Debrowser <True/False>
* --fastqc <True/False> 
* --spades <True/False> 

### Example command: 
`python slowqc.py --Paired True --trimmomatic "-phred33 -threads 40 ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75" --kallisto "-i ~/bin/homo_sapiens/transcriptome.idx --threads 40" --debrowser True --spades True` 
