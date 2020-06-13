/*
 * Define the parameters used in the processes below
 */

/*
 * Define the processes used in this workflow
 */


process trim_files {

    container "quay.io/biocontainers/bbmap:38.76--h516909a_0"

    input:
      tuple val(prefix), file(r1), file(r2)

    output:
      tuple file("${r1}"), file("${r2}")

    script:
      """
#!/bin/bash

set -e

# For logging and debugging, list all of the files in the working directory
ls -lahtr

# Get the sample name from the file name
sample_name=\$(echo ${r1} | sed 's/.R1.fastq.gz//')
echo "Processing \$sample_name"
# head ${r1} 
# echo ${r2}
# Rename the input file to make sure we don't use it as the output


echo "Masking ${r1}"
bbduk.sh \
    in=${r1} \
	in2=${r2} \
    out=${r1}.trimmed.fastq.gz \
	out2=${r2}.trimmed.fastq.gz \
    ref=${baseDir}/bin/adapters.fa 
    

    mv ${r1}.trimmed.fastq.gz ${r1}
	mv ${r2}.trimmed.fastq.gz ${r2}

"""
}


process kallisto_qc { 
	container "quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1"

    input:
      tuple file(r1), file(r2)

    output:
      tuple file("${r1}"), file("${r2}"), file("*.QC")
 script:

      """
#!/bin/bash
	folderName=`basename ${r1} "_1.fastq.gz"`
	echo \$folderName
	kallisto quant -i ${baseDir}/bin/qc_idx.idx -t ${task.cpus} -o \$folderName ${r1} ${r2}

	# find and rename the kallisto files
	for i in `find . -name *.tsv`; do  temp=`echo \$i | cut -d / -f2`; newfilename="\$temp"".tsv";echo \$newfilename; echo \$temp; cp \$temp/abundance.tsv \$newfilename.QC; done


"""
}


process kallisto_human { 
	//container "quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1"

    input:
      tuple file(r1), file(r2)
	  file kallistoIndex
	  val kallistoArgs

    output:
      tuple file("${r1}"), file("${r2}")
	  file "*.tsv"
 script:

      """
#!/bin/bash
	folderName=`basename ${r1} "_1.fastq.gz"`
	echo \$folderName
	kallisto quant -i ${kallistoIndex} ${kallistoArgs} -t ${task.cpus} -o \$folderName ${r1} ${r2}

	# find and rename the kallisto files
	for i in `find . -name *.tsv`; do  temp=`echo \$i | cut -d / -f2`; newfilename="\$temp"".tsv";echo \$newfilename; echo \$temp; cp \$temp/abundance.tsv \$newfilename; done
"""
}

process DEbrowser { 
    container "quay.io/vpeddu/rgeneratesummary:latest"

    input:
    file(tsv_list)

    output:
      file "DEBrowser_input.txt"
 script:

      """
	#!/bin/bash
	rscript --vanilla ${baseDir}/bin/generate_debrowser.r . ${baseDir}/bin/transcripts_to_genes.txt

"""
}
