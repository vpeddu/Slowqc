/*
 * Define the parameters used in the processes below
 */

/*
 * Define the processes used in this workflow
 */

process trim_files_SE {

    container "quay.io/biocontainers/bbmap:38.76--h516909a_0"

    input:
      file r1 

    output:
		file "${r1}"

    script:
      """
#!/bin/bash

set -e

# For logging and debugging, list all of the files in the working directory
ls -lahtr

# Get the sample name from the file name
sample_name=\$(echo ${r1} | sed 's/.R1.fastq.gz//')
echo "Processing \$sample_name"
# Rename the input file to make sure we don't use it as the output


echo "Masking ${r1}"
bbduk.sh \
    in=${r1} \
    out=${r1}.trimmed.fastq.gz \
    ref=${baseDir}/bin/adapters.fa 
    

    mv ${r1}.trimmed.fastq.gz ${r1}

"""
}

process trim_files_PE {

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
-Xmx4g \
    in=${r1} \
	in2=${r2} \
    out=${r1}.trimmed.fastq.gz \
	out2=${r2}.trimmed.fastq.gz \
    ref=${baseDir}/bin/adapters.fa 
    

    mv ${r1}.trimmed.fastq.gz ${r1}
	mv ${r2}.trimmed.fastq.gz ${r2}

"""
}
process kallisto_qc_SE { 
	container "quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1"

    input:
      file r1
	  val KALLISTO_ARGS

    output:
      file "*.QC"
 script:

      """
#!/bin/bash
	folderName=`basename ${r1} ".fastq.gz"`
	echo \$folderName
	kallisto quant --single -i ${baseDir}/bin/qc_idx ${KALLISTO_ARGS} -t ${task.cpus} -o \$folderName ${r1}

	# find and rename the kallisto files
	for i in `find . -name *.tsv`; do  temp=`echo \$i | cut -d / -f2`; newfilename="\$temp"".tsv";echo \$newfilename; echo \$temp; cp \$temp/abundance.tsv \$newfilename.QC; done


"""
}

process kallisto_qc_PE { 
	container "quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1"

    input:
      tuple file(r1), file(r2)

    output:
      file "*.QC"
 script:

      """
#!/bin/bash
	folderName=`basename ${r1} "_1.fastq.gz"`
	echo \$folderName
	kallisto quant -i ${baseDir}/bin/qc_idx -t ${task.cpus} -o \$folderName ${r1} ${r2}

	# find and rename the kallisto files
	for i in `find . -name *.tsv`; do  temp=`echo \$i | cut -d / -f2`; newfilename="\$temp"".tsv";echo \$newfilename; echo \$temp; cp \$temp/abundance.tsv \$newfilename.QC; done


"""
}


process kallisto_human_SE { 
	container "quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1"

    input:
      file r1
	  file kallistoIndex
	  val kallistoArgs

    output:
      file "${r1}" 
	  file "*.tsv"
 script:

      """
	#!/bin/bash
	folderName=`basename ${r1} ".fastq.gz"`
	echo \$folderName
	kallisto quant --single -i ${kallistoIndex} ${kallistoArgs} -t ${task.cpus} -o \$folderName ${r1}

	# find and rename the kallisto files
	for i in `find . -name *.tsv`; do  temp=`echo \$i | cut -d / -f2`; newfilename="\$temp"".tsv";echo \$newfilename; echo \$temp; cp \$temp/abundance.tsv \$newfilename; done
	"""
}

process kallisto_human_PE { 
	container "quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1"

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

	ls -lah

	Rscript --vanilla ${baseDir}/bin/generate_debrowser.r . ${baseDir}/bin/transcripts_to_genes.txt

"""
}

process Final_SE{ 
	container "quay.io/wshands/fastqc"

	input: 
		file r1 
		file qcfile
		file(DEBrowser_input)
		file(tsv_list)
		
	output: 
		file "sample_counts.csv"
		file "trimmed_files"
		file "fastqc_output"
		file "DEBrowser_input.txt"

	publishDir "${params.OUTDIR}"
	script:
	"""
	#!/bin/bash

	ls -lah 

	echo "Sample,Mitochondrial,Ribosomal,Human, Total" >> sample_counts.csv

	for i in *.tsv; do  temp=`basename \$i .tsv` ; temp_fastq="\$temp""_1.fastq.gz"
	echo -n "\$temp," >> sample_counts.csv
	echo -n `cut -f1,4  \$temp.tsv.QC  | grep MN046426.1 | cut -f2` >> sample_counts.csv ; echo -n , >> sample_counts.csv
	echo -n `cut -f1,4  \$temp.tsv.QC  | grep "NR"  | cut -f2 | awk '{s+=\$1} END {print s}'`>> sample_counts.csv ; echo -n , >> sample_counts.csv
	echo -n `cut -f4  \$temp.tsv  | awk '{s+=\$1} END {print s}'` >> sample_counts.csv ; echo -n , >> sample_counts.csv
	zgrep -c "^+\$" \$temp_fastq >> sample_counts.csv
	done

	
	mkdir fastqc_output
	for i in *.fastq.gz* ; do fastqc \$i ; done
	mv *_fastqc* fastqc_output


	mkdir trimmed_files
	mv *.gz trimmed_files

	"""

}

process Final_PE { 
	container "quay.io/wshands/fastqc"

	input: 
		tuple file(r1), file(r2)
		file qcfile
		file(DEBrowser_input)
		file(tsv_list)
		
	output: 
		file "sample_counts.csv"
		file "trimmed_files"
		file "fastqc_output"
		file "DEBrowser_input.txt"

	publishDir "${params.OUTDIR}"

	script:
	"""
	#!/bin/bash

	ls -lah 

	echo "Sample,Mitochondrial,Ribosomal,Human, Total" >> sample_counts.csv

	for i in *.tsv; do  temp=`basename \$i .tsv` ; temp_fastq="\$temp""_1.fastq.gz"
	echo -n "\$temp," >> sample_counts.csv
	echo -n `cut -f1,4  \$temp.tsv.QC  | grep MN046426.1 | cut -f2` >> sample_counts.csv ; echo -n , >> sample_counts.csv
	echo -n `cut -f1,4  \$temp.tsv.QC  | grep "NR"  | cut -f2 | awk '{s+=\$1} END {print s}'`>> sample_counts.csv ; echo -n , >> sample_counts.csv
	echo -n `cut -f4  \$temp.tsv  | awk '{s+=\$1} END {print s}'` >> sample_counts.csv ; echo -n , >> sample_counts.csv
	zgrep -c "^+\$" \$temp_fastq >> sample_counts.csv
	done

	
	mkdir fastqc_output
	for i in *_1.fastq.gz* ; do fastqc \$i ; done
	mv *_fastqc* fastqc_output


	mkdir trimmed_files
	mv *.gz trimmed_files

	"""

}
