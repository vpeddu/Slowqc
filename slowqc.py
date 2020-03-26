#!/bin/python

import os
import sys
import csv
import glob 
import time
import gzip
import pysam
import fnmatch
import logging
import argparse
import subprocess
import urllib.request
import multiprocessing
from Bio import Entrez
from functools import reduce
Entrez.email = "vpeddu@uw.edu"


#Create log file
logging.basicConfig(filename='slowqc.log',level=logging.DEBUG, format='%(asctime)s %(message)s')

#List all R1 fastq files in the current folder 
r1_fastqs = fnmatch.filter(os.listdir(), '*R1*')
r1_fastqs.sort()
r2_fastqs = fnmatch.filter(os.listdir(), '*R2*')
r2_fastqs.sort()

logging.info(r1_fastqs)
logging.info(r2_fastqs)

#Parse arguments 
parser = argparse.ArgumentParser(description = "input trimmomatic flags")
parser.add_argument("--paired", dest = 'paired', type = bool, help = "Option for paired end reads")
parser.add_argument("--spades", dest = 'spades' ,type=bool,  help = "Input all flags after SPADES.py enclosed in double quotes here")
parser.add_argument("--trimmomatic", dest='trimmomatic',nargs = '?' ,type=str, 
	help="Input flags after 'trimmomatic PE' enclosed in double quotes \
	excluding the filenames. Adapter file must be named adapters.fa. This will \
	automatically download from github")
parser.add_argument("--kallisto", dest = 'kallisto',type = str, help = "Input all flags after 'kallisto_quant' enclosed in double quotes here \
	IMPORTANT: if single end you must supply -l and -s arguments for kallisto here- the pipeline will break if you do not")
parser.add_argument("--debrowser", dest = 'debrowser', type = bool, help = "Takes kallisto output and rewrites as debrowser ready input file")
parser.add_argument("--fastqc", dest = 'fastqc', type = bool, help = "Runs fastqc on all trimmed files")

args = parser.parse_args()

#Run Trimmomatic
if(args.trimmomatic):
	#download adapters.fa from github
	url = "https://github.com/vpeddu/Slowqc/raw/master/bin/adapters.fa"
	filename, headers = urllib.request.urlretrieve(url, filename="adapters.fa")
	trimmed_read_count = []
	trimmed_read_count.append('Trimmed read count')
	print('Running Trimmomatic with arguments: ' + args.trimmomatic)
	for fq in range(len(r1_fastqs)): 
		fq_name = (r1_fastqs[fq].split('_R')[0])
		print('Trimming ' + fq_name)
		trim_r1_name = 'trimmed.' + r1_fastqs[fq]
		if(args.paired):
			trim_r2_name = 'trimmed.' + r2_fastqs[fq]
			trimmomatic_cmd = 'trimmomatic PE ' + ' ' + r1_fastqs[fq] + ' ' + r2_fastqs[fq] + ' ' + trim_r1_name + ' /dev/null/ ' + trim_r2_name + ' /dev/null/ ' + args.trimmomatic
		else:
			trimmomatic_cmd = 'trimmomatic SE ' + ' ' + r1_fastqs[fq] + ' ' + trim_r1_name  + ' ' + args.trimmomatic
		subprocess.call(trimmomatic_cmd, shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
		logging.info(trimmomatic_cmd)
	#Redefining fastq list so everything below works off of the trimmed files 
	trimmed_files = fnmatch.filter(os.listdir(), 'trimmed.*')
	r1_fastqs = fnmatch.filter(trimmed_files, '*R1*')
	r2_fastqs = fnmatch.filter(trimmed_files, '*R2*')


#Run Kallisto
if(args.kallisto):
	kallisto_mapped_read_count = []
	kallisto_mapped_read_count.append('Kallisto mapped read count')
	print('Running Kallisto with arguments: ' + args.kallisto)
	for fq in range(len(r1_fastqs)): 
		fq_name = (r1_fastqs[fq].split('_R')[0])
		fq_name = fq_name.split('trimmed.')[1]
		print('Running Kallisto for ' + fq_name)
		if args.paired: 
			kallisto_cmd = 'kallisto quant ' + args.kallisto + " -o " + fq_name + ' ' + r1_fastqs[fq] + ' ' + r2_fastqs[fq]
		else: 
			kallisto_cmd = 'kallisto quant --single ' + args.kallisto + " -o " + fq_name + ' ' + r1_fastqs[fq] 
		logging.info(kallisto_cmd)
		subprocess.call(kallisto_cmd, shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
		count_cmd = "awk '{total = total + $4}END{print total}' " + fq_name + '/abundance.tsv'
		kallisto_mapped_read_count.append((subprocess.check_output(count_cmd, shell = True).decode('utf-8').strip()))

#Run SPADES
if(args.spades):
	print('Running SPADES')
	for fq in range(len(r1_fastqs)): 
		fq_name = (r1_fastqs[fq].split('_R')[0])
		fq_name = fq_name.split('trimmed.')[1]
		print('Running SPADES for ' + fq_name)
		if args.paired: 
			spades_cmd = 'spades.py ' + " -1 " + r1_fastqs[fq] + ' -2 ' + r2_fastqs[fq] +" -o " + fq_name 
			subprocess.call(spades_cmd, shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
		else: 
			spades_cmd = 'spades.py ' + " --s1 " + r1_fastqs[fq] + " -o " + fq_name 
			subprocess.call(spades_cmd, shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
		logging.info(spades_cmd)

print('Counting reads in fastq files')

#Function to count lines in a file
def file_len_zipped(fname):
    with gzip.open(fname, 'rb') as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def file_len_unzipped(fname):
    with gzip.open(fname, 'rb') as f:
        for i, l in enumerate(f):
            pass
    return i + 1

read_counts=[]
read_counts.append('read counts')

#counts lines in fastq file and divides by 4 to get total read count
for fq in r1_fastqs: 
	if(fq.lower().endswith('.gz')):
		total_lines = file_len_zipped(fq)
	else:
		total_lines = file_len_unzipped(fq)
	if args.paired:
		reads = 2 * (total_lines / 4)
	if not args.paired:
		reads = (total_lines / 4)
	read_counts.append(reads)


cores_available = str(multiprocessing.cpu_count())
print('Counting mitochondrial reads')

#Mitochondrial reference genome. If you want to change this change the id argument to your favorite reference
mitochondria = Entrez.efetch(db="nucleotide", id="NC_012920.1", rettype="fasta", retmode="text")
mitochondria_file = open("mitochondria.fasta", "w")
mitochondria_file.write(mitochondria.read())
mitochondria_file.close()

#Builds Bowtie2 index for the mitochondrial reference
subprocess.call("bowtie2-build mitochondria.fasta mitochondria", shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)

mito_count=[]
mito_count.append('Mitochondrial reads')

#Runs Bowtie2 for the fastq against mitochondria
for fq in range(len(r1_fastqs)): 
	fq_name = (r1_fastqs[fq].split('_R')[0])
	out =  fq_name + '.mitochondria.sam'
	if args.paired: 
		cmd = 'bowtie2 -x mitochondria --no-unal -1 ' + r1_fastqs[fq] + ' -2 ' + r2_fastqs[fq] + ' -p ' + cores_available + ' -S ' + out
	else:
		cmd = 'bowtie2 -x mitochondria --no-unal -U ' + r1_fastqs[fq] +  ' -p ' + cores_available + ' -S ' + out
	logging.info(cmd)
	subprocess.call(cmd, shell = True, stderr = subprocess.DEVNULL)
	mito_count.append(int(subprocess.check_output('cat '+ out + ' | wc -l ', shell = True).decode('utf-8').strip())-3)

print('Counting 5.8s reads')

#18s reference genome. If you want to change this change the id argument to your favorite reference
ribosome_5 = Entrez.efetch(db="nucleotide", id="NR_146147.1", rettype="fasta", retmode="text")

ribosome_5_file = open("5s.fasta", "w")
ribosome_5_file.write(ribosome_5.read())
ribosome_5_file.close()

#Builds Bowtie2 index for the 18s reference
subprocess.call("bowtie2-build 5s.fasta 5s", shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)

ribosome_5_count=[]
ribosome_5_count.append('5s counts')

#Runs Bowtie2 for the fastq against 18s
for fq in range(len(r1_fastqs)): 
	fq_name = (r1_fastqs[fq].split('_R')[0])
	out =  fq_name + '.5s.sam'
	if args.paired:
		cmd = 'bowtie2 -x 5s --no-unal -1 ' + r1_fastqs[fq] + ' -2 ' + r2_fastqs[fq] + ' -p ' + cores_available + ' -S ' + out
	else: 
		cmd = 'bowtie2 -x 5s --no-unal -U ' + r1_fastqs[fq] +  ' -p ' + cores_available + ' -S ' + out
	logging.info(cmd)
	subprocess.call(cmd, shell = True, stderr = subprocess.DEVNULL)
	ribosome_5_count.append(int(subprocess.check_output('cat '+ out + ' | wc -l ', shell = True).decode('utf-8').strip())-3)

print('Counting 18s reads')

#18s reference genome. If you want to change this change the id argument to your favorite reference
ribosome_18 = Entrez.efetch(db="nucleotide", id="NR_146146.1", rettype="fasta", retmode="text")

ribosome_18_file = open("18s.fasta", "w")
ribosome_18_file.write(ribosome_18.read())
ribosome_18_file.close()

#Builds Bowtie2 index for the 18s reference
subprocess.call("bowtie2-build 18s.fasta 18s", shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)

ribosome_18_count=[]
ribosome_18_count.append('18s counts')

#Runs Bowtie2 for the fastq against 18s
for fq in range(len(r1_fastqs)): 
	fq_name = (r1_fastqs[fq].split('_R')[0])
	out =  fq_name + '.18s.sam'
	if args.paired:
		cmd = 'bowtie2 -x 18s --no-unal -1 ' + r1_fastqs[fq] + ' -2 ' + r2_fastqs[fq] + ' -p ' + cores_available + ' -S ' + out
	else:
		cmd = 'bowtie2 -x 18s --no-unal -U ' + r1_fastqs[fq] +  ' -p ' + cores_available + ' -S ' + out
	logging.info(cmd)
	subprocess.call(cmd, shell = True, stderr = subprocess.DEVNULL)
	ribosome_18_count.append(int(subprocess.check_output('cat '+ out + ' | wc -l ', shell = True).decode('utf-8').strip())-3)


print('Counting 28s reads')

#28s reference genome. If you want to change this change the id argument to your favorite reference
ribosome_28 = Entrez.efetch(db="nucleotide", id="NR_146148.1", rettype="fasta", retmode="text")
ribosome_28_file = open("28s.fasta", "w")
ribosome_28_file.write(ribosome_28.read())
ribosome_28_file.close()

#Builds Bowtie2 index for the 28s reference
subprocess.call("bowtie2-build 28s.fasta 28s", shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)

ribosome_28_count=[]
ribosome_28_count.append('28s counts')

#Runs Bowtie2 for the fastq against 28s
for fq in range(len(r1_fastqs)): 
	fq_name = (r1_fastqs[fq].split('_R')[0])
	out =  fq_name + '.28s.sam'
	if args.paired:
		cmd = 'bowtie2 -x 28s --no-unal -1 ' + r1_fastqs[fq] + ' -2 ' + r2_fastqs[fq] + ' -p ' + cores_available + ' -S ' + out
	else:
		cmd = 'bowtie2 -x 28s --no-unal -U ' + r1_fastqs[fq] +  ' -p ' + cores_available + ' -S ' + out
	logging.info(cmd)
	subprocess.call(cmd, shell = True, stderr = subprocess.DEVNULL)
	ribosome_28_count.append(int(subprocess.check_output('cat '+ out + ' | wc -l ', shell = True).decode('utf-8').strip())-3)

fq_column_header = 'fastq file'
fastqs =[fq_column_header] + r1_fastqs
if(args.kallisto):
	final = zip(fastqs, read_counts, kallisto_mapped_read_count ,mito_count, ribosome_5_count, ribosome_18_count, ribosome_28_count)
else:
	final = zip(fastqs, read_counts, mito_count, ribosome_5_count, ribosome_18_count, ribosome_28_count)

print('Writing read counts file')

#writing CSV file
with open('Read_counts.csv', "w") as f:
    writer = csv.writer(f)
    for row in final:
        writer.writerow(row)

	
#Tidy everything up  
subprocess.call('mkdir alignments', shell = True)
subprocess.call('mv *.sam alignments; mv *.bt2 alignments', shell = True)
if args.trimmomatic: 
	subprocess.call('mkdir trimmed_files; mv trimmed.* trimmed_files', shell = True)

#Runs kallisto_to_debrowser.r 
if args.debrowser:
	print('Creating debrowser formatted input')
	url = "https://github.com/vpeddu/Bioinformatics-scripts/raw/master/kallisto_to_debrowser.r"
	filename, headers = urllib.request.urlretrieve(url, filename="kallisto_to_debrowser.r")
	tsv_cmd = 'mkdir kallisto_tsvs ; for i in `find . -name *.tsv`; do  temp=`echo $i | cut -d / -f2`; newfilename="$temp"".tsv"; cp $temp/abundance.tsv kallisto_tsvs/$newfilename; done'
	subprocess.call(tsv_cmd, shell = True, stderr = subprocess.DEVNULL)
	subprocess.call('rscript --vanilla kallisto_to_debrowser.r kallisto_tsvs/', shell = True , stderr = subprocess.DEVNULL)
	subprocess.call("sed 's/\"//g' DEBrowser_input.txt", shell = True)




#Runs FASTQC
if args.fastqc:
	print('Running FASTQC') 
	subprocess.call('FASTQC *R1*', shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
	subprocess.call('mkdir fastqc_files', shell = True)
	subprocess.call('mv *fastqc* fastqc_files', shell = True, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
