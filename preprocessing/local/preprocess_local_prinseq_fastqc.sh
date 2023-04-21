#!/bin/bash

# load fastqc
module load fastqc/0.11.5

# Get the system time
now="$(date +"%c")"

# set variable
SRADIR="/lustre/alice3/scratch/spectre/(UID, e.g. y/yyl23)/sradata"

# Read accession numbers from the file and store in an array
mapfile -t ACCESSIONS < $SRADIR/acc_list.txt

# Unzip 
tar -xvzf prinseq-lite-0.20.4.tar.gz

# Create a new directory to pre-process dataset and copy the fastq files into it
mkdir -p $SRADIR/preprocess
echo "Copying fastqdump to new location. Please wait."
cp ./fastqdump/*fastq ./preprocess

cd $SRADIR/preprocess

# Loop through each accession number and download the corresponding file
for ACCESSION in "${ACCESSIONS[@]}"
do
  echo "PrinSeq-ing: $ACCESSION"
  perl ../prinseq-lite-0.20.4/prinseq-lite.pl -verbose -fastq "${ACCESSION}_1.fastq" -fastq2 "${ACCESSION}_2.fastq" -out_format 3 -out_bad null -min_len 30 -trim_left 10 -trim_qual_right 25 -lc_method entropy -lc_threshold 65
done
echo "PrinSeq Completed" 

# testing fastqc operativity
# fastqc --help

# make a directory to store fastqc reports
mkdir -p $SRADIR/fastqcoutput

# run fastqc to screen for overrepresented sequences and generates the analysis output in a text file. Alternatively SRR19*good_!(*singletons*) parenthized non-glob-match could be used
fastqc --threads 6 --outdir="$SRADIR/fastqcoutput" SRR19*good_????.fastq

# Unzip outfast to get .txt for parser to work
cd $SRADIR/fastqcoutput
unzip '*.zip'

echo "fastqc Completed" | mail -s "fastqc completed" (UID)@student.le.ac.uk

