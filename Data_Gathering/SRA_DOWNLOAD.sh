#!/bin/bash

#Load Modules
module load sratoolkit/2.11.1

# Navigate to the directory containing the script
cd "$(dirname "$0")"



# Read accession numbers from the file and store in an array
mapfile -t ACCESSIONS < Ascension_list.txt

cd $SCRATCHDIR

# Create a new directory for storing downloaded files
mkdir -p sradata

cd sradata
mkdir -p fastqdump
# Loop through each accession number and download the corresponding file
for ACCESSION in "${ACCESSIONS[@]}"
do
  echo "Downloading file for accession: $ACCESSION"
  prefetch $ACCESSION
  fasterq-dump $ACCESSION -O ./fastqdump/
done
echo "SRA download complete" | mail -s "SRA download complete" your_email@example.com
