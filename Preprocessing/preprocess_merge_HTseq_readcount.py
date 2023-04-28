#!/usr/bin/env python

"""
***********************************************************************************************************
Authors: Friday
Modified Date: 16-03-2023 02:00 BST
license:
Version = "1.0"

# Accepted argument format, e.g. python3 merge_readcount.py 'SRR197*' acc_list.txt > GSE67835_processed_matrix.tsv, where the 3rd argument contains the name to the files accepting wildcard and 4th argument is the name to the accession list, in which the accession number.(filenameExtension) are speparated in different rows. If a list of filename would be input, the variable sample counts and subsequent 2 lines would have to be adjusted accordingly, e.g. count_filename = sample_counts.split(",") (where the filenames are fed as a comma-seperated single argument to the function.)

The script is written in python 3.9.12. It was tested and meant to run smoothly under the same environment.
https://www.reformattext.com/sequential-number-generator.html
for f in *.csv; do mv "$f" "${f/_*.csv/.csv}" ; done
************************************************************************************************************
"""
import sys, glob

sample_counts = sys.argv[1]
count_files = glob.glob(sample_counts)   # to allow for quoted wildcard input, e.g. python3 merge_readcount.py 'SRR*'
count_filename = sorted(count_files)

acc_list_input = sys.argv[2]
with open(acc_list_input) as acl:
    acc_list = acl.read().splitlines() # alternatively [acl.rstrip("\n") for acl in open(acc_list_input)]
acl.close()

# define a dictionary
count_dict = {}

for count_file in count_filename:
    with open(count_file) as count:
        for line in count:
            if line.startswith("_"):   # skipping the the descritor at the tail
                continue
            line = line.rstrip("\n")
            element = line.split("\t")
            #print(element)
            if element[0] in count_dict:
                count_dict[element[0]].append(element[1])
            else:
                count_dict[element[0]] = [element[1]]

acc_list = [accession.split('.')[0] for accession in acc_list]

print("gene_id\t"+"\t".join(acc_list))
for gene_id in count_dict:
    print(gene_id + "\t" + "\t".join(count_dict[gene_id]))
