import pandas as pd
import numpy as np
import csv

def get_ens_dict(file_path):
    with open(file_path) as f:
        gtf = list(f)

    gtf = [x for x in gtf if not x.startswith("#")]
    gtf = [x for x in gtf if 'gene_id "' in x and 'gene_name "' in x]
    if len(gtf) == 0:
        print('you need to change gene_id " and gene_name " formats')

    gtf = list(map(lambda x: (x.split('gene_id "')[1].split('"')[0], x.split('gene_name "')[1].split('"')[0]), gtf))
    print("\nThere are " + str(len(gtf)) + " numbers of gene id/ transcript accession found in the reference genome.\n")

    # There are multiple entries and isoforms in gtf therefore wrapped as tuples, but not list of list.
    gtf = dict(set(gtf))   # to give unique element of the list

    print("There are " + str(len(gtf)) + " gene id entries in the gtf file.\n")
    unique_gene_counts = len(list(set(gtf.values())))   # Convert the value of a dict to a set to remove the duplicate and pass it back to a list
    print("There are " + str(unique_gene_counts) + " unique gene names in the gtf file, which includes different biotypes, such as protein coding, long non-coding RNA(lincRNA), ncRNA, pseudogene,etc..\n")
    return gtf

gtf_dict = get_ens_dict("Homo_sapiens.GRCh37.75.gtf")
csv_file = "HSGRCh37_75_gtf_dict.csv"
try:
    with open(csv_file, "w") as csvfile:
        for key in gtf_dict.keys():
            csvfile.write("%s\t%s\n" % (key, gtf_dict[key]))
except IOError:
    print("I/O error")

# We manipute the original CPM transformed count matrix as the base10 logarithm count matrix is harder to manipute and
# we cant do simple summation because log of the same base can be added together by multiplying their arguments: log(xy) = log(x) + log(y)
df = pd.read_csv("GSE67835_CPM.tsv", delimiter="\t")

# Mapping of gene names to their id, then dropped the id column and set sorted name as index.
df["Gene_Name"] = df["gene_id"].map(gtf_dict)
df = df.drop("gene_id", axis=1)
df.sort_values(by="Gene_Name", inplace=True) # sort gene name by numeric-alphabetic order
print("There are " + str(int(df["Gene_Name"].nunique())) + " genes mapped to the GRCh37.75 reference genome.\n")
df = df.groupby(df["Gene_Name"]).aggregate('sum')  # groupby function automatically set_index("Gene_Name", inplace=True)

# To drop rows with All zeros, alternatively use df.loc[(df!=0).any(1)]
df1 = df.loc[~(df==0).all(axis=1)]
dropped_row_count = len(df.index) - len(df1.index)
print("There are " + str(dropped_row_count) + " gene names not mapped to any reads, and, hence, dropped.")
#df1.to_csv("GSE67835_CPM_named_clean.tsv", index=True, sep="\t")

with np.errstate(divide='ignore'):
    df1 = np.log10(df1).replace([np.inf, -np.inf], 0)
df1.to_csv("GSE67835_transformed_log_notrounded.tsv", index=True, sep="\t")

with np.errstate(divide='ignore'):   # to ignore the runtime warnings with logging zeros
    df1 = np.round(np.log10(df1).replace([np.inf, -np.inf], 0)).astype(int)
df1.to_csv("GSE67835_transformed_rounded.tsv", index=True, sep="\t")
