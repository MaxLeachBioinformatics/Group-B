#!/usr/bin/env python

"""
***********************************************************************************************************
Authors: Friday
Modified Date: 23-03-2023 09:00 BST
license:
Version = "1.0"

This script reads Log.final.out merged by the STAR awk script convert the total number of reads into a
pandas series then use the total number of read for individual cell to normalise the count matrix generated
by HTseq. Then the CPM value is converted to a log base 10 with infinite value replaced by 0.

The script is written in python 3.9.12. It was tested and meant to run smoothly under the same environment.
************************************************************************************************************
"""

import pandas as pd
import numpy as np


# Read the STAR Log.final.out merged file and read into pandas
df = pd.read_csv("GSE67835_processed_Log.final.out", delimiter=";", index_col=0)
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
df.columns = df.columns.str.replace("./", "", regex=False)
df = df.drop("Unnamed: 467", axis=1)
df.columns = df.columns.str.replace("_goodLog", "")
df = df.rename(columns=lambda x: x.split('.')[0]).loc["Number of input reads"]   # Number of input reads at iloc[4]
df = pd.to_numeric(df)

# Read the count matrix into pandas
df1 = pd.read_csv("GSE67835_processed_matrix.tsv", delimiter="\t", index_col=0)
df1.astype(int)

#The counts of all genes for any given cell where converted to counts per million
#(CPM) by diving with the total number of reads and multiplying by 10^6 followed by
# conversion to a log base 10.
out = np.log10(df1*1000000/df)
out.replace([np.inf, -np.inf], 0, inplace=True)
print(out.head(10))

#genreate an output
out.to_csv("GSE67835_CPM_log10.tsv", sep='\t')