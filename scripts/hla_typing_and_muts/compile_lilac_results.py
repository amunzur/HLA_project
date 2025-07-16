import os
import pandas as pd

"""
Compiles LILAC results from multiple samples into one csv file.
"""

def list_all_files(directory):
    all_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            all_files.append(os.path.join(root, file))
    return all_files

DIR_lilac="/groups/wyattgrp/users/amunzur/hla_pipeline/results/lilac/"

all_files=list_all_files(DIR_lilac)
df_list=[]

for f in all_files:
    if "WGS" in f and f.endswith(".lilac.tsv"):
        sample_name=os.path.basename(f).replace(".lilac.tsv", "")
        df=pd.read_csv(f, sep="\t")
        df=df[['Allele', 'RefTotal', 'RefUnique', 'RefWild']]
        df["Sample"]=sample_name
        df=df[['Sample', 'Allele', 'RefTotal', 'RefUnique']]
        df_list.append(df)

combined_df=pd.concat(df_list).reset_index(drop=True)
combined_df.to_csv("/groups/wyattgrp/users/amunzur/hla_pipeline/results/lilac/ALL_WGS_hla_types_from_LILAC.csv", index=False)
