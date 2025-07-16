import os
import pandas as pd

"""
Calculates median depth in hla alleles - allele specific doesn't lump alelle1 and 2 into one single metric.
"""

DIR_allele_depth="/groups/wyattgrp/users/amunzur/hla_pipeline/results/metrics/depth_hla_allele_specific"
DIR_read_counts="/groups/wyattgrp/users/amunzur/hla_pipeline/results/metrics/raw_fastq_read_counts"
path_output="/groups/wyattgrp/users/amunzur/hla_pipeline/results/metrics/compiled_depth/hla_allele_specific_depth.tsv"

all_files=[os.path.join(DIR_allele_depth, f) for f in os.listdir(DIR_allele_depth)]

mylist=[]
for f in all_files:
    sample_name=os.path.basename(f)
    df = pd.read_csv(f, sep="\t", names=["allele", "pos", "depth"])
    df = df.groupby("allele")["depth"].median().reset_index()
    
    # Normalize by number of reads
    path_read_counts=os.path.join(DIR_read_counts, sample_name)
    df_read_counts=pd.read_csv(path_read_counts, header=None, sep="\t")
    nreads=df_read_counts[1].values[0]
    nreads_mils=nreads/10E6
    df["depth_norm"]=round(df["depth"]/nreads_mils)
    
    # Reorder cols
    df["Sample_name"]=sample_name.replace(".txt", "")
    df=df[["Sample_name", "allele", "depth", "depth_norm"]]
    mylist.append(df)
    
df_major=pd.concat(mylist).reset_index(drop=True)
df_major.to_csv(path_output, sep="\t", index=False)