import os
import pandas as pd

"""
Calculates nonhla depth.
"""

dir_metrics="/groups/wyattgrp/users/amunzur/hla_pipeline/results/metrics/depth_nonhla_panel_target"
path_panel_regions="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/hla_panel/target_regions_nonHLA.bed"
DIR_read_counts="/groups/wyattgrp/users/amunzur/hla_pipeline/results/metrics/raw_fastq_read_counts"
path_output_per_gene="/groups/wyattgrp/users/amunzur/hla_pipeline/results/metrics/compiled_depth/nonhla_depth_per_gene.tsv"
path_output_general="/groups/wyattgrp/users/amunzur/hla_pipeline/results/metrics/compiled_depth/nonhla_depth.tsv"

file_list = [os.path.join(dir_metrics, f) for f in os.listdir(dir_metrics)]

def annotate_genes(df, path_panel_regions):
    """
    Adds gene ids based on regions in the df.
    """
    panel_df =pd.read_csv(path_panel_regions, sep="\t", header = None, names=['chrom', 'start', 'end', "gene"])
    gene_min_pos = panel_df.groupby("gene")["start"].min().reset_index()
    gene_max_pos = panel_df.groupby("gene")["end"].max().reset_index()
    gene_chrom=panel_df[["chrom", "gene"]].drop_duplicates()
    
    gene_loc = gene_chrom.merge(gene_min_pos).merge(gene_max_pos)
    
    df['gene'] = None
    
    # Iterate over each region in gene_loc to annotate df
    for _, row in gene_loc.iterrows():
        chrom, gene, start, end = row['chrom'], row['gene'], row['start'], row['end']
        
        # Annotate rows in df where chrom and pos overlap with the current gene region
        df.loc[(df['chrom'] == chrom) & (df['pos'] >= start) & (df['pos'] <= end), 'gene'] = gene
    
    return df

    
df_list=[]
for f in file_list:
    sample_name = os.path.basename(f).replace(".txt", "")
    print(sample_name)
    df = pd.read_csv(f, sep = "\t", header = None, names=['chrom', 'pos', 'depth'])
    
    # Add genes to depth df
    df = annotate_genes(df, path_panel_regions)
    
    # Calculate median depth
    df = df.groupby("gene")["depth"].median().reset_index()
    
    # Normalize by read counts
    path_read_counts=os.path.join(DIR_read_counts, sample_name+".txt")
    df_read_counts=pd.read_csv(path_read_counts, header=None, sep="\t")
    nreads=df_read_counts[1].values[0]
    nreads_mils=nreads/10E6
    df["depth_norm"]=round(df["depth"]/nreads_mils)
    
    # Aes
    df["Sample_name"] = sample_name
    df = df[["Sample_name", "gene", "depth", "depth_norm"]]
    
    df_list.append(df)

concat_df = pd.concat(df_list).reset_index(drop = True)
concat_df.to_csv(path_output_per_gene, sep="\t", index=False)

###################################3
# Also do the general statistic, not per gene.

mydict={}
for f in file_list:
    sample_name = os.path.basename(f).replace(".txt", "")
    print(sample_name)
    df = pd.read_csv(f, sep = "\t", header = None, names=['chrom', 'pos', 'depth'])
    
    # Add genes to depth df
    median_val=df["depth"].median()
    
    mydict[sample_name]=median_val
        
    # Normalize by read counts
    path_read_counts=os.path.join(DIR_read_counts, sample_name+".txt")
    df_read_counts=pd.read_csv(path_read_counts, header=None, sep="\t")
    nreads=df_read_counts[1].values[0]
    nreads_mils=nreads/10E6
    depth_norm=int(median_val/nreads_mils)
    
    mydict[sample_name]={"depth": median_val, "depth_norm": depth_norm}

df2=pd.DataFrame(mydict).T.reset_index()
df2.to_csv(path_output_general, sep="\t", index=False)
