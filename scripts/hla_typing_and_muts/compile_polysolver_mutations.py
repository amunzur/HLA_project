import os
import pandas as pd

dir_muts="/groups/wyattgrp/users/amunzur/hla_pipeline/results/polysolver"
path_sample_list="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test.tsv"

samples=pd.read_csv(path_sample_list, sep="\t")

df_mutect=[]
df_strelka=[]
for tumor_sample in os.listdir(dir_muts):
    print(tumor_sample)
    path_mutations=os.path.join(dir_muts, tumor_sample, "hla_mutations", f"{tumor_sample}.mutect.filtered.nonsyn.annotated")
    path_strelka=os.path.join(dir_muts, tumor_sample, "hla_mutations", f"{tumor_sample}.strelka_indels.filtered.annotated")
    print(f"Loading {path_mutations}")
    muts_df=pd.read_csv(path_mutations, sep="\t", index_col=False)
    indels_df=pd.read_csv(path_strelka, sep="\t")
    
    df_mutect.append(muts_df)
    df_strelka.append(indels_df)


concatted_muts=pd.concat(df_mutect)
concatted_indels=pd.concat(df_strelka)

concatted_muts.to_csv("/groups/wyattgrp/users/amunzur/hla_pipeline/results/variant_calling/hla_mutations_polysolver.csv", index=False)