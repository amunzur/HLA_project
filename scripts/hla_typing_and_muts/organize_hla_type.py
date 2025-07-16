"""
Organizes HLA types from Polysolver, Optitype and LILAC into one file.
"""

import os
import pandas as pd
import numpy as np

path_sample_list_targeted="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test.tsv"
path_sample_list_wes="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/WES_samples.tsv"
path_sample_list_wes="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/WES_samples.tsv"

path_select_samples_lilac_wgs="/groups/wyattgrp/users/amunzur/hla_pipeline/results/lilac/WGS_hla_types_from_LILAC.csv"

lilac_dir="/groups/wyattgrp/users/amunzur/hla_pipeline/results/lilac"
optitype_dir="/groups/wyattgrp/users/amunzur/hla_pipeline/results/optitype"
polysolver_dir="/groups/wyattgrp/users/amunzur/hla_pipeline/results/polysolver/hla_types"

def list_all_files(directory):
    """
    Lists all files in dir.
    """
    all_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            all_files.append(os.path.join(root, file))
    return all_files

def return_lilac_results(sample_name, lilac_dir):
    """
    Given a single sample name, return the result from LILAC.
    """
    sample_dir=os.path.join(lilac_dir, sample_name)
    all_files=list_all_files(sample_dir)
    lilac_calls=[f for f in all_files if f.endswith("lilac.tsv")]
    if len(lilac_calls)==1:
        df=pd.read_csv(lilac_calls[0], sep="\t")
        if not df.empty:
            calls=df["Allele"]
            calls.index=['A', 'A', 'B', 'B', 'C', 'C']
            calls = calls.groupby(level=0, group_keys=False).apply(lambda x: x.sort_values())
            calls.index=['A1', 'A2', 'B1', 'B2', 'C1', 'C2']
        else:
            calls=pd.Series(np.nan, index=['A1', 'A2', 'B1', 'B2', 'C1', 'C2'])
    else:
        calls=pd.Series(np.nan, index=['A1', 'A2', 'B1', 'B2', 'C1', 'C2'])
    return(calls)

def return_optitype_results(sample_name, optitype_dir):
    """
    Given a single sample name, return the result from Optitype.
    """
    sample_dir=os.path.join(optitype_dir, sample_name)
    all_files=list_all_files(sample_dir)
    optitype_calls=[f for f in all_files if f.endswith("result.tsv")]
    if len(optitype_calls)==1:
        calls=pd.read_csv(optitype_calls[0], sep="\t").iloc[0, :][['A1', 'A2', 'B1', 'B2', 'C1', 'C2']]
        calls.index=['A', 'A', 'B', 'B', 'C', 'C']
        calls = calls.groupby(level=0, group_keys=False).apply(lambda x: x.sort_values())
        calls.index=['A1', 'A2', 'B1', 'B2', 'C1', 'C2']
    else:
        calls=pd.Series(np.nan, index=['A1', 'A2', 'B1', 'B2', 'C1', 'C2'])
    return(calls)

def return_polysolver_results(sample_name, polysolver_dir):
    """
    Given a single sample name, return the result from Polysolver.
    """
    path_hla_calls=os.path.join(polysolver_dir, sample_name, "winners.hla.txt")
    if os.path.exists(path_hla_calls):
        df=pd.read_csv(path_hla_calls, sep="\t", header=None, names=["Gene", "Allele1", "Allele2"])
        melted_df = pd.melt(df, id_vars="Gene", value_vars=["Allele1", "Allele2"],var_name="AlleleType",value_name="Call")
        calls=melted_df["Call"].to_frame()
        calls.index=['A', 'B', 'C', 'A', 'B', 'C']
        calls=calls.groupby(level=0, group_keys=False).apply(lambda x: x.sort_values(by="Call"))
        
        # Extract allele name in A*XX:YY format
        calls_extracted = calls["Call"].str.extract(r'hla_([abc])_([\d_]+)')
        def format_hla(row):
            hla = f"{row.iloc[0].upper()}*"
            numeric_parts = row.iloc[1:].dropna().astype(str).tolist()
            hla += ":".join(part.replace("_", ":") for part in numeric_parts)
            return hla
        
        calls_formatted = calls_extracted.apply(format_hla, axis=1).replace("_", ":")
        calls_formatted.name = "Allele"
        calls_formatted.index=['A1', 'A2', 'B1', 'B2', 'C1', 'C2']
    else:
        calls_formatted=pd.Series(np.nan, index=['A1', 'A2', 'B1', 'B2', 'C1', 'C2'])
    
    return(calls_formatted)

def check_hla_concordance(df):
    """
    Checks if HLA calls are in perfect concordance across Polysolver, LILAC, and Optitype, considering allele swaps.
    """
    df_mod=df.copy()
    if not all(df_mod["Polysolver"].isna()):
        df_mod["Polysolver"]=df_mod["Polysolver"].apply(lambda x: ":".join(x.split(":")[:2]))
    df["Concordance"]=df_mod.nunique(axis=1) == 1
    return(df)

def compile_caller_results(path_sample_list, sequencing_modality):
    """
    Compiles results from all three callers.
    """
    results=[]
    wbc_samples=pd.read_csv(path_sample_list, sep="\t")["WBC_name"]
    for wbc in wbc_samples:
        lilac=return_lilac_results(wbc, lilac_dir).to_frame(name="Allele").assign(Caller="LILAC").reset_index()
        optitype=return_optitype_results(wbc, optitype_dir).to_frame(name="Allele").assign(Caller="Optitype").reset_index()
        polysolver=return_polysolver_results(wbc, polysolver_dir).to_frame(name="Allele").assign(Caller="Polysolver").reset_index()
        merged=polysolver.merge(lilac, on="index").merge(optitype, on="index").set_index("index")
        merged.columns = ["Polysolver", "Polysolver_Caller", "LILAC", "LILAC_Caller", "Optitype", "Optitype_Caller"]
        merged = merged[["Polysolver", "LILAC", "Optitype"]]  # Keeping only allele calls
        merged=check_hla_concordance(merged)
        merged["Sample"]=wbc
    
        results.append(merged)
    
    df=pd.concat(results).reset_index(drop=True)
    df["Sequencing"]=sequencing_modality
    
    # Generate allele column
    df["Base"] = df["Polysolver"].combine_first(df["LILAC"]).combine_first(df["Optitype"]).apply(lambda x: x.split("*")[0])  
    df["Number"] = df.groupby("Base").cumcount() % 2 + 1
    df["Allele"] = df["Base"] + df["Number"].astype(str)
    df = df.drop(columns=["Base", "Number"])
    return(df)

# Targeted seq
targeted=compile_caller_results(path_sample_list_targeted, sequencing_modality="Targeted").drop_duplicates()
wes=compile_caller_results(path_sample_list_wes, sequencing_modality="WES").drop_duplicates()
wes["Sample"]=wes["Sample"].str.replace("_WES", "")

targeted["Patient"] = targeted["Sample"].str.replace(r"_WBC.*", "", regex=True)
wes["Patient"] = wes["Sample"].str.replace(r"_WBC.*", "", regex=True)

targeted.columns=['Polysolver_targeted','LILAC_targeted','Optitype_targeted','Concordance_targeted','Sample','Sequencing',"Allele","Patient"]
targeted=targeted[['Patient','Allele','Polysolver_targeted','LILAC_targeted','Optitype_targeted','Concordance_targeted']]

wes.columns=['Polysolver_wes','LILAC_wes','Optitype_wes','Concordance_wes','Sample','Sequencing','Allele','Patient']
wes=wes[['Patient','Allele','Polysolver_wes','LILAC_wes','Optitype_wes','Concordance_wes']]

# Compare targeted to wes
merged=targeted.merge(wes, on=["Patient", "Allele"], how="outer")
merged[["Polysolver_wes", "LILAC_wes", "Optitype_wes", "Concordance_wes", ]]=merged[["Polysolver_wes", "LILAC_wes", "Optitype_wes", "Concordance_wes", ]].fillna("NO WES")

# Add data from WGS
# 1. LILAC
lilac_df_wgs=pd.read_csv(path_select_samples_lilac_wgs)
lilac_df_wgs["Patient"]=lilac_df_wgs["Sample"].str.replace(r"_WBC.*", "", regex=True)
lilac_df_wgs=lilac_df_wgs[['Patient', 'Allele']].rename(columns={"Allele":"LILAC_wgs"})
lilac_df_wgs["Base"] = lilac_df_wgs["LILAC_wgs"].apply(lambda x: x.split("*")[0])  
lilac_df_wgs.drop_duplicates(["Patient", "Base"])
lilac_df_wgs = lilac_df_wgs.sort_values(by=['Patient', 'Base', 'LILAC_wgs'])
lilac_df_wgs["Number"] = lilac_df_wgs.groupby("Base").cumcount() % 2 + 1
lilac_df_wgs["Allele"] = lilac_df_wgs["Base"] + lilac_df_wgs["Number"].astype(str)
lilac_df_wgs = lilac_df_wgs.drop(columns=["Base", "Number"])

merged=merged.merge(lilac_df_wgs, how="left")
merged.to_csv("/groups/wyattgrp/users/amunzur/hla_pipeline/results/compiled_hla_types.csv", index=False)