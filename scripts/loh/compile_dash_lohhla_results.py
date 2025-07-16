#!/bin/bash
import os 
import pandas as pd
import numpy as np
import glob

# This script combines outputs from dash from many samples and makes one "LOH" output file.
DIR_dash="/groups/wyattgrp/users/amunzur/hla_pipeline/results/dash"
DIR_lohhla="/groups/wyattgrp/users/amunzur/hla_pipeline/results/lohhla"
DIR_out="/groups/wyattgrp/users/amunzur/hla_pipeline/results/dash_lohhla_compiled"
path_compiled_hla_types="/groups/wyattgrp/users/amunzur/hla_pipeline/results/compiled_hla_types.csv"
path_sample_info="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test.tsv"

# SETP 1. COMPILE DASH RESULTS
# Compile DASH
summary_LOH=[] 
detailed_LOH=[]
sample_info_list=[]
for subdir in os.listdir(DIR_dash):
    sub_path = os.path.join(DIR_dash, subdir)
    if os.path.isdir(sub_path):
        if os.path.isfile(os.path.join(sub_path, "DASH.output.txt")): # check if dash output exists
            dash_raw=pd.read_csv(os.path.join(sub_path, "DASH.output.txt"), "\t")[["Alleles", "DASH_deletion", "Purities", "Ploidies"]]
            dash=dash_raw.copy()[["Alleles", "DASH_deletion"]]
            dash.set_index('Alleles', inplace=True)
            dash=dash.T
            gene_alleles = {
                'HLA-A': list(dash.columns[0:2]),
                'HLA-B': list(dash.columns[2:4]),
                'HLA-C': list(dash.columns[4:6])
                }
            # SUMMARY LOH
            df = dash.copy()
            summary_results = []
            for gene, alleles in gene_alleles.items():
                # Get the DASH_deletion values for the alleles
                allele_values = [df.loc['DASH_deletion', allele] for allele in alleles]
                if not alleles[0] == alleles[1]: # Check if the alleles are different
                    if all(allele_values): # if both are true (deep deletion)
                        result = 2
                    elif any(allele_values): # if one of them is true (LOH)
                        result = 1
                    else: # both are false, which means no deletion
                        result = 0
                else:
                    # The alleles are homozygous, so we can't have LOH.
                    result = "Homozygous"
                # Add the result to the list of results
                summary_results.append(result)
            # Create a DataFrame from the results
            result_df = pd.DataFrame({'Sample': subdir, 'Gene': list(gene_alleles.keys()), 'Deletion': summary_results})
            result_df = result_df[["Sample", "Gene", "Deletion"]]
            result_df.sort_values(by=['Gene'], inplace=True)
            result_df = result_df.pivot(index='Sample', columns='Gene', values='Deletion')
            # Sample information
            sample_info=pd.DataFrame({'Sample': subdir, "Ploidy":dash_raw["Ploidies"].unique(), "Purity":dash_raw["Purities"].unique()}).set_index("Sample")
            result_df=pd.merge(result_df, sample_info, left_index=True, right_index=True)
            # add back to the results df
            summary_LOH.append(result_df)
            # LOH IN DETAIL
            df = pd.DataFrame(np.vstack([df.columns, df]), columns = ["HLA_A_1", "HLA_A_2", "HLA_B_1", "HLA_B_2", "HLA_C_1", "HLA_C_2"]).T
            df.columns = ["Alleles", "Deletion"]
            df["Sample"] = subdir
            df = df[["Sample", "Alleles", "Deletion"]]
            detailed_LOH.append(df)
        else: 
            print(sub_path + " has no DASH output.")

summary_LOH = pd.concat(summary_LOH)
detailed_LOH = pd.concat(detailed_LOH)

summary_LOH.to_csv(os.path.join(DIR_out, "summary_LOH.csv"), index = True)  
detailed_LOH.to_csv(os.path.join(DIR_out, "detailed_LOH.csv"), index = True)

detailed_LOH["Patient"]=detailed_LOH["Sample"].str.replace(r"(_WBC|_cfDNA|_FiT).*", "", regex=True)
detailed_LOH=detailed_LOH.rename(columns={"Deletion": "Deletion DASH"})

# STEP 2. 
# Compile LOHHLA
lohhla_results=[]

files = glob.glob(os.path.join(DIR_lohhla, "**", "*DNA.HLAlossPrediction_CI.xls"), recursive=True)
for file in files:
    try:
        # Annotate allele loss in DASH style
        df = pd.read_csv(file, sep="\t")
        df=df[["region", "LossAllele", "KeptAllele"]].drop_duplicates()
        df_melted=df.melt(id_vars=["region"], value_vars=["LossAllele", "KeptAllele"])
        df_melted["Deletion"]=df_melted["variable"]=="LossAllele"
        df_melted.drop("variable", axis=1, inplace=True)
        df_melted.columns=["Sample", "Alleles", "Deletion"]
        
        lohhla_results.append(df_melted)
    except pd.errors.EmptyDataError:
        print(f"Skipping empty file: {file}")

# Merge all lohhla calls
lohhla_merged=pd.concat(lohhla_results, ignore_index=True).sort_values(by=["Sample", "Alleles"]).rename(columns={"Deletion": "Deletion Lohhla"})
lohhla_merged["Patient"]=lohhla_merged["Sample"].str.replace(r"(_WBC|_cfDNA|_FiT).*", "", regex=True)
lohhla_merged.to_csv("/groups/wyattgrp/users/amunzur/hla_pipeline/results/loh/lohhla_results.csv", index=False)
# Merge lohhla and dash calls
# detailed_LOH.reset_index(drop=True).merge(lohhla_merged, how="outer")

# STEP 3.
# COMBINE DASH, LOHHLA AND HLA TYPES.
hla_types=pd.read_csv(path_compiled_hla_types)[["Patient", "Polysolver_targeted"]]
hla_types["Polysolver_targeted"] = "hla_"+ (
    hla_types["Polysolver_targeted"]
    .str.replace("*", "_", regex=False)
    .str.replace(":", "_")
    .apply(lambda x: x[0].lower() + x[1:] if x else x)  # Change only first letter to lowercase
)
hla_types.columns=["Patient", "Alleles"]

sample_info=pd.read_csv(path_sample_info, sep="\t").rename(columns={"Tumor_name": "Sample"}).rename(columns={"Patient_ID": "Patient"})[["Patient", "Sample"]]
hla_types=hla_types.merge(sample_info)
hla_types["Gene"]=hla_types["Alleles"].str.extract(r"(hla_[a-z]+)")

grouped = hla_types.groupby(["Patient", "Sample", "Gene"])
hla_types["Homozygous"] = grouped["Alleles"].transform(lambda x: x.nunique() == 1)
hla_types=hla_types[["Patient", "Sample", "Gene", "Alleles", "Homozygous"]]

hla_types["Allele_Number"] = (hla_types.groupby(["Patient", "Sample", "Gene"]).cumcount() + 1)
hla_types["Allele_Number"] = (hla_types["Gene"] + "_" + hla_types["Allele_Number"].astype(str))



loh_df=hla_types.merge(lohhla_merged, how="left", suffixes=["", "_lohhla"]).merge(detailed_LOH, how="left").drop_duplicates()


loh_df["Sample included in targeted test run?"] = loh_df["Patient"].isin(sample_info["Patient"])

# Fill in the Sample column
# del loh_df["Sample"]
# loh_df=loh_df.merge(sample_info, left_on="Patient", right_on="Patient_ID").drop("Patient_ID", axis=1)
# loh_df=loh_df[["Patient", "Sample", "Alleles", "Deletion Lohhla", "Deletion DASH", "Sample included in targeted test run?"]]

loh_df.to_csv("/groups/wyattgrp/users/amunzur/hla_pipeline/results/loh/lohhla_dash_compiled.csv", index=False)

#####################
# PLOTTING
#####################
DIR_dash_output="/groups/wyattgrp/users/amunzur/hla_pipeline/results/dash_compiled"
DIR_figures="/groups/wyattgrp/users/amunzur/hla_pipeline/results/figures/dash"

summ=summary_LOH.copy()
det=detailed_LOH.copy()

def replace_words(df):
    # Beautify the data frame for plotting
    for col_name in ["HLA-A", "HLA-B", "HLA-C"]:
        df[col_name].replace(0, "Intact", inplace=True)
        df[col_name].replace(1, "LOH", inplace=True)
        df[col_name].replace(2, "Deep deletion", inplace=True)
    return(df)

def return_plotting_df(df, gene_to_plot):
    values=pd.DataFrame(df[gene_to_plot].value_counts())
    color_df=pd.DataFrame({"Intact":"limegreen", "Homozygous":"grey", "LOH":"cornflowerblue", "Deep deletion":"blue"}, index=[0]).T
    color_df.columns=["Color"]
    df_plotting=pd.merge(values, color_df, left_index=True, right_index=True)
    df_explode=pd.DataFrame({"Intact":0, "Homozygous":0, "LOH":0.1, "Deep deletion":0}, index=[0]).T
    df_explode.columns=["Explode_values"]
    df_plotting=pd.merge(df_plotting, df_explode, left_index=True, right_index=True)
    return(df_plotting)

def plot_gene(df_plotting, plot_title, figure_name, DIR_figures):
    fig, ax = plt.subplots()
    ax.pie(df_plotting.iloc[:,0], autopct= lambda x: '{:.0f}'.format(x*df_plotting.iloc[:,0].sum()/100), labels=df_plotting.index, colors=df_plotting.iloc[:,1], explode = df_plotting.iloc[:,2])
    ax.set_title(plot_title)
    fig.savefig(os.path.join(DIR_figures, figure_name), dpi=199)

summ=replace_words(summ)

for gene in ["HLA-A", "HLA-B", "HLA-C"]:
    # All purities
    figure_ending=".png"
    df_plotting=return_plotting_df(summ, gene)
    plot_gene(df_plotting, gene, gene+figure_ending, DIR_figures)
    # Purity >= 20
    figure_ending="_purity0.20.png"
    summ_pur=summ[summ['Purity'] > 0.20]
    df_plotting=return_plotting_df(summ_pur, gene)
    plot_gene(df_plotting, gene, gene+figure_ending, DIR_figures)
    




loh_df[(loh_df["Homozygous"]==False) & (loh_df["Deletion DASH"]==True)]