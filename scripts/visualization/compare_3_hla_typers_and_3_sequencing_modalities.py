import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os
import re
import numpy as np
from scipy.stats import spearmanr
import subprocess
import matplotlib.gridspec as gridspec

"""
Oncoprint style, compares the concordance of the calls made by Optitype, Polysolver and LILAC on targeted data.
"""
text_color_dict={True: "black", False: "red"}

path_hla_types_compiled="/groups/wyattgrp/users/amunzur/hla_pipeline/results/compiled_hla_types.csv"

hla_types=pd.read_csv(path_hla_types_compiled)
hla_types=hla_types[hla_types["Polysolver_wes"]!="NO WES"]
hla_types["Polysolver_targeted"] = hla_types["Polysolver_targeted"].str.replace(r"(:[^:]+):.*", r"\1", regex=True)
hla_types["Polysolver_wes"] = hla_types["Polysolver_wes"].str.replace(r"(:[^:]+):.*", r"\1", regex=True)

hla_types = hla_types.applymap(lambda x: re.sub(r".*\*", "", x) if isinstance(x, str) else x)
hla_types["LILAC_wgs"]=hla_types["LILAC_wgs"].fillna("No WGS")

def set_up_grid():
    global allele_num_ax, targ_poly_ax, targ_opt_ax, targ_lilac_ax, targ_conc_ax
    global wes_poly_ax, wes_opt_ax, wes_lilac_ax, wes_conc_ax, wgs_lilac_ax
    
    fig = plt.figure(figsize=(12, 6))
    outer_gs = gridspec.GridSpec(1, 6, width_ratios=[0.2, 1, 0.2, 1, 0.2, 0.4], hspace = 0, wspace = 0.2)
    targeted_gs = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer_gs[1], width_ratios=[1, 1, 1], hspace = 0, wspace = 0)
    wes_gs=gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer_gs[3], width_ratios=[1, 1, 1], hspace = 0, wspace = 0)
    
    allele_num_ax=plt.subplot(outer_gs[0])
    
    targ_poly_ax=plt.subplot(targeted_gs[0])
    targ_opt_ax=plt.subplot(targeted_gs[1])
    targ_lilac_ax=plt.subplot(targeted_gs[2])
    targ_conc_ax=plt.subplot(outer_gs[2])
    
    wes_poly_ax=plt.subplot(wes_gs[0])
    wes_opt_ax=plt.subplot(wes_gs[1])
    wes_lilac_ax=plt.subplot(wes_gs[2])
    wes_conc_ax=plt.subplot(outer_gs[4])
    
    wgs_lilac_ax=plt.subplot(outer_gs[5])
    
    # Set limits so the text fits within the figure
    for ax, title in zip(
        [allele_num_ax, targ_poly_ax, targ_opt_ax, targ_lilac_ax, targ_conc_ax, wes_poly_ax, wes_opt_ax, wes_lilac_ax, wes_conc_ax, wgs_lilac_ax],
        ["Allele", "Polysolver", f"Targeted\nOptitype", "Lilac", "Targeted\nconcordance?", "Polysolver", "WES\nOptitype", "Lilac", "WES\nconcordance?", "WGS\nLilac"]):
        ax.set_xlim(0, 1)
        ax.set_ylim(-1, len(df))  # Adjust based on the number of rows in df
        # ax.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
        ax.set_yticks(range(0, len(df)))
        ax.set_xlim((0.45, 0.55))
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_xticklabels([])
    
    targ_poly_ax.spines['right'].set_visible(False)
    targ_opt_ax.spines[['right', 'left']].set_visible(False)
    targ_lilac_ax.spines['left'].set_visible(False)
    
    wes_poly_ax.spines['right'].set_visible(False)
    wes_opt_ax.spines[['right', 'left']].set_visible(False)
    wes_lilac_ax.spines['left'].set_visible(False)
    
    for ax in [targ_poly_ax, targ_opt_ax, targ_lilac_ax, targ_conc_ax, wes_poly_ax, wes_opt_ax, wes_lilac_ax, wes_conc_ax, wgs_lilac_ax]:
        ax.tick_params(axis="y", labelleft=False, left=False)
    
    return(fig, outer_gs)

def plot_hla_types(df):
    for i, row in df.iterrows():
        if i % 2 == 0:
            facecolor="gray"
        else:
            facecolor="white"
        allele_num_ax.text(0.5, i, row["Allele"], fontsize=10, color='black', ha='center', va='center', bbox=dict(facecolor=facecolor, edgecolor='none', alpha=0.3, pad=0))
        
        targ_poly_ax.text(0.5, i, row["Polysolver_targeted"], fontsize=10, color='black', ha='center', va='center', bbox=dict(facecolor=facecolor, edgecolor='none', alpha=0.3, pad=0))
        targ_opt_ax.text(0.5, i, row["Optitype_targeted"], fontsize=10, color='black', ha='center', va='center', bbox=dict(facecolor=facecolor, edgecolor='none', alpha=0.3, pad=0))
        targ_lilac_ax.text(0.5, i, row["LILAC_targeted"], fontsize=10, color='black', ha='center', va='center', bbox=dict(facecolor=facecolor, edgecolor='none', alpha=0.3, pad=0))
        
        color_targ="black" if row["Concordance_targeted"] else "red"
        targ_conc_ax.text(0.5, i, row["Concordance_targeted"], fontsize=10, color=color_targ, ha='center', va='center', bbox=dict(facecolor=facecolor, edgecolor='none', alpha=0.3, pad=0))
        
        wes_poly_ax.text(0.5, i, row["Polysolver_wes"], fontsize=10, color='black', ha='center', va='center', bbox=dict(facecolor=facecolor, edgecolor='none', alpha=0.3, pad=0))
        wes_opt_ax.text(0.5, i, row["Optitype_wes"], fontsize=10, color='black', ha='center', va='center', bbox=dict(facecolor=facecolor, edgecolor='none', alpha=0.3, pad=0))
        wes_lilac_ax.text(0.5, i, row["LILAC_wes"], fontsize=10, color='black', ha='center', va='center', bbox=dict(facecolor=facecolor, edgecolor='none', alpha=0.3, pad=0))
        
        color_wes = "black" if row["Concordance_wes"] == "True" else "red"
        wes_conc_ax.text(0.5, i, row["Concordance_wes"], fontsize=10, color=color_wes, ha='center', va='center', bbox=dict(facecolor=facecolor, edgecolor='none', alpha=0.3, pad=0))
        wgs_lilac_ax.text(0.5, i, row["LILAC_wgs"], fontsize=10, color='black', ha='center', va='center', bbox=dict(facecolor=facecolor, edgecolor='none', alpha=0.3, pad=0))
    
    allele_num_ax.set_yticklabels(df["Patient"])

for gene in ["A", "B", "C"]:
    df=hla_types[hla_types['Allele'].str.contains(gene)].sort_values(by=["Patient", "Allele"], ascending=False).reset_index(drop=True)
    fig, outer_gs = set_up_grid()
    plot_hla_types(df)
    fig.suptitle(f"HLA {gene}/HLA typing concordance across 3 callers and 3 sequencing modalities")
    fig.savefig(f"/groups/wyattgrp/users/amunzur/hla_pipeline/results/figures/hla_types/{gene}_caller_and_sequencing_modality_concordance.png")
    fig.savefig(f"/groups/wyattgrp/users/amunzur/hla_pipeline/results/figures/hla_types/{gene}_caller_and_sequencing_modality_concordance.pdf")