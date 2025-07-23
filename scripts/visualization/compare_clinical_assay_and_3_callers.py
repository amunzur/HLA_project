import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os
import re
import numpy as np
import subprocess
import matplotlib.gridspec as gridspec

"""
Oncoprint style, compares the concordance of the calls made by Optitype, Polysolver and LILAC on targeted data.
"""
text_color_dict={True: "black", False: "red"}

# hla_types = pd.read_csv('https://docs.google.com/spreadsheets/d/19Gm3X6io2Fgcrui8qILxeSzpwNW517siLnhyKETCquM/export?format=csv&gid=734748589')
hla_types=pd.read_csv("/groups/wyattgrp/users/amunzur/hla_project/data/HLA_typing/HLA - WES and WGS list - Clinical HLA typing .csv")
hla_types=hla_types[['Patient', 'Allele', "Clinical assay", 'Polysolver_targeted', 'LILAC_targeted', 'Optitype_targeted', 'Polysolver_wes', 'LILAC_wes', 'Optitype_wes', "LILAC_wgs"]]

# Just first two digits only
hla_types["Polysolver_targeted"] = hla_types["Polysolver_targeted"].str.replace(r"(:[^:]+):.*", r"\1", regex=True)
hla_types["Polysolver_wes"] = hla_types["Polysolver_wes"].str.replace(r"(:[^:]+):.*", r"\1", regex=True)

hla_types = hla_types.applymap(lambda x: re.sub(r".*\*", "", x) if isinstance(x, str) else x)
hla_types["LILAC_wgs"]=hla_types["LILAC_wgs"].fillna("NO WGS")

# Determine match with clinical assay for each sequencing modality
hla_types["clinical assay two dig"]=hla_types["Clinical assay"].str.replace(r'^((?:[^:]*):[^:]*):.*', r"\1", regex=True)
for hla_caller in ["Polysolver", "LILAC", "Optitype"]:
    for seq_type in ["targeted", "wes"]:
        hla_types[f"{hla_caller}_{seq_type}_match"]=hla_types[f"{hla_caller}_{seq_type}"]==hla_types["clinical assay two dig"]
    
    hla_types.loc[hla_types[f"{hla_caller}_wes"]=="NO WES", f"{hla_caller}_wes_match"] = ""

hla_types["LILAC_wgs_match"]=hla_types["LILAC_wgs"]==hla_types["clinical assay two dig"]
hla_types.loc[hla_types["LILAC_wgs"]=="NO WGS", "LILAC_wgs_match"] = ""
hla_types.to_csv("/groups/wyattgrp/users/amunzur/hla_project/data/HLA_typing/clinical_grade_assay_vs_three_callers.csv", index=False)


def set_up_grid(df):
    global allele_num_ax, clinical_assay_ax, targ_poly_ax, targ_opt_ax, targ_lilac_ax, targ_conc_ax, targ_conc_poly_ax, targ_conc_opt_ax, targ_conc_lilac_ax
    global wes_poly_ax, wes_opt_ax, wes_lilac_ax, wes_conc_poly_ax, wes_conc_opt_ax, wes_conc_lilac_ax
    global wgs_lilac_ax, wgs_conc_lilac_ax
    
    fig = plt.figure(figsize=(20, 6))
    outer_gs = gridspec.GridSpec(1, 8, width_ratios=[0.2, 0.5, 1, 1, 1, 1, 0.4, 0.4], hspace = 0, wspace = 0.2)
    
    targeted_gs = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer_gs[2], width_ratios=[1, 1, 1], hspace = 0, wspace = 0)
    targ_poly_ax=plt.subplot(targeted_gs[0])
    targ_opt_ax=plt.subplot(targeted_gs[1])
    targ_lilac_ax=plt.subplot(targeted_gs[2])
    
    targeted_conc_gs=gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer_gs[3], width_ratios=[1, 1, 1], hspace = 0, wspace = 0)
    targ_conc_poly_ax=plt.subplot(targeted_conc_gs[0])
    targ_conc_opt_ax=plt.subplot(targeted_conc_gs[1])
    targ_conc_lilac_ax=plt.subplot(targeted_conc_gs[2])   
    
    wes_gs=gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer_gs[4], width_ratios=[1, 1, 1], hspace = 0, wspace = 0)
    wes_poly_ax=plt.subplot(wes_gs[0])
    wes_opt_ax=plt.subplot(wes_gs[1])
    wes_lilac_ax=plt.subplot(wes_gs[2])
    
    wes_conc_gs=gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer_gs[5], width_ratios=[1, 1, 1], hspace = 0, wspace = 0)
    wes_conc_poly_ax=plt.subplot(wes_conc_gs[0])
    wes_conc_opt_ax=plt.subplot(wes_conc_gs[1])
    wes_conc_lilac_ax=plt.subplot(wes_conc_gs[2])   
    
    allele_num_ax=plt.subplot(outer_gs[0])
    clinical_assay_ax=plt.subplot(outer_gs[1])
    
    wgs_lilac_ax=plt.subplot(outer_gs[6])
    wgs_conc_lilac_ax=plt.subplot(outer_gs[7])
    
    # Set limits so the text fits within the figure
    for ax, title in zip(
        [allele_num_ax, 
        clinical_assay_ax, 
        
        targ_poly_ax, 
        targ_opt_ax, 
        targ_lilac_ax, 
        
        targ_conc_poly_ax,
        targ_conc_opt_ax,
        targ_conc_lilac_ax,
        
        wes_poly_ax, 
        wes_opt_ax, 
        wes_lilac_ax,
        
        wes_conc_poly_ax,
        wes_conc_opt_ax,
        wes_conc_lilac_ax,
        
        wgs_lilac_ax, 
        
        wgs_conc_lilac_ax],
        
        ["Allele", 
        "Clinical\nassay", 
        
        "Polysolver", 
        f"Targeted\nOptitype", 
        "Lilac", 
        
        "Polysolver",
        "tNGS match\nOptitype",
        "Lilac",
        
        "Polysolver",
        "WES\nOptitype", 
        "Lilac",
        
        "Polysolver",
        "WES match\nOptitype",
        "Lilac",

        "WGS\nLilac",
        
        "WGS match\nLilac"]):
        
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
    
    targ_conc_poly_ax.spines['right'].set_visible(False)
    targ_conc_opt_ax.spines[['right', 'left']].set_visible(False)
    targ_conc_lilac_ax.spines['left'].set_visible(False)
    
    wes_conc_poly_ax.spines['right'].set_visible(False)
    wes_conc_opt_ax.spines[['right', 'left']].set_visible(False)
    wes_conc_lilac_ax.spines['left'].set_visible(False)
    
    for ax in [clinical_assay_ax, targ_poly_ax, targ_opt_ax, targ_lilac_ax, targ_conc_poly_ax, targ_conc_opt_ax, targ_conc_lilac_ax, wes_poly_ax, wes_opt_ax, wes_lilac_ax, wes_conc_poly_ax, wes_conc_opt_ax, wes_conc_lilac_ax, wgs_lilac_ax, wgs_conc_lilac_ax]:
        ax.tick_params(axis="y", labelleft=False, left=False)
    
    return(fig, outer_gs)


def set_up_grid(df):
    global allele_num_ax, clinical_assay_ax, targ_poly_ax, targ_opt_ax, targ_lilac_ax
    global targ_conc_poly_ax, targ_conc_opt_ax, targ_conc_lilac_ax
    global wes_poly_ax, wes_opt_ax, wes_lilac_ax
    global wes_conc_poly_ax, wes_conc_opt_ax, wes_conc_lilac_ax
    global wgs_lilac_ax, wgs_conc_lilac_ax
    
    fig = plt.figure(figsize=(20, 6))
    outer_gs = gridspec.GridSpec(1, 8, width_ratios=[0.2, 0.5, 1, 1, 1, 1, 0.4, 0.4], hspace=0, wspace=0.2)
    
    # Helper to create subplots from a GridSpec
    def create_subplots(spec, n, axes_names):
        gs = gridspec.GridSpecFromSubplotSpec(1, n, subplot_spec=spec, wspace=0)
        return [plt.subplot(gs[i]) for i in range(n)]
    
    # Define subplot structure
    (targ_poly_ax, targ_opt_ax, targ_lilac_ax) = create_subplots(outer_gs[2], 3, ["targ_poly_ax", "targ_opt_ax", "targ_lilac_ax"])
    (targ_conc_poly_ax, targ_conc_opt_ax, targ_conc_lilac_ax) = create_subplots(outer_gs[3], 3, [])
    (wes_poly_ax, wes_opt_ax, wes_lilac_ax) = create_subplots(outer_gs[4], 3, [])
    (wes_conc_poly_ax, wes_conc_opt_ax, wes_conc_lilac_ax) = create_subplots(outer_gs[5], 3, [])
    
    allele_num_ax = plt.subplot(outer_gs[0])
    clinical_assay_ax = plt.subplot(outer_gs[1])
    wgs_lilac_ax = plt.subplot(outer_gs[6])
    wgs_conc_lilac_ax = plt.subplot(outer_gs[7])
    
    all_axes = [
        (allele_num_ax, "Allele"),
        (clinical_assay_ax, "Clinical\nassay"),
        
        (targ_poly_ax, "Polysolver"),
        (targ_opt_ax, "Targeted\nOptitype"),
        (targ_lilac_ax, "Lilac"),
        
        (targ_conc_poly_ax, "Polysolver"),
        (targ_conc_opt_ax, "tNGS match\nOptitype"),
        (targ_conc_lilac_ax, "Lilac"),
        
        (wes_poly_ax, "Polysolver"),
        (wes_opt_ax, "WES\nOptitype"),
        (wes_lilac_ax, "Lilac"),
        
        (wes_conc_poly_ax, "Polysolver"),
        (wes_conc_opt_ax, "WES match\nOptitype"),
        (wes_conc_lilac_ax, "Lilac"),
        
        (wgs_lilac_ax, "WGS\nLilac"),
        (wgs_conc_lilac_ax, "WGS match\nLilac")
    ]
    
    for ax, title in all_axes:
        ax.set_xlim((0.45, 0.55))
        ax.set_ylim(-1, len(df))
        ax.set_yticks(range(len(df)))
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_title(title)
    
    # Define spines to hide per axis group
    def hide_spines(ax, sides):
        for side in sides:
            ax.spines[side].set_visible(False)
    
    spine_instructions = {
        targ_poly_ax: ['right'],
        targ_opt_ax: ['right', 'left'],
        targ_lilac_ax: ['left'],
        targ_conc_poly_ax: ['right'],
        targ_conc_opt_ax: ['right', 'left'],
        targ_conc_lilac_ax: ['left'],
        wes_poly_ax: ['right'],
        wes_opt_ax: ['right', 'left'],
        wes_lilac_ax: ['left'],
        wes_conc_poly_ax: ['right'],
        wes_conc_opt_ax: ['right', 'left'],
        wes_conc_lilac_ax: ['left'],
    }
    
    for ax, sides in spine_instructions.items():
        hide_spines(ax, sides)
    
    # Hide y-axis ticks and labels for all except the first column
    skip_ticks = [clinical_assay_ax, targ_poly_ax, targ_opt_ax, targ_lilac_ax,
                  targ_conc_poly_ax, targ_conc_opt_ax, targ_conc_lilac_ax,
                  wes_poly_ax, wes_opt_ax, wes_lilac_ax,
                  wes_conc_poly_ax, wes_conc_opt_ax, wes_conc_lilac_ax,
                  wgs_lilac_ax, wgs_conc_lilac_ax]
    
    for ax in skip_ticks:
        ax.tick_params(axis="y", labelleft=False, left=False)
    
    return fig, outer_gs


def plot_hla_types(df):
    def match_color(val):
        return 'red' if val == False or val == "False" else 'forestgreen'
    
    for i, row in df.iterrows():
        facecolor = "gray" if i % 2 == 0 else "white"
        bbox_args = dict(facecolor=facecolor, edgecolor='none', alpha=0.3, pad=0)
        
        allele_num_ax.text(0.5, i, row["Allele"], fontsize=10, color='black', ha='center', va='center', bbox=bbox_args)
        clinical_assay_ax.text(0.5, i, row["Clinical assay"], fontsize=10, color='black', ha='center', va='center', bbox=bbox_args)
        
        # Targeted calls
        targ_poly_ax.text(0.5, i, row["Polysolver_targeted"], fontsize=10, color='black', ha='center', va='center', bbox=bbox_args)
        targ_opt_ax.text(0.5, i, row["Optitype_targeted"], fontsize=10, color='black', ha='center', va='center', bbox=bbox_args)
        targ_lilac_ax.text(0.5, i, row["LILAC_targeted"], fontsize=10, color='black', ha='center', va='center', bbox=bbox_args)
        
        # Targeted match results
        for ax, col in [
            (targ_conc_poly_ax, "Polysolver_targeted_match"),
            (targ_conc_opt_ax, "Optitype_targeted_match"),
            (targ_conc_lilac_ax, "LILAC_targeted_match")
        ]:
            ax.text(0.5, i, row[col], fontsize=10, color=match_color(row[col]), ha='center', va='center', bbox=bbox_args)
        
        # WES calls
        wes_poly_ax.text(0.5, i, row["Polysolver_wes"], fontsize=10, color='black', ha='center', va='center', bbox=bbox_args)
        wes_opt_ax.text(0.5, i, row["Optitype_wes"], fontsize=10, color='black', ha='center', va='center', bbox=bbox_args)
        wes_lilac_ax.text(0.5, i, row["LILAC_wes"], fontsize=10, color='black', ha='center', va='center', bbox=bbox_args)
        
        # WES match results
        for ax, col in [
            (wes_conc_poly_ax, "Polysolver_wes_match"),
            (wes_conc_opt_ax, "Optitype_wes_match"),
            (wes_conc_lilac_ax, "LILAC_wes_match")
        ]:
            ax.text(0.5, i, row[col], fontsize=10, color=match_color(row[col]), ha='center', va='center', bbox=bbox_args)
        
        # WGS calls and matches
        wgs_lilac_ax.text(0.5, i, row["LILAC_wgs"], fontsize=10, color='black', ha='center', va='center', bbox=bbox_args)
        wgs_conc_lilac_ax.text(0.5, i, row["LILAC_wgs_match"], fontsize=10, color=match_color(row["LILAC_wgs_match"]), ha='center', va='center', bbox=bbox_args)
    
    allele_num_ax.set_yticklabels(df["Patient"])

for gene in ["A", "B", "C"]:
    df=hla_types[hla_types['Allele'].str.contains(gene)].sort_values(by=["Patient", "Allele"], ascending=False).reset_index(drop=True)
    fig, outer_gs = set_up_grid(df)
    plot_hla_types(df)
    fig.suptitle(f"HLA {gene}/HLA typing concordance across HLA callers and clinical grade HLA typing assay.")
    fig.savefig(f"/groups/wyattgrp/users/amunzur/hla_project/figures/hla_types/{gene}_clinical_assay_concordance.png")
    fig.savefig(f"/groups/wyattgrp/users/amunzur/hla_project/figures/hla_types/{gene}_clinical_assay_concordance.pdf")