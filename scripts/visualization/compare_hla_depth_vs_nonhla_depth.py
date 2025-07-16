import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
from matplotlib.gridspec import GridSpec
from scipy.stats import spearmanr
import subprocess
import matplotlib.gridspec as gridspec

"""
Tries to answer the question, does having high HLA depth deplete depth from other regions?
Expecting to see a correlation here, to make the point that there is no depletion.

Uses these two files: 
path_hla_depth="/groups/wyattgrp/users/amunzur/hla_pipeline/results/metrics/compiled_depth/hla_allele_specific_depth.tsv"
path_nonhla_depth="/groups/wyattgrp/users/amunzur/hla_pipeline/results/metrics/compiled_depth/nonhla_depth.tsv"

These files generate these two files:
/groups/wyattgrp/users/amunzur/hla_pipeline/workflow/depth_scripts/calculate_hla_allele_specific_depth.py
/groups/wyattgrp/users/amunzur/hla_pipeline/workflow/depth_scripts/calculate_non_hla_depth.py

Saves figure to: 
/groups/wyattgrp/users/amunzur/hla_pipeline/results/figures/coverage_analysis/hla_vs_nonhla_depth_scatter.png

"""

color_dict={"cfDNA": "orangered", "FiT": "deepskyblue", "WBC": "limegreen"}
color_dict_ancestry={"SAS": "#1f77b4", "EAS": "#ff7f0e", "EUR": "#2ca02c", "AMR": "#d62728"}

path_DNA_concentration="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/Library_DNA_concentrations.csv"
path_ancestry="/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/crpc_bams/ancestry/somalier-ancestry.somalier-ancestry_copy.tsv"
path_hla_depth="/groups/wyattgrp/users/amunzur/hla_pipeline/results/metrics/compiled_depth/hla_allele_specific_depth.tsv"
path_nonhla_depth="/groups/wyattgrp/users/amunzur/hla_pipeline/results/metrics/compiled_depth/nonhla_depth.tsv"
dir_figures="/groups/wyattgrp/users/amunzur/hla_pipeline/results/figures/coverage_analysis"
color_dict={"cfDNA": "orangered", "FiT": "deepskyblue", "WBC": "limegreen"}

hla_depth=pd.read_csv(path_hla_depth, sep = "\t")
hla_depth["gene"]=hla_depth["allele"].str.extract(r"(hla_a|hla_b|hla_c)")
hla_depth.columns=["sample_name", "allele", "hla_depth", "hla_depth_norm", "gene"]

nonhla_depth=pd.read_csv(path_nonhla_depth, sep = "\t")
nonhla_depth.columns=["sample_name", "nonhla_depth", "nonhla_depth_norm"]
combined=hla_depth.merge(nonhla_depth)

combined["sample_type"]=combined["sample_name"].str.extract(r"(cfDNA|WBC|FiT)")
combined["color"]=combined["sample_type"].map(color_dict)
combined["Patient_id"]=combined["sample_name"].str.replace(r"[_-]?(WBC|FiT|cfDNA).*", "", regex=True)

# Load DNA concentration information
dna_conc_df=pd.read_csv(path_DNA_concentration).rename(columns={"Sample name": "sample_name"})[["sample_name", "DNA concentration"]]

# Load ancestry information
ancestry_df=pd.read_csv(path_ancestry, sep="\t")[["#sample_id", "predicted_ancestry"]]
ancestry_df.columns=["Patient_id", "Ancestry"]
ancestry_df["Patient_id"] = ancestry_df["Patient_id"].str.replace(r"[_-]WBC.*", "", regex=True)
ancestry_df["Ancestry_color"]=ancestry_df["Ancestry"].map(color_dict_ancestry)

combined=combined.merge(ancestry_df).merge(dna_conc_df)

fig = plt.figure(figsize=(8, 6))
# fig.suptitle("Correlating nprobes with allele depth - ideally we won't observe a positive correlation.")
outer_gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1], hspace = 0.5, wspace = 0.5)


def STEP1_set_up_axes(i, gene):
    gene_gs = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[1, 1], subplot_spec=outer_gs[i], wspace=0, hspace = 0.5)
    allele1_ax=plt.subplot(gene_gs[0])
    allele2_ax=plt.subplot(gene_gs[1], sharex=allele1_ax)
    
    ax_title=gene.upper().replace("_", " ")
    allele1_ax.set_title(f"{ax_title}\nAllele 1")
    allele2_ax.set_title("Allele 2")
    allele2_ax.set_xlabel("Normalized non-HLA depth")
    
    ax_list=[allele1_ax, allele2_ax]
    return(ax_list)

def STEP2_manage_ax_aesthetics(ax_list, gene):
    for ax in ax_list:
        ax.plot([0, 4000], [0, 6000], linestyle='--', color='red', label='y=x')
        if gene == "hla_a":
            ax.set_ylabel("Normalized HLA depth")
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # ax.set_aspect('equal')
        ax.set_xlim((0, 6000))
        ax.set_xticks([0, 3000, 6000])
        ax.set_xticklabels(["0", "3K", "6K"]) 
        ax.set_ylim((0, 6000))
        ax.set_yticks([0, 3000, 6000])
        ax.set_yticklabels(["0", "3K", "6K"]) 
        ax.set_aspect("equal")
    
    return(ax_list)

def STEP3_add_legend(ax_list, color_dict):
    for ax in ax_list:
        legend_colors = color_dict.values()
        legend_labels = color_dict.keys()
        legend_handles = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=4, linestyle='') for color, label in zip(legend_colors, legend_labels)]
        ax.legend(handles=legend_handles, loc="lower right", frameon=False, fontsize = 8, handlelength=2, handletextpad=0.1)
    
    return(ax_list)

def main(df, color_col, legend_dict, path_figure, add_spearmans_r=True):
    for i, (gene, group) in enumerate(df.groupby("gene")):
        
        # For spearmans later
        hla_vals_list=[[], []]
        nonhla_vals_list=[[], []]
        
        ax_list=STEP1_set_up_axes(i, gene)
        ax_list=STEP2_manage_ax_aesthetics(ax_list, gene)
        ax_list=STEP3_add_legend(ax_list, legend_dict)
        
        for j, pt_group in group.groupby("sample_name"):
            if pt_group.shape[0]>2:
                pt_group=pt_group[0:2]
            for k, row in pt_group.reset_index(drop=True).iterrows():
                hla_vals_list[k].append(row["hla_depth_norm"])
                nonhla_vals_list[k].append(row["nonhla_depth_norm"])
                ax_list[k].scatter(row["nonhla_depth_norm"], row["hla_depth_norm"], edgecolor="None", s=10, color=row[color_col])
        
        if add_spearmans_r:
            rs_1, p_value = spearmanr(hla_vals_list[0], nonhla_vals_list[0])
            rs_2, p_value = spearmanr(hla_vals_list[1], nonhla_vals_list[1])
            
            ax_list[0].text(0.05, 0.95, f"Spearman's r:\n{rs_1:.2f}", transform=ax_list[0].transAxes, ha='left', va='top', fontsize=10, bbox=dict(boxstyle="round", facecolor="white", edgecolor="none"))
            ax_list[1].text(0.05, 0.95, f"Spearman's r:\n{rs_2:.2f}", transform=ax_list[1].transAxes, ha='left', va='top', fontsize=10, bbox=dict(boxstyle="round", facecolor="white", edgecolor="none"))
        
        # Save
        outer_gs.tight_layout(fig)
        fig.savefig(path_figure)

main(combined, color_col="color", legend_dict=color_dict, path_figure=os.path.join(dir_figures, f"hla_vs_nonhla_depth.png"), add_spearmans_r=True)
main(combined, color_col="color", legend_dict=color_dict, path_figure=os.path.join(dir_figures, f"hla_vs_nonhla_depth.pdf"), add_spearmans_r=True)

main(combined, color_col="color", legend_dict=color_dict, path_figure=os.path.join(dir_figures, f"hla_vs_nonhla_depth.png"), add_spearmans_r=True)
main(combined, color_col="color", legend_dict=color_dict, path_figure=os.path.join(dir_figures, f"hla_vs_nonhla_depth.pdf"), add_spearmans_r=True)

cfDNA_entries=combined[combined["sample_name"].str.contains("cfDNA")]
main(cfDNA_entries, color_col="Ancestry_color", legend_dict=color_dict_ancestry, path_figure=os.path.join(dir_figures, f"hla_vs_nonhla_depth_ancestry_cfDNA.png"), add_spearmans_r=False)
main(cfDNA_entries, color_col="Ancestry_color", legend_dict=color_dict_ancestry, path_figure=os.path.join(dir_figures, f"hla_vs_nonhla_depth_ancestry_cfDNA.pdf"), add_spearmans_r=False)

WBC_entries=combined[combined["sample_name"].str.contains("WBC")]
main(WBC_entries, color_col="Ancestry_color", legend_dict=color_dict_ancestry, path_figure=os.path.join(dir_figures, f"hla_vs_nonhla_depth_ancestry_WBC.png"), add_spearmans_r=False)
main(WBC_entries, color_col="Ancestry_color", legend_dict=color_dict_ancestry, path_figure=os.path.join(dir_figures, f"hla_vs_nonhla_depth_ancestry_WBC.pdf"), add_spearmans_r=False)

FiT_entries=combined[combined["sample_name"].str.contains("FiT")]
main(FiT_entries, color_col="Ancestry_color", legend_dict=color_dict_ancestry, path_figure=os.path.join(dir_figures, f"hla_vs_nonhla_depth_ancestry_FiT.png"), add_spearmans_r=False)
main(FiT_entries, color_col="Ancestry_color", legend_dict=color_dict_ancestry, path_figure=os.path.join(dir_figures, f"hla_vs_nonhla_depth_ancestry_FiT.pdf"), add_spearmans_r=False)


