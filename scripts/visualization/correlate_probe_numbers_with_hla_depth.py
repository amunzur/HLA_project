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
import matplotlib as mpl
from sklearn.preprocessing import MinMaxScaler
import matplotlib.colors as mcolors
import matplotlib.cm as cm

"""
Tries to answer the question, do samples with more probes aligning to their HLA alleles have more sequencing depth?
Ideally we would see no positive correlation here.

Uses "/groups/wyattgrp/users/amunzur/hla_pipeline/results/metrics/compiled_depth/hla_allele_specific_depth.tsv" file.
Script calculate_hla_allele_specific_depth.py makes this file ^.
Saves to "/groups/wyattgrp/users/amunzur/hla_pipeline/results/figures/coverage_analysis"
"""

path_DNA_concentration="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/Library_DNA_concentrations.csv"
path_ancestry="/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/crpc_bams/ancestry/somalier-ancestry.somalier-ancestry_copy.tsv"
dir_hla_fasta="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/lohhla_fasta"
path_probe_sequences="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/hla_panel/HLA_Panel_design_FOR_ROCHE_known_probes_HLA_genes_formatted.fa"
dir_alignments="/groups/wyattgrp/users/amunzur/hla_pipeline/results/metrics/probe_hlafasta_alignment"
path_allele_depth="/groups/wyattgrp/users/amunzur/hla_pipeline/results/metrics/compiled_depth/hla_allele_specific_depth.tsv"
probe_seq=pd.read_csv(path_probe_sequences)
dir_figures="/groups/wyattgrp/users/amunzur/hla_pipeline/results/figures/coverage_analysis"

color_dict={"cfDNA": "orangered", "FiT": "deepskyblue", "WBC": "limegreen"}
color_dict_ancestry={"SAS": "#1f77b4", "EAS": "#ff7f0e", "EUR": "#2ca02c", "AMR": "#d62728"}

def align_hla_probes_to_patient_fasta(path_fasta, path_probe_sequences, dir_alignments, sample_name):
    """
    Given path to a hla_fasta file aligns probe sequences to it. Fixes mates, sorts the bam file, and indexes.
    """
    # Step 2: Align probe sequences to the reference fasta using BWA MEM
    path_aligned_sam = os.path.join(dir_alignments, f"{sample_name}_aligned.sam")
    bwa_mem_cmd = f"bwa mem {path_fasta} {path_probe_sequences} > {path_aligned_sam}"
    subprocess.run(bwa_mem_cmd, shell=True, check=True)

def calculate_probe_matches_to_gene(path_aligned_sam, hla_gene):
    """
    Given the name of an HLA gene calculates how many probes aligned to the fasta file with a user defined number of exact matches.
    """
    hla_gene=hla_gene.lower()
    
    sam=pd.read_csv(path_aligned_sam, sep="\t", comment="@", header=None, on_bad_lines="skip")
    sam_gene=sam[sam[2].str.contains(hla_gene)]
    sam_gene_filtered = sam_gene[sam_gene[5].str.contains(r'10[0-9]M|11[0-9]M|120M')]
    n_probes=sam_gene_filtered[0].drop_duplicates().shape[0]
    
    return(n_probes)

def locate_WBC_hla_fasta(wbc_sample, dir_hla_fasta):
    """
    Given a WBC sample name, locates its fasta file with the sequences of all HLA alleles.
    """
    path_fasta=os.path.join(dir_hla_fasta, wbc_sample, "hla_fasta.fa")
    return(path_fasta)

def find_matching_wbc(sample_name):
    """
    Returns the sample name of the WBC match.
    """
    path_info="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test.tsv"
    sample_info=pd.read_csv(path_info, sep="\t")
    sample_info_subset=sample_info[sample_info["Tumor_name"]==sample_name]
    if sample_info_subset.shape[0]!=1:
        raise ValueError(f"{sample_name} didn't match to exactly one WBC sample.")
    else:
        wbc_name=sample_info[sample_info["Tumor_name"]==sample_name]["WBC_name"].values[0]
    return(wbc_name)

def return_polysolver_winners(wbc_name):
    path_winners=os.path.join("/groups/wyattgrp/users/amunzur/hla_pipeline/results/polysolver/hla_types", wbc_name, "winners.hla.txt")
    winners=pd.read_csv(path_winners, sep="\t", header=None)
    winners.columns=["hla_gene", "allele1", "allele2"]
    winners["hla_gene"]=winners["hla_gene"].replace({"HLA-A": "hla_a", "HLA-B": "hla_b", "HLA-C": "hla_c"})
    
    return(winners)

def plot1_nprobe_depth_scatter(gene_df, color_dict, ax, add_legend = True, subset_sample_type=None):
    """
    Main plotting function.
    """
    gene_df["Sample_type"]=gene_df["Sample_name"].str.extract(r"(cfDNA|WBC|FiT)")
    gene_df["color"]=gene_df["Sample_type"].map(color_dict)
    ax.scatter(gene_df["nprobes"], gene_df["depth_norm"], edgecolor="None", s=10, color=gene_df["color"])
    
    ax.set_xlabel("Number of probes")
    ax.set_ylabel("Normalized depth")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # ax.set_aspect('equal')
    # ax.set_ylim((0, 3600))
    # ax.set_yticks([0, 1000, 2000, 3000])
    # ax.set_yticklabels(["0", "1000", "2000", "3000"])
    
    if add_legend:
        if subset_sample_type is None:
            legend_colors = color_dict.values()
            legend_labels = color_dict.keys()
        else:
            legend_colors = [color_dict[subset_sample_type]]
            legend_labels = [subset_sample_type]
        
        legend_handles = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=4, linestyle='') for color, label in zip(legend_colors, legend_labels)]
        ax.legend(handles=legend_handles, loc="lower left", frameon=False, fontsize = 8, handlelength=2, handletextpad=0.1)
    
    return(ax)

def plot2_nprobe_depth_scatter_with_DNA_conc(plotting_df, scatter_color, ax, add_colorbar = True):
    """
    Second plotting function.
    """
    scaler = MinMaxScaler(feature_range=(0.2, 1.0))
    plotting_df["alpha"]=scaler.fit_transform(plotting_df["DNA concentration"].values.reshape(-1, 1))
    plotting_df["color"]=scatter_color
    ax.scatter(plotting_df["nprobes"], plotting_df["depth_norm"], s=10, color=plotting_df["color"], alpha=plotting_df["alpha"], edgecolor="black")
    
    ax.set_xlabel("Number of probes")
    ax.set_ylabel("Normalized depth")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    if add_colorbar:
        norm = mcolors.Normalize(vmin=plotting_df["DNA concentration"].min(), vmax=plotting_df["DNA concentration"].max())
        
        custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", ["white", scatter_color])
        sm = cm.ScalarMappable(cmap=custom_cmap, norm=norm)
        sm.set_array([])  # No actual data needed
        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label("DNA Concentration")
    return(ax)

def plot3_nprobe_depth_scatter_with_ancestry(plotting_df, color_col, ax, add_legend = True, legend_dict=None):
    """
    Third plotting function.
    """
    scaler = MinMaxScaler(feature_range=(0.2, 1.0))
    ax.scatter(plotting_df["nprobes"], plotting_df["depth_norm"], s=10, color=plotting_df[color_col], edgecolor="None")
    
    ax.set_xlabel("Number of probes")
    ax.set_ylabel("Normalized depth")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    if add_legend and legend_dict is not None:
        legend_colors = legend_dict.values()
        legend_labels = legend_dict.keys()
        legend_handles = [plt.Line2D([0], [0], marker='o', color=color, label=label, markersize=4, linestyle='') for color, label in zip(legend_colors, legend_labels)]
        ax.legend(handles=legend_handles, loc="upper left", frameon=False, fontsize = 8, handlelength=2, handletextpad=0.1)
    return(ax)

path_sample_information="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test.tsv"
sample_info=pd.read_csv(path_sample_information, sep="\t")
all_samples = pd.concat([sample_info["WBC_name"].drop_duplicates(), sample_info["Tumor_name"]]).dropna().unique()

sample_dict={}
for sample_name in all_samples:
    print(f"Working on {sample_name}.")
    if "WBC" not in sample_name:
        wbc_name = find_matching_wbc(sample_name)
    else:
        wbc_name=sample_name
    
    # Locate patient fasta
    path_fasta=locate_WBC_hla_fasta(wbc_name, dir_hla_fasta)
    
    # Align probes to patient fasta
    path_aligned_sam = os.path.join(dir_alignments, f"{sample_name}_aligned.sam")
    if not os.path.exists(path_aligned_sam):
        align_hla_probes_to_patient_fasta(path_fasta, path_probe_sequences, dir_alignments, sample_name)
    
    # Determine how many probes aligned with 100M for each allele
    hla_types=return_polysolver_winners(wbc_name)
    hla_types=pd.concat([hla_types['allele1'], hla_types['allele2']], ignore_index=True).tolist()
    
    hla_dict={}
    for hla_type in hla_types:        
        n_probes=int(calculate_probe_matches_to_gene(path_aligned_sam, hla_type))
        hla_dict[hla_type]=n_probes
    
    sample_dict[sample_name]=hla_dict            

data = []
for sample, alleles in sample_dict.items():
    for allele, count in alleles.items():
        data.append([sample, allele, count])

# Create a DataFrame from the list
nprobes_df = pd.DataFrame(data, columns=["Sample_name", "allele", "nprobes"])

# Now combine this with depth information
depth_df=pd.read_csv(path_allele_depth, sep="\t")
combined_df = depth_df.merge(nprobes_df, how = "inner")

# Load DNA concentration information
dna_conc_df=pd.read_csv(path_DNA_concentration).rename(columns={"Sample name": "Sample_name"})[["Sample_name", "DNA concentration"]]

# Load ancestry information
ancestry_df=pd.read_csv(path_ancestry, sep="\t")[["#sample_id", "predicted_ancestry"]]
ancestry_df.columns=["Sample_name", "Ancestry"]
ancestry_df["Sample_name"] = ancestry_df["Sample_name"].str.replace(r"[_-]WBC.*", "", regex=True)
ancestry_df["Ancestry_color"]=ancestry_df["Ancestry"].map(color_dict_ancestry)

#########################################################
# PLOTTING 1. SCATTER, NPROBES VS ALLELE DEPTH. All sample types shown
fig = plt.figure(figsize=(10, 4))
fig.suptitle("Correlating nprobes with allele depth - ideally we won't observe a positive correlation.")
outer_gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1], hspace = 0.5, wspace = 0.5) # outer most gs with 3 cols

for i, hla_gene in enumerate(["hla_a", "hla_b", "hla_c"]):
    gene_df=combined_df[combined_df["allele"].str.contains(hla_gene)]
    ax=plt.subplot(outer_gs[i])
    ax=plot1_nprobe_depth_scatter(gene_df, color_dict, ax)
    
    ax_title=hla_gene.upper().replace("_", " ")
    ax.set_title(ax_title)

outer_gs.tight_layout(fig)
fig.savefig(os.path.join(dir_figures, f"allele_depth_vs_nprobes.png"))
fig.savefig(os.path.join(dir_figures, f"allele_depth_vs_nprobes.pdf"))

#########################################################
# PLOTTING 2. SCATTER, NPROBES VS ALLELE DEPTH. Individual plots for sample types. Determine alpha based on DNA concentration.
fig = plt.figure(figsize=(10, 8))
fig.suptitle("Correlating nprobes with allele depth - scatter transparency determined by DNA concentration.")
outer_gs = gridspec.GridSpec(3, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1, 1], hspace = 0.5, wspace = 0.5) # outer most gs with 3 cols

for i, sample_type in enumerate(['WBC', 'cfDNA', 'FiT']):
    for j, hla_gene in enumerate(["hla_a", "hla_b", "hla_c"]):
        scatter_color=color_dict[sample_type]
        plotting_df=combined_df[(combined_df["allele"].str.contains(hla_gene)) & (combined_df["Sample_name"].str.contains(sample_type))]
        plotting_df=plotting_df.merge(dna_conc_df)
        ax=plt.subplot(outer_gs[i, j])
        ax=plot2_nprobe_depth_scatter_with_DNA_conc(plotting_df, scatter_color, ax, add_colorbar = True)
        
        ax_title=hla_gene.upper().replace("_", " ")+"\n"+sample_type
        ax.set_title(ax_title)
        
        if i!=2:
            ax.set_xlabel("") 
        if j!=0:
            ax.set_ylabel("")
        if j==1:
            ax.set_title(ax_title)

outer_gs.tight_layout(fig)
fig.savefig(os.path.join(dir_figures, f"allele_depth_vs_nprobes_DNA_conc.png"))
fig.savefig(os.path.join(dir_figures, f"allele_depth_vs_nprobes_DNA_conc.pdf"))

#########################################################
# PLOTTING 3. SCATTER, NPROBES VS ALLELE DEPTH. Individual plots for sample types. Coloring samples based on ancestry.
fig = plt.figure(figsize=(10, 8))
fig.suptitle("Correlating nprobes with allele depth - coloring scatter based on ancestry.")
outer_gs = gridspec.GridSpec(3, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1, 1], hspace = 0.5, wspace = 0.5) # outer most gs with 3 cols

for i, sample_type in enumerate(['WBC', 'cfDNA', 'FiT']):
    for j, hla_gene in enumerate(["hla_a", "hla_b", "hla_c"]):
        if i==0 and j==0:
            add_legend=True
        else:
            add_legend=False
        scatter_color=color_dict[sample_type]
        plotting_df=combined_df[(combined_df["allele"].str.contains(hla_gene)) & (combined_df["Sample_name"].str.contains(sample_type))]
        plotting_df["Sample_name"]=plotting_df["Sample_name"].str.replace(r"[_-]?(WBC|FiT|cfDNA).*", "", regex=True)
        plotting_df=plotting_df.merge(ancestry_df)
        
        ax=plt.subplot(outer_gs[i, j])
        ax=plot3_nprobe_depth_scatter_with_ancestry(plotting_df, color_col="Ancestry_color", ax=ax, add_legend = add_legend, legend_dict=color_dict_ancestry)
        
        ax_title=hla_gene.upper().replace("_", " ")+"\n"+sample_type
        ax.set_title(ax_title)
        
        if i!=2:
            ax.set_xlabel("") 
        if j!=0:
            ax.set_ylabel("")
        if j==1:
            ax.set_title(ax_title)

outer_gs.tight_layout(fig)
fig.savefig(os.path.join(dir_figures, f"allele_depth_vs_nprobes_ancestry.png"))
fig.savefig(os.path.join(dir_figures, f"allele_depth_vs_nprobes_ancestry.pdf"))

# PLOTTING 3. SIMPLE BAR CHART OF N PROBES PER ALLELE IN FORM OF FANCY ONCOPRINT STAYLA
# Normalize nprobes values to [0, 1] range
norm = plt.Normalize(vmin=nprobes_df["nprobes"].min(), vmax=nprobes_df["nprobes"].max())
cmap = plt.cm.Reds

nprobes_df_wbc=nprobes_df[nprobes_df["Sample_name"].str.contains("WBC")]
nprobes_df_wbc["Sample_name"] = nprobes_df_wbc["Sample_name"].str.extract(r"(.*WBC)")

nprobes_df_wbc["color"] = nprobes_df_wbc["nprobes"].apply(lambda x: cmap(norm(x)))
nprobes_df_wbc=nprobes_df_wbc.sort_values(by="nprobes", ascending=False)

fig = plt.figure(figsize=(3, 7))
fig.suptitle("probe-allele al. 100M. x:homozygous.")
outer_gs = gridspec.GridSpec(2, 1, height_ratios=[1, 0.01], hspace = 0, wspace = 0)
heatmap_gs = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=outer_gs[0], width_ratios=[3, 1, 1, 1], hspace = 0, wspace = 0.1)

hla_a_ax=plt.subplot(heatmap_gs[1])
hla_b_ax=plt.subplot(heatmap_gs[2])
hla_c_ax=plt.subplot(heatmap_gs[3])
legend_ax=plt.subplot(outer_gs[1])

hla_a_ax.set_xlabel("HLA A")
hla_a_ax.xaxis.set_label_position("top")

hla_b_ax.set_xlabel("HLA B")
hla_b_ax.xaxis.set_label_position("top")

hla_c_ax.set_xlabel("HLA C")
hla_c_ax.xaxis.set_label_position("top")

for i, group in nprobes_df_wbc.groupby("Sample_name"):
    for gene, ax in zip(["hla_a", "hla_b", "hla_c"], [hla_a_ax, hla_b_ax, hla_c_ax]):
        gene_subset=group[group["allele"].str.contains(gene)].reset_index(drop=True)
                
        # Plotting
        for j, row in gene_subset.iterrows():
            ax.scatter(j, i, edgecolor="None", s=60, marker="s", color=row["color"])
        
        if gene_subset.shape[0]==1:
            ax.scatter(1, i, edgecolor="None", s=60, marker="x", color="black")
    
        # AES
        ax.set_xlim((-0.5, 1.5))
        ax.spines[['top', 'right', 'left', 'bottom']].set_visible(False)
        ax.set_xticks([])
        ax.set_xticklabels([])
        
        if gene!="hla_a":
            ax.set_yticks([])
            ax.set_yticklabels([])

# Color bar
# Add color bar to legend_ax
cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=legend_ax, orientation="horizontal")
cbar.set_label("Number of Probes")

outer_gs.tight_layout(fig)
fig.savefig(os.path.join(dir_figures, f"nprobes_per_sample_per_allele_heatmap.png"))
fig.savefig(os.path.join(dir_figures, f"nprobes_per_sample_per_allele_heatmap.pdf"))


