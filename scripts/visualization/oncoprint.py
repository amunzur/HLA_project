import matplotlib.gridspec as gridspec
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
from matplotlib.patches import Patch

"""
This script generates a multi-panel oncoprint-style visualization summarizing HLA and non-HLA genomic features 
for a cohort of tumor samples, using cfDNA and FiT sequencing data.

The output figure includes:
1. Tumor fraction (ctDNA%) per sample.
2. Sample type (cfDNA vs. FiT).
3. HLA gene zygosity (homozygous vs. heterozygous).
4. HLA loss of heterozygosity (LOH) status.
5. Non-HLA somatic mutation presence and type.
6. Variant allele frequency (VAF) of mutations (log scale).

Data sources:
- Sample metadata
- HLA LOH output (LOHHLA)
- Curated somatic mutation calls
- Ploidy and ctDNA estimates

Legends are automatically generated and placed in the top-right corner. 
Output is saved as both PNG and PDF.
"""

path_sample_list="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test.tsv"
path_loh="/groups/wyattgrp/users/amunzur/hla_pipeline/results/loh/lohhla_dash_compiled.csv"
path_ctdna_mutations="/groups/wyattgrp/users/amunzur/hla_pipeline/results/variant_calling/SOMATIC_SSCS2_curated.csv"
path_ctdna_frac_and_ploidy="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_ploidy.tsv"

mut_dict = {
    "Missense": '#79B443',
    "Stop gain": '#BD4398',
    "Frameshift InDel": '#FFC907',
    "Splicing": "darkorange"}

renaming_dict={
    "missense": 'Missense',
    "stop_gain": 'Stop gain',
    "frameshift": 'Frameshift InDel',
    "splicing": "Splicing"   
}

sample_info=pd.read_csv(path_sample_list, sep="\t")

fig = plt.figure(figsize=(8, 7))
outer_gs = gridspec.GridSpec(5, 2, height_ratios=[0.3, 0.06, 0.24, 0.24, 1], width_ratios=[1, 0.2], hspace=0.05, wspace=0)

ax_tc=ax=plt.subplot(outer_gs[0, 0])
ax_sample_type=plt.subplot(outer_gs[1, 0], sharex=ax_tc)
ax_zygosity=plt.subplot(outer_gs[2, 0], sharex=ax_tc)
ax_loh=plt.subplot(outer_gs[3, 0], sharex=ax_tc)
ax_muts=plt.subplot(outer_gs[4, 0], sharex=ax_tc)
ax_vaf=plt.subplot(outer_gs[4, 1], sharey=ax_muts)

# STEP 1. Plot the ctDNA fraction
ctdna_df=pd.read_csv(path_ctdna_frac_and_ploidy, sep="\t")
ctdna_df["ctdna"]=ctdna_df["ctdna"]*100
ctdna_df["sample_type"]=ctdna_df["sample"].str.extract("(cfDNA|FiT)")
ctdna_df["Patient"]=ctdna_df["sample"].str.replace(r'_(cfDNA|FiT).*$', '', regex=True)
ctdna_df=ctdna_df.sort_values(by=["sample_type", "ctdna"]).reset_index(drop=True)
ctdna_df=ctdna_df[ctdna_df["ctdna"]>0].reset_index(drop=True)
ax_tc.bar(ctdna_df.index, ctdna_df["ctdna"], color="black")

ax_tc.set_xlim((-1, 42))
ax_tc.set_ylim((0, 100))
ax_tc.set_yticks([0, 50, 100])
ax_tc.set_yticklabels(['0', '50', '100'])
ax_tc.set_ylabel("Tumor %", rotation=0, ha="right", va="center")
ax_tc.spines[['top', 'right']].set_visible(False)
ax_tc.tick_params(bottom=False, labelbottom=False)

samples_enumerated_df=ctdna_df.reset_index()[["sample", "index"]]
samples_enumerated_df["Patient"] = samples_enumerated_df["sample"].str.replace(r'_(cfDNA|FiT).*$', '', regex=True)

# STEP 2. Plot sample type
sample_type_color={"cfDNA": "orangered", "FiT": "deepskyblue"}
sample_type_df=samples_enumerated_df.copy()
sample_type_df["type"]=sample_type_df["sample"].str.extract("(cfDNA|FiT)")
sample_type_df["color"]=sample_type_df["type"].map(sample_type_color)
ax_sample_type.bar(sample_type_df["index"], np.repeat(1, len(sample_type_df)), color=sample_type_df["color"])

ax_sample_type.set_ylabel("Sample type", rotation=0, ha="right", va="center")
ax_sample_type.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
ax_sample_type.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)

# STEP 3. HLA gene homozygous/heterozygous state for each HLA gene (binary)
homozygosity_dict={True: "limegreen", False: "orange"}
hla_gene_ypos={"hla_a":2, "hla_b":1, "hla_c":0}
loh_main=pd.read_csv(path_loh)
zygosity_df=loh_main[["Patient", "Gene", "Homozygous"]]
zygosity_df=zygosity_df.merge(samples_enumerated_df).reset_index(drop=True)[["Patient", "Gene", 'Homozygous', "index"]].drop_duplicates()
zygosity_df["color"]=zygosity_df["Homozygous"].map(homozygosity_dict)

bottom=0
for hla_gene in ["hla_c", "hla_b", "hla_a"]:
    subset_zygosity_df=zygosity_df[zygosity_df["Gene"]==hla_gene]
    ax_zygosity.bar(
        subset_zygosity_df["index"], 
        np.repeat(1, len(subset_zygosity_df)), 
        bottom=bottom, 
        color=subset_zygosity_df["color"], 
        edgecolor="white")
    bottom+=1

ax_zygosity.set_ylabel("Zygosity", rotation=0, ha="right", va="center")
ax_zygosity.set_yticks([0.5, 1.5, 2.5])
ax_zygosity.set_yticklabels(["HLA-C", "HLA-B", "HLA-A"])
ax_zygosity.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
ax_zygosity.tick_params(left=False, bottom=False, labelbottom=False)

# STEP 4. LOH status
loh=loh_main.copy()
loh=loh[["Patient", "Sample", "Gene", "Deletion DASH", "Homozygous"]].drop_duplicates().rename(columns={"Sample": "sample"})
loh=loh.merge(samples_enumerated_df[["sample", "index"]])
loh.loc[loh["Homozygous"]==True, "Deletion DASH"]="homozygous"
deletion_dict={True: "Red", False: "blue", "homozygous": "white", np.nan: "grey"} # np.nan means fail
loh["LOH color"]=loh["Deletion DASH"].map(deletion_dict)

bottom=0
for hla_gene in ["hla_c", "hla_b", "hla_a"]:
    subset_loh=loh[loh["Gene"]==hla_gene]
    ax_loh.bar(subset_loh["index"], np.repeat(1, len(subset_loh)), bottom=bottom,color=subset_loh["LOH color"], edgecolor="white")
    
    # Mark homozygous genes
    homo_df=subset_loh[subset_loh["Deletion DASH"]=="homozygous"]
    ax_loh.scatter(homo_df["index"], np.repeat(bottom+0.5, len(homo_df)), marker="x", s=20, edgecolor="none", color="black")
    bottom+=1

ax_loh.set_ylabel("LOH", rotation=0, ha="right", va="center")
ax_loh.set_yticks([0.5, 1.5, 2.5])
ax_loh.set_yticklabels(["HLA-C", "HLA-B", "HLA-A"])
ax_loh.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
ax_loh.tick_params(left=False, bottom=False, labelbottom=False)

# STEP 5 - Plot nonHLA mutations
muts_main=pd.read_csv(path_ctdna_mutations)
muts_main["Sample_type"]=muts_main["Sample_name_t"].str.extract("(cfDNA|FiT)")
muts=muts_main.drop(muts_main[(muts_main["Sample_type"] == "FiT") & (muts_main["VAF_t"] < 0.05)].index)
muts=muts[["Patient_id", "Sample_name_t", "Sample_type", "Gene", "VAF_t", "Effects", "Protein_annotation"]]
muts["VAF_t"]=muts["VAF_t"]*100

# Plot the gray grid
genes_order=muts["Gene"].drop_duplicates().reset_index(drop=True).reset_index()
n_genes = len(genes_order)
n_samples = len(samples_enumerated_df)

x = np.repeat(np.arange(n_samples), n_genes)      # sample index (j)
y = np.tile(np.arange(n_genes), n_samples)        # gene index (i)
ax_muts.bar(x, 1, bottom=y, color="whitesmoke", edgecolor="white")

# Plot the mutations on the grid
muts=muts.merge(genes_order).rename(columns={"index": "gene_pos"})
muts=muts.merge(samples_enumerated_df, left_on="Sample_name_t", right_on="sample").rename(columns={"index": "sample_pos"})
muts["Effects_renamed"]=muts["Effects"].map(renaming_dict)
muts["color"]=muts["Effects_renamed"].map(mut_dict)

for (sample_pos, gene_pos), group in muts.groupby(["sample_pos", "gene_pos"]):
    n = len(group)
    center = gene_pos + 0.5
    
    if n == 1:
        ax_muts.scatter(sample_pos, center, s=10, color=group.iloc[0]["color"], marker='s')
    else:
        # Vertical offsets within the bar (bar spans 1 unit from gene_pos to gene_pos+1)
        offsets = np.linspace(center - 0.3, center + 0.3, n)  # stay inside the bar
        for offset, (_, row) in zip(offsets, group.iterrows()):
            ax_muts.scatter(sample_pos, offset, s=10, color=row["color"], marker='s')

ax_muts.set_yticks(genes_order["index"]+0.5)
ax_muts.set_yticklabels(genes_order["Gene"], fontstyle="italic")
ax_muts.spines[['top', 'right', 'left', 'bottom']].set_visible(False)
ax_muts.tick_params(left=False)
ax_muts.set_xticks(samples_enumerated_df["index"])

# Set x ticklabels
samples_enumerated_df["tick label"] = (samples_enumerated_df["Patient"] + " " + samples_enumerated_df["sample"].str.extract(r"(Baseline|OnTreatment)").squeeze())
samples_enumerated_df["tick label"] = samples_enumerated_df["tick label"].fillna(samples_enumerated_df["Patient"])
ax_muts.set_xticklabels(samples_enumerated_df["tick label"], rotation=90)

# STEP 6. Plot mutation vafs on the right panel 
muts["vaf log"]=np.log10(muts["VAF_t"])

jitter_strength = 0.08
y_jitter = muts["gene_pos"] + 0.5 + np.random.uniform(-jitter_strength, jitter_strength, size=len(muts))
ax_vaf.scatter(muts["vaf log"], y_jitter, marker="o", s=3, color="black")

ax_vaf.set_xlim((np.log10(0.4), np.log10(35)))
ax_vaf.set_xticks([np.log10(0.5), np.log10(2), np.log10(10), np.log10(30)])
ax_vaf.set_xticklabels(["0.5", "2", "10", "30"])
ax_vaf.tick_params(left=False, labelleft=False)
ax_vaf.spines[['top', 'right']].set_visible(False)
ax_vaf.set_xlabel("Mutation VAF% (log)")

# Add legends to the top right corner

# 1. Sample Type Legend
sample_legend_handles = [
    Patch(facecolor=sample_type_color["cfDNA"], label="cfDNA"),
    Patch(facecolor=sample_type_color["FiT"], label="FiT")
]

# 2. Zygosity Legend
zygosity_legend_handles = [
    Patch(facecolor=homozygosity_dict[True], label="Homozygous"),
    Patch(facecolor=homozygosity_dict[False], label="Heterozygous")
]

# 3. Mutation Effect Legend
mutation_legend_handles = [
    Patch(facecolor=color, label=label) for label, color in mut_dict.items()
]

# 4. LOH Status Legend
loh_legend_handles = [
    Patch(facecolor="red", label="LOH (+)"),
    Patch(facecolor="blue", label="LOH (–)"),
    Patch(facecolor="white", edgecolor="black", label="Homozygous"),
    Patch(facecolor="grey", label="Failed LOH call")
]
spacer = Patch(facecolor='none', edgecolor='none', label='')  # invisible spacer
# Combine all legends
legend_handles = (
    sample_legend_handles +
    [spacer] +
    zygosity_legend_handles +
    [spacer] +
    mutation_legend_handles +
    [spacer] +
    loh_legend_handles
)
# Add legend to top-right of the figure
fig.legend(handles=legend_handles, loc="upper right", frameon=False, ncol=1, fontsize=8)









outer_gs.tight_layout(fig)
fig.savefig("/groups/wyattgrp/users/amunzur/hla_pipeline/results/figures/oncoprints/op.png")
fig.savefig("/groups/wyattgrp/users/amunzur/hla_pipeline/results/figures/oncoprints/op.pdf")


# STEP 4. HLA LOH status for each gene (binary)



