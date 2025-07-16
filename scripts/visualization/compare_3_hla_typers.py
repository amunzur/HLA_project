import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
from scipy.stats import spearmanr
import subprocess
import matplotlib.gridspec as gridspec
import matplotlib as mpl

"""
Oncoprint style, compares the concordance of the calls made by Optitype, Polysolver and LILAC on targeted data.
"""

def calculate_median_depth_from_polysolver_sam(path_sam):
    """
    Given path to polsolver bam calculates median depth.
    """
    if os.path.exists(path_sam):
        depth_path=path_sam.replace(".sam", ".depth")
        subprocess.run(["/home/amunzur/anaconda3/envs/snakemake/bin/samtools", "index", path_sam])
        subprocess.run(["/home/amunzur/anaconda3/envs/snakemake/bin/samtools", "depth", "a -d0 -g", path_sam, ">", depth_path])
        depth_df=pd.read_csv(depth_path, sep="\t", header=None)
        depth=depth_df[2].median()
    else:
        depth=np.nan
    
    return(depth)

def return_hla_types_for_samples(wbc_samples_list):
    """
    Given a list of samples returns their HLA type.
    They must all be WBC samples.
    """
    dir_polysolver_main="/groups/wyattgrp/users/amunzur/hla_pipeline/results/polysolver/hla_types/"
    
    sample_dict={}
    for wbc in wbc_samples_list:
        dir_polysolver=os.path.join(dir_polysolver_main, wbc)
        path_winners=os.path.join(dir_polysolver, "winners.hla.txt")
        if os.path.exists(path_winners):
            winners=pd.read_csv(path_winners, sep="\t", header=None).T
            winners.columns = winners.iloc[0]  # Set the first row as header
            winners = winners[1:].reset_index(drop=True)
            alleles=set(pd.concat([winners["HLA-A"], winners["HLA-B"], winners["HLA-C"]]).tolist())
            sample_dict[wbc]=alleles
    
    return(sample_dict)

path_info="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test.tsv"
sample_info=pd.read_csv(path_info, sep="\t")
targeted_samples=sample_info["WBC_name"].drop_duplicates().tolist()
targeted_types=return_hla_types_for_samples(wbc_samples)

path_info_wes="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/WES_samples.tsv"
sample_info_wes=pd.read_csv(path_info_wes, sep="\t")
wes_samples=sample_info_wes["WBC_name"].drop_duplicates().tolist()

extra_wes=["GU-16-005_WBC_31Oct2016_WES", "GU-17-046_WBC_20Mar2017_WES", "GU-17-105_WBC_01Jun2017_WES", "GU-17-295_WBC_08Nov2017_WES", "GU-18-161_WBC_11Apr2018_WES", "GU-18-296_WBC_09Jul2018_WES"]
wes_samples.extend(extra_wes)

wes_types=return_hla_types_for_samples(wes_samples)

# Figure set up
fig = plt.figure(figsize=(6, 3))
outer_gs = gridspec.GridSpec(1, 4, width_ratios=[1, 1, 1, 1], hspace = 0, wspace = 0)

hla_a_ax=plt.subplot(outer_gs[1])
hla_b_ax=plt.subplot(outer_gs[2], sharey=hla_a_ax)
hla_c_ax=plt.subplot(outer_gs[3], sharey=hla_a_ax)
legend_ax=plt.subplot(outer_gs[0])

hla_a_ax.set_xlabel("HLA A")
hla_a_ax.xaxis.set_label_position("top")

hla_b_ax.set_xlabel("HLA B")
hla_b_ax.xaxis.set_label_position("top")

hla_c_ax.set_xlabel("HLA C")
hla_c_ax.xaxis.set_label_position("top")

zygosity_color={1: "deepskyblue", 2: "crimson"} # 1 homo 2 het

# Plotting function
sample_counter=0
names_list=[]
for wbc in targeted_types.keys():
    wes_name=wbc+"_WES"
    if wes_name in wes_types:
        names_list.append(wes_name)
        for gene, ax in zip(["hla_a", "hla_b", "hla_c"], [hla_a_ax, hla_b_ax, hla_c_ax]):
            t_al=set(sorted([x for x in targeted_types[wbc] if gene in x]))
            wes_al=set(sorted([x for x in wes_types[wes_name] if gene in x]))
            
            t_al_len=len(t_al)
            wes_al_len=len(wes_al)
            
            for xpos, zygosity in enumerate([t_al_len, wes_al_len]):
                color=zygosity_color[zygosity]
                ax.scatter(xpos, sample_counter, edgecolor="None", s=60, marker="s", color=color)
            
            ax.set_xlim((-0.7, 1.7))
            ax.spines[['top', 'right', "bottom", "left"]].set_visible(False)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["Targeted", "WES"], rotation=90)
        
        sample_counter+=1

hla_a_ax.set_yticks(range(0, sample_counter))
hla_a_ax.set_yticklabels([x.split("_WBC")[0] for x in names_list])

hla_b_ax.tick_params(axis="y", labelleft=False, left=False)
hla_c_ax.tick_params(axis="y", labelleft=False, left=False)

# Add legend
legend_colors = zygosity_color.values()
legend_labels = ["Hom.", "Het."]
legend_handles = [plt.Line2D([0], [0], marker='s', color=color, label=label, markersize=4, linestyle='') for color, label in zip(legend_colors, legend_labels)]
legend_ax.legend(handles=legend_handles, loc="lower left", frameon=False, fontsize = 8, handlelength=2, handletextpad=0.1)

legend_ax.spines[['top', 'right', "bottom", "left"]].set_visible(False)
legend_ax.set_yticks([])
legend_ax.set_yticklabels([])
legend_ax.set_xticks([])
legend_ax.set_xticklabels([])

outer_gs.tight_layout(fig)
fig.suptitle("")
fig.savefig(os.path.join(dir_figures, f"wes_vs_targeted_zygosity_heatmap.png"))
fig.savefig(os.path.join(dir_figures, f"wes_vs_targeted_zygosity_heatmap.pdf"))


# FIGURE 2. 
# Second heatmap. 
# Figure set up
fig = plt.figure(figsize=(6, 5))
outer_gs=gridspec.GridSpec(2, 1, height_ratios=[1, 0.25], hspace = 0.2, wspace = 0)
heatmap_gs = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer_gs[0], width_ratios=[1, 1, 1], hspace = 0, wspace = 0.1)

hla_a_ax=plt.subplot(heatmap_gs[0])
hla_b_ax=plt.subplot(heatmap_gs[1], sharey=hla_a_ax)
hla_c_ax=plt.subplot(heatmap_gs[2], sharey=hla_a_ax)
legend_ax=plt.subplot(outer_gs[1])

hla_a_ax.set_xlabel("HLA A")
hla_a_ax.xaxis.set_label_position("top")

hla_b_ax.set_xlabel("HLA B")
hla_b_ax.xaxis.set_label_position("top")

hla_c_ax.set_xlabel("HLA C")
hla_c_ax.xaxis.set_label_position("top")


sample_counter=0
common_pts=[]
for wbc in targeted_types.keys():
    wes_name=wbc+"_WES"
    if wes_name in wes_types:
        pt=wes_name.split("_WBC")[0]
        common_pts.append(pt)

for i, pt in enumerate(common_pts):
    targeted_als={key: value for key, value in targeted_types.items() if pt in key}
    wes_als={key: value for key, value in wes_types.items() if pt in key}
    
    first_match = None
    for key, value in targeted_als.items():
        if pt in key:
            first_match = {key: value}
            break
    targeted_als=first_match
    
    for hla_gene, ax in zip(["hla_a", "hla_b", "hla_c"], [hla_a_ax, hla_b_ax, hla_c_ax]):
        targeted_als_gene=sorted(list({v for value in targeted_als.values() for v in value if hla_gene in v}))
        wes_als_gene=sorted(list({v for value in wes_als.values() for v in value if hla_gene in v}))
        
        if len(targeted_als_gene)==len(wes_als_gene):
                if len(targeted_als_gene) == 1:  # If only one element in the list
                    allele_max = 0
                    allele_offset_max = -0.5
                else:
                    allele_max = 1
                    allele_offset_max = 0.5
                for allele, allele_offset in zip([0, allele_max], [-0.5, allele_offset_max]):
                    targeted_al=targeted_als_gene[allele].split("_")[2:]
                    wes_al=wes_als_gene[allele].split("_")[2:]
                    
                    targeted_al += ["99"] * (4 - len(targeted_al))
                    wes_al += ["99"] * (4 - len(wes_al))
                    for num, offset in zip([0, 1, 2, 3], [-0.3, -0.1, 0.1, 0.3]):
                        targeted_pos=targeted_al[num]
                        wes_pos=wes_al[num]
                        if targeted_pos==wes_pos and targeted_pos!="99" and wes_pos!="99":
                            color="limegreen"
                        elif targeted_pos=="99" and wes_pos!="99":
                            color="gray"
                        elif wes_pos=="99" and targeted_pos!="99":
                            color="blue"
                        elif targeted_pos!=wes_pos and targeted_pos!="99" and wes_pos!="99":
                            color="red"
                        ax.scatter(allele_offset+offset, i, marker="s", color=color, s=100, edgecolor="None")
                        
                        if len(targeted_als_gene) == 1:
                            ax.scatter([x+0.5 for x in [-0.3, -0.1, 0.1, 0.3]], np.repeat(i, 4), marker="x", color="black", s=20, edgecolor="None")
        else:
            # Determine the allele the singular one is most similar to
            homo_original=[targeted_als_gene[0] if len(targeted_als_gene)==1 else wes_als_gene[0]]
            het_original=targeted_als_gene if len(targeted_als_gene)>1 else wes_als_gene
            
            homo=[a.split("_")[2:] for a in homo_original]
            het=[a.split("_")[2:] for a in het_original]
            
            [het[i].extend(["99"] * (4 - len(het[i]))) for i in range(len(het))]
            [homo[i].extend(["99"] * (4 - len(homo[i]))) for i in range(len(homo))]
            
            similarity_scores=[0, 0]
            
            for pos in [0, 1, 2, 3]:
                for het_allele_n in [0, 1]:
                    homo_pos=homo[0][pos]
                    het_pos=het[het_allele_n][pos]
                    if homo_pos==het_pos:
                        similarity_scores[het_allele_n]+=1
            
            most_similar_allele_n=similarity_scores.index(max(similarity_scores))
            most_similar_allele_het=het[most_similar_allele_n]
            
            # Now determine which one was targeted and which one was wes
            targeted_al=None
            wes_al=None
            
            targeted_al = homo[0] if any("_".join(homo[0]) in x for x in targeted_als_gene) else most_similar_allele_het
            wes_al = homo if homo is None else most_similar_allele_het
            
            # Plotting
            allele_max = 0
            allele_offset_max = -0.5
            for allele, allele_offset in zip([0, allele_max], [-0.5, allele_offset_max]):
                for num, offset in zip([0, 1, 2, 3], [-0.3, -0.1, 0.1, 0.3]):
                    targeted_pos=targeted_al[num]
                    wes_pos=wes_al[num]
                    if targeted_pos==wes_pos and targeted_pos!="99" and wes_pos!="99":
                        color="limegreen"
                    elif targeted_pos=="99" and wes_pos!="99":
                        color="gray"
                    elif wes_pos=="99" and targeted_pos!="99":
                        color="blue"
                    elif targeted_pos!=wes_pos and targeted_pos!="99" and wes_pos!="99":
                        color="red"
                    ax.scatter(allele_offset+offset, i, marker="s", color=color, s=100, edgecolor="None")
            
            # Hashed boxes in zygosity mismatch
            allele_pos = 0.5
            for offset in [-0.3, -0.1, 0.1, 0.3]:
                ax.scatter(allele_pos+offset, i, marker="s", color="white", s=100, edgecolor="Black")

                    # if len(targeted_als_gene) == 1:
                    #     ax.scatter([x+0.5 for x in [-0.3, -0.1, 0.1, 0.3]], np.repeat(i, 4), marker="x", color="black", s=20, edgecolor="None")

for ax in [hla_a_ax, hla_b_ax, hla_c_ax]:
    ax.spines[['top', 'right', "bottom", "left"]].set_visible(False)
    ax.set_xticks([-0.5, 0.5])
    ax.set_xticklabels(["Allele1", "Allele2"])
    ax.tick_params(axis="x", bottom=False)

hla_a_ax.tick_params(axis="y", left=False)
hla_b_ax.tick_params(axis="y", labelleft=False, left=False)
hla_c_ax.tick_params(axis="y", labelleft=False, left=False)

hla_a_ax.set_yticks(range(0, len(common_pts)))
hla_a_ax.set_yticklabels(common_pts)

# Add legend
legend_colors = ["limegreen", "gray", "blue", "red", "black", "white"]
legend_edgecolors=["None", "None", "None", "None", "black", "black"]
legend_labels= ["Exact match", "Extra field in WES", "Extra field in targeted", "Mismatch in same field", "Homozygous in both targeted and WES", "Zygosity mismatch between targeted and WES"]
legend_markers=['s', 's', 's', 's', 'x', 's']
legend_handles = [plt.Line2D([0], [0], marker=marker, color=color, markeredgecolor= edgecolor, label=label, markersize=7, linestyle='') for marker, color, edgecolor, label in zip(legend_markers, legend_colors, legend_edgecolors, legend_labels)]
legend_ax.legend(handles=legend_handles, loc="lower left", frameon=False, fontsize = 8, handlelength=2, handletextpad=0.1)
legend_ax.spines[['top', 'right', "bottom", "left"]].set_visible(False)
legend_ax.tick_params(axis="both", labelleft=False, left=False, labelbottom=False, bottom=False)
legend_ax.set_title("")

outer_gs.tight_layout(fig)
# fig.suptitle("Only showing alleles with same zygosity in targeted vs WES")
fig.savefig(os.path.join(dir_figures, f"HLA_type_OP.png"))
fig.savefig(os.path.join(dir_figures, f"HLA_type_OP.pdf"))

