"""
To give an idea about the ctDNA fraction and ploidy disribution
"""

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
import re

path_ct="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_ploidy.tsv"
dir_figures="/groups/wyattgrp/users/amunzur/hla_pipeline/results/figures/coverage_analysis"

df = pd.read_csv(path_ct, sep="\t")
df["Sample_type"]=df["sample"].str.extract(r"(cfDNA|FiT)")
df.loc[df["ctdna"]<0, "ctdna"]=np.nan
df=df[~df["ctdna"].isna()]
df=df.groupby("Sample_type").apply(lambda x: x.sort_values('ctdna')).reset_index(drop=True)


bar_colors={"FiT": "chocolate", "cfDNA": "hotpink"}
df["color"]=df["Sample_type"].map(bar_colors)

# Plotting
fig, ax = plt.subplots(figsize=(8, 3))
ax.bar(df.index, df["ctdna"], color=df["color"])
twin_ax=ax.twinx()
twin_ax.scatter(df.index, df["ploidy"], color="black", s=15)

# AES
ax.spines['top'].set_visible(False)
ax.set_yticks([0, 0.25, 0.50, 0.75, 1])
ax.set_yticklabels(["0", "25", "50", "75", "100"])
ax.set_ylabel("Copy number based ctDNA%\nBAR")
ticklabels=[x.split("_")[0] for x in df["sample"]]
ax.set_xticks(range(0, len(ticklabels), 1))
ax.set_xticklabels(ticklabels, rotation=90)

twin_ax.spines['top'].set_visible(False)
twin_ax.set_yticks([1, 2, 3, 4, 5])
twin_ax.set_yticklabels(["1", "2", "3", "4", "5"])
twin_ax.set_ylabel("Sample ploidy\nSCATTER")

# Add legend
legend_colors = bar_colors.values()
legend_labels= bar_colors.keys()
legend_handles = [plt.Line2D([0], [0], marker="s", color=color, label=label, markersize=4, linestyle='') for marker, color, edgecolor, label in zip(legend_markers, legend_colors, legend_edgecolors, legend_labels)]
ax.legend(handles=legend_handles, loc="best", frameon=False, fontsize = 8, handlelength=2, handletextpad=0.1)

fig.tight_layout()
fig.savefig(os.path.join(dir_figures, f"ctDNA_and_ploidy.png"))
fig.savefig(os.path.join(dir_figures, f"ctDNA_and_ploidy.pdf"))


# Calculate median values
df[df["Sample_type"]=="cfDNA"]["ctdna"].median()
df[df["Sample_type"]=="cfDNA"]["ctdna"].min()
df[df["Sample_type"]=="cfDNA"]["ctdna"].max()

# Calculate median values
df[df["Sample_type"]=="FiT"]["ctdna"].median()
df[df["Sample_type"]=="FiT"]["ctdna"].min()
df[df["Sample_type"]=="FiT"]["ctdna"].max()
