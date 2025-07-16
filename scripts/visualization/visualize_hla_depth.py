# Plots the WES coverage against HLA coverage
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')  # Use the non-interactive backend (Agg)
import matplotlib.pyplot as plt
import os
from sklearn.metrics import r2_score

df = pd.read_csv("/groups/wyattgrp/users/amunzur/hla_pipeline/results/metrics/sequencing_quality_metrics.txt")
df['Sample_type'] = df['Sample_name'].str.extract(r'(cfDNA|Tissue|WBC)')
df['Sample_type'] = df['Sample_type'].fillna('Unknown')
df = df[["Sample_name", "Sample_type", "HLA-A_median", "HLA-B_median", "HLA-C_median", "WES_median"]]
# df.columns = ['Gene', 'Sample_name', 'Sample_Type', 'mismatch', 'Median_hla', 'WES_median', 'Sample_type']

mask = df["Sample_type"] != "Unknown"
df = df[mask]

dir_plots = "/groups/wyattgrp/users/amunzur/hla_pipeline/results/figures/depth"
fig, axs = plt.subplots(1, 3, figsize=(15, 5)) # Create a figure with 3 subplots in one row

# Define colors for each sample type
colors = {'cfDNA': 'blue', 'Tissue': 'green', 'WBC': 'orange'}

# Plot wes median vs hla A median and save
def make_HLAvsWES_depth_plots(ax, df, gene):
    wes_median = df["WES_median"]
    hla_median = df[gene]
    sample_type = df['Sample_type']
    
    # Plot the scatter points with colors based on sample type
    ax.scatter(wes_median, hla_median, c=sample_type.map(colors))
    ax.set_xlabel("WES Median")
    ax.set_ylabel(gene)
    # ax.set_title("WES Median vs " + gene)
    
    # Add best fit line and equation to plot
    slope, intercept = np.polyfit(wes_median, hla_median, 1)
    best_fit_line = slope * wes_median + intercept
    ax.plot(wes_median, best_fit_line, linestyle="--", color="red")
    equation = "y = {:.2f}x + {:.2f}".format(slope, intercept)
    r2 = r2_score(wes_median, hla_median)
    r_squared = "R^2 = {:.2f}".format(r2)
    ax.text(max(wes_median) * 0.75, max(hla_median) * 0.70, equation + "\n" + r_squared, color="red")
    
    # Create a separate legend for the sample type points
    handles = [plt.scatter([], [], c=color, label=label) for label, color in colors.items()]
    ax.legend(handles=handles, title="Sample Type", loc='upper left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def make_HLA_depth_plots(df, path_save): 
    df_subset = df[["HLA-A_median", "HLA-B_median", "HLA-C_median", "WES_median"]]
    fig, ax = plt.subplots()
    
    x_values = {
        "HLA-A_median": [0.9, 1.1],
        "HLA-B_median": [1.9, 2.1],
        "HLA-C_median": [2.9, 3.1],
        "WES_median": [3.9, 4.1]}

    for index, row in df_subset.iterrows():
        for a in ["HLA-A_median", "HLA-B_median", "HLA-C_median", "WES_median"]:
            x_range = x_values[a]
            x = np.random.uniform(x_range[0], x_range[1], size=1)
            y = row[a]
            ax.scatter(x, y, color='gray', alpha=0.7)

    boxplot_positions = [1, 2, 3, 4]
    
    for idx, depth_col in enumerate([df_subset["HLA-A_median"], df_subset["HLA-B_median"], df_subset["HLA-C_median"], df_subset["WES_median"]]):
        ax.boxplot(depth_col, positions=[boxplot_positions[idx]], widths=0.3, showfliers=False, medianprops=dict(linewidth=2.0, color="black"), boxprops=dict(linewidth=2.0), capprops=dict(linewidth=2.0), whiskerprops=dict(linewidth=2.0))

    ax.set_xlim(right=4.5)
    ax.set_xlim(left=-0.005)
    ax.set_xticks(boxplot_positions)
    ax.set_xticklabels(["HLA-A", "HLA-B", "HLA-C", "WES"])
    ax.set_title('Showing HLA and WES median depth')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig(path_save)


make_HLAvsWES_depth_plots(axs[0], df, "HLA-A_median")
make_HLAvsWES_depth_plots(axs[1], df, "HLA-B_median")
make_HLAvsWES_depth_plots(axs[2], df, "HLA-C_median")

plt.tight_layout()
plt.savefig(os.path.join(dir_plots, "test.png"))

make_HLA_depth_plots(df, os.path.join(dir_plots, "test.png"))