#!/bin/bash

while IFS=$'\t' read -r Patient_ID WBC_name Tumor_name; do
    python /groups/wyattgrp/users/amunzur/toolkit/hla/SLURM_run_DASH.py \
    --normal_fastq /groups/wyattgrp/users/amunzur/hla_project/snakemake_hla/results/data/fq/trimmed_combined/${WBC_name}.fq.gz \
    --tumor_fastq /groups/wyattgrp/users/amunzur/hla_project/snakemake_hla/results/data/fq/trimmed_combined/${Tumor_name}.fq.gz \
    --output_dir /groups/wyattgrp/users/amunzur/hla_project/data/dash_lilac \
    --DIR_batch_scripts /groups/wyattgrp/users/amunzur/hla_project/scripts/batch_scripts/dash_with_lilac 
done < /groups/wyattgrp/users/amunzur/hla_project/snakemake_hla/resources/sample_lists/sample_list_panel.tsv