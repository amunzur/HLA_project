#!/bin/bash

# Running
while IFS=$'\t' read -r Patient_ID WBC_name Tumor_name; do
    /groups/wyattgrp/users/amunzur/toolkit/SLURM_run_optitype.py \
    --dir_fastqs /groups/wyattgrp/users/amunzur/hla_pipeline/results/data/fq/trimmed \
    --sample_name ${WBC_name} \
    --dir_output_main /groups/wyattgrp/users/amunzur/hla_pipeline/results/optitype \
    --dir_batch_scripts /groups/wyattgrp/users/amunzur/hla_pipeline/workflow/batch_scripts/optitype \
    --dir_logs /groups/wyattgrp/users/amunzur/hla_pipeline/results/logs_slurm/optitype \
    --threads 16 \
    --memory 64G \
    --timelimit 6:00:00 \
    --jobname_keyword OPT
done < /groups/wyattgrp/users/amunzur/hla_pipeline/resources/all_samples.tsv