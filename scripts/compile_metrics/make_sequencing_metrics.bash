#!/bin/bash

# Calculates following depth metrics for samples:
# 10%, median, 90% depth for HLA-A, HLA-B, HLA-C.
# 10%, median, 90% depth for nonHLA regions in the targeted panel.
# 10%, median, 90% depth for exome (regardless of whether sample is WES)
# Number of reads
# Duplicate fraction
# If any of these is missing, adds a NA to that column. 

out="/groups/wyattgrp/users/amunzur/hla_pipeline/results/metrics/sequencing_quality_metrics_2025_try2.csv"
rm ${out}

DIR_metrics="/groups/wyattgrp/users/amunzur/hla_pipeline/results/metrics"
DIR_bams="/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/sorted"

DIR_median_depth_hla="${DIR_metrics}/depth_hla_MEDIAN"
DIR_median_depth_nonHLAdepth="${DIR_metrics}/depth_panel_nonHLA_MEDIAN"
DIR_median_depth_exome="${DIR_metrics}/depth_exome_MEDIAN"

DIR_markdup="${DIR_metrics}/duplicate_fraction"
DIR_readcounts="${DIR_metrics}/raw_fastq_read_counts"

output="${DIR_metrics}/seq_quality_metrics/panel_test_seq_quality_metrics.csv"

for samp in $(ls ${DIR_bams} | grep "bam$"); do
    samp=$(basename ${samp/.bam/})
    
    # Check for each file and assign values or NA
    if [ -e "${DIR_median_depth_hla}/${samp}_A.txt" ]; then
        A_depth=$(cat ${DIR_median_depth_hla}/${samp}_A.txt | cut -f2 -d ":" | sed 's/^ //; s/ /,/g' | paste -sd ",")
    else
        A_depth="NA,NA,NA"
    fi
    
    if [ -e "${DIR_median_depth_hla}/${samp}_B.txt" ]; then
        B_depth=$(cat ${DIR_median_depth_hla}/${samp}_B.txt | cut -f2 -d ":" | sed 's/^ //; s/ /,/g' | paste -sd ",")
    else
        B_depth="NA,NA,NA"
    fi
    
    if [ -e "${DIR_median_depth_hla}/${samp}_C.txt" ]; then
        C_depth=$(cat ${DIR_median_depth_hla}/${samp}_C.txt | cut -f2 -d ":" | sed 's/^ //; s/ /,/g' | paste -sd ",")
    else
        C_depth="NA,NA,NA"
    fi
    
    if [ -e "${DIR_median_depth_nonHLAdepth}/${samp}.txt" ]; then
        nonhladepth=$(cat ${DIR_median_depth_nonHLAdepth}/${samp}.txt | cut -f2 -d ":" | sed 's/^ //; s/ /,/g' | paste -sd ",")
    else
        nonhladepth="NA,NA,NA"
    fi
    
    if [ -e "${DIR_median_depth_exome}/${samp}.txt" ]; then
        exome_depth=$(cat ${DIR_median_depth_nonHLAdepth}/${samp}.txt | cut -f2 -d ":" | sed 's/^ //; s/ /,/g' | paste -sd ",")
    else
        exome_depth="NA,NA,NA"
    fi
    
    if [ -e "${DIR_readcounts}/${samp}.txt" ]; then
        read_counts=$(cat ${DIR_readcounts}/${samp}.txt | cut -f2)
    else
        read_counts="NA"
    fi
    
    if [ -e "${DIR_markdup}/${samp}.txt" ]; then
        duplicate_fraction=$(cat ${DIR_markdup}/${samp}.txt | cut -f2)
    else
        duplicate_fraction="NA"
    fi
    
    # Append the results to the output file
    echo "${samp},${A_depth},${B_depth},${C_depth},${nonhladepth},${exome_depth},${read_counts},${duplicate_fraction}" >> ${out}
    echo "${samp},${A_depth},${B_depth},${C_depth}"
done

sort -k1 -o ${out} ${out}

cat <(echo -e "Sample_name,HLA-A_10%,HLA-A_median,HLA-A_90%,HLA-B_10%,HLA-B_median,HLA-B_90%,HLA-C_10%,HLA-C_median,HLA-C_90%,nonHLA_10%,nonHLA_median,nonHLA_90%,exome_10%,exome_median,exome_90%,Total_reads,Duplicate_fraction") <(cat ${out}) > "${out}_temp"
mv "${out}_temp"  ${out}