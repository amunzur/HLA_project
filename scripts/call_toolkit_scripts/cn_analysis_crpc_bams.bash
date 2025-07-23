DIR_hla="/groups/wyattgrp/users/amunzur/hla_project"
DIR_bams_cfDNA="${DIR_hla}/data/crpc_panel_bams/CRPC2022_cfDNA_FiT/"
DIR_vardict_wbc_snps="${DIR_hla}/copynum/wbc_snps_in_snp_grid"
PATH_pon="${DIR_hla}/copynum/pooled_normals/ctdna_neg_evolution_snp_grid_10_samples.cnn"

# Generate coverage files CRPC panel bam files, using cfDNA samples
for file in ${DIR_bams}/*cfDNA*.bam; do
    /groups/wyattgrp/users/amunzur/toolkit/copy_number/SLURM_run_cnvkit_coverage.py \
    --path_bam ${file} \
    --dir_output ${DIR_hla}/copynum/crpc_panel/coverage/snps \
    --path_panel /groups/wyattgrp/users/jbacon/reference/bed_files/all/SNP_spikein_probes.bed \
    --dir_batch_scripts ${DIR_hla}/scripts/batch_scripts/cnvkit_coverage_snps \
    --dir_logs ${DIR_hla}/logs/cnvkit_coverage_snps

# WBC snps using Vardict
for file in /groups/wyattgrp/users/amunzur/hla_project/data/crpc_panel_bams/CRPC2020V2_WBC/*bam \
            /groups/wyattgrp/users/amunzur/hla_project/data/crpc_panel_bams/CRPC2022_WBC/*bam; do
            
            /groups/wyattgrp/users/amunzur/toolkit/variant_calling/SLURM_run_vardict.py \
                --dir_batch_scripts ${DIR_hla}/scripts/batch_scripts/vardict_wbc_snps_in_snp_grid \
                --path_bam ${file} \
                --dir_logs ${DIR_hla}/logs/run_vardict \
                --path_hg38 /groups/wyattgrp/users/jbacon/reference/hg38_emseq/GRCh38_no_alt_plus_hs38d1.fa \
                --threshold_min_vaf 0.40 \
                --min_alt_reads 5 \
                --path_bed /groups/wyattgrp/users/jbacon/reference/bed_files/all/SNP_spikein_probes.bed \
                --dir_output ${DIR_hla}/copynum/crpc_panel/wbc_snps
done

# Run cnvkit fix. 
for file in ${DIR_hla}/copynum/crpc_panel/coverage/snps/*; do

    /groups/wyattgrp/users/amunzur/toolkit/copy_number/SLURM_run_cnvkit_fix.py \
        --path_cnn ${file} \
        --path_pooled_reference_normal "${DIR_hla}/copynum/pooled_normals/ctdna_neg_evolution_snp_grid_10_samples.cnn" \
        --dir_wbc_vcf ${DIR_hla}/copynum/wbc_snps_in_snp_grid \
        --dir_logs ${DIR_hla}/logs/cnvkit_fix \
        --dir_batch_scripts ${DIR_hla}/scripts/batch_scripts/cnvkit_fix \
        --dir_output ${DIR_hla}/copynum/crpc_panel/fix/snps \
        --dir_batch_scripts ${DIR_hla}/scripts/batch_scripts/cnvkit_fix

done

