DIR_hla="/groups/wyattgrp/users/amunzur/hla_project"



# Generate coverage files from WGS and WES, using the cfDNA samples
cut -f5,7 ${DIR_hla}/sample_lists/panel_wes_wgs_samples.tsv | tail -n +2 | \
tr '\t' '\n' | grep -v '^$' | while read -r sample_name; do

    /groups/wyattgrp/users/amunzur/toolkit/copy_number/SLURM_run_cnvkit_coverage.py \
    --path_bam ${DIR_hla}/snakemake_hla/results/data/bam/sorted/${sample_name}.bam \
    --dir_output ${DIR_hla}/copynum/coverage/snps \
    --path_panel /groups/wyattgrp/users/jbacon/reference/bed_files/all/SNP_spikein_probes.bed \
    --dir_batch_scripts ${DIR_hla}/scripts/batch_scripts/cnvkit_coverage_snps \
    --dir_logs ${DIR_hla}/logs/cnvkit_coverage_snps

done

# Run vardict for germline snps across the snp grid, using WBC samples. First WES samples.
for type in WES WGS; do
    if [ "$type" == "WES" ]; then
        col=4
        bam_dir=${DIR_hla}/snakemake_hla/results/data/bam/sorted
    else
        col=6
        bam_dir=${DIR_hla}/data/bams_without_alt_contigs/alignments
    fi
    
    tail -n +2 ${DIR_hla}/sample_lists/panel_wes_wgs_samples.tsv | cut -f4 | grep -v '^$' | \
    while read -r wbc_name; do
        /groups/wyattgrp/users/amunzur/toolkit/variant_calling/SLURM_run_vardict.py \
            --dir_batch_scripts ${DIR_hla}/scripts/batch_scripts/vardict_wbc_snps_in_snp_grid \
            --path_bam --path_bam ${bam_dir}/${wbc_name}.bam \
            --dir_logs ${DIR_hla}/logs/run_vardict \
            --path_hg38 /groups/wyattgrp/users/jbacon/reference/hg38_emseq/GRCh38_no_alt_plus_hs38d1.fa \
            --threshold_min_vaf 0.40 \
            --min_alt_reads 5 \
            --path_bed /groups/wyattgrp/users/jbacon/reference/bed_files/all/SNP_spikein_probes.bed \
            --dir_output ${DIR_hla}/copynum/wbc_snps_in_snp_grid
    done

# Before running cnvkit fix we generate a panel of normals using ctdna negative bams
for bam in ${DIR_hla}/copynum/ctdna_neg_bams/*.bam; do
    /groups/wyattgrp/users/amunzur/toolkit/copy_number/SLURM_run_cnvkit_coverage.py \
        --path_bam "${bam}" \
        --dir_output "${DIR_hla}/copynum/coverage/snps" \
        --path_panel /groups/wyattgrp/users/jbacon/reference/bed_files/all/SNP_spikein_probes.bed \
        --dir_batch_scripts "${DIR_hla}/scripts/batch_scripts/cnvkit_coverage_snps" \
        --dir_logs "${DIR_hla}/logs/cnvkit_coverage_snps"
done

/groups/wyattgrp/users/amunzur/toolkit/copy_number/SLURM_generate_pooled_normal.py \
--path_cnns "${DIR_hla}/copynum/pooled_normals/ctdna_neg_path.txt" \
--path_hg38 /groups/wyattgrp/users/jbacon/reference/matti_hg38/hg38_masked.fa \
--dir_logs "${DIR_hla}/logs/make_pooled_normal" \
--path_output "${DIR_hla}/copynum/pooled_normals/ctdna_neg_evolution_snp_grid_10_samples.cnn" \
--path_sbatch "${DIR_hla}/scripts/batch_scripts/pooled_normal/pooled_norm_ctdna_neg_evolution_snp_grid_10_samples.batch" \

# Run cnvkit fix on WGS and WES, using pooled normal. 
cut -f5,7 ${DIR_hla}/sample_lists/panel_wes_wgs_samples.tsv | tail -n +2 | \
tr '\t' '\n' | grep -v '^$' | while read -r sample_name; do

    /groups/wyattgrp/users/amunzur/toolkit/copy_number/SLURM_run_cnvkit_fix.py \
        --path_cnn "${DIR_hla}/copynum/coverage/snps/${sample_name}.cnn" \
        --path_pooled_reference_normal "${DIR_hla}/copynum/pooled_normals/ctdna_neg_evolution_snp_grid_10_samples.cnn" \
        --dir_wbc_vcf ${DIR_hla}/copynum/wbc_snps_in_snp_grid \
        --dir_logs ${DIR_hla}/logs/cnvkit_fix \
        --dir_batch_scripts ${DIR_hla}/scripts/batch_scripts/cnvkit_fix \
        --dir_output ${DIR_hla}/copynum/fix/snps \
        --dir_batch_scripts ${DIR_hla}/scripts/batch_scripts/cnvkit_fix

done

