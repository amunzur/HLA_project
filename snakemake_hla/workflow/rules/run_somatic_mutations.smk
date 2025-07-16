def get_wbc_path(wildcard, DIR_bams):
    paired_samples = pd.read_csv("/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test.tsv", sep="\t")
    paired_samples = paired_samples[["WBC_name", "Tumor_name"]]
    
    mask = paired_samples["Tumor_name"] == wildcard
    wbcname = paired_samples["WBC_name"][mask].tolist()[0]
    wbcpath=DIR_bams + "/sorted/" + wbcname + ".bam"
    
    return(wbcpath)

rule run_VarDict_somatic:
    input:
        cfDNA=DIR_bams + "/sorted/{wildcard}.bam",
        wbc=lambda wildcards: get_wbc_path(wildcards.wildcard, "/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam"),
    output:
        temp(DIR_results+ "/variant_calling_somatic_raw/Vardict/{wildcard}.vcf"),
    conda:
        "/groups/wyattgrp/users/amunzur/lu_chip/workflow/envs/vardict_env.yaml"
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
        THRESHOLD_VarFreq=0.005,
        sample_name_cfDNA="{wildcard}",
        sample_name_wbc=lambda wildcards: get_wbc_name(wildcards.wildcard),
    threads: 12
    shell:
        """
        /home/amunzur/VarDictJava/build/install/VarDict/bin/VarDict \
            -G {params.PATH_hg38} \
            -b '{input.cfDNA}|{input.wbc}' \
            -f {params.THRESHOLD_VarFreq} \
            -N {params.sample_name_cfDNA} \
            -k 1 \
            -c 1 \
            -S 2 \
            -E 3 \
            -g 4 \
            -r 4 \
            --nosv \
            -th {threads} \
            {params.PATH_bed} | \
            ~/VarDictJava/build/install/VarDict/bin/testsomatic.R | \
            ~/VarDictJava/build/install/VarDict/bin/var2vcf_paired.pl \
        -N '{params.sample_name_cfDNA}|{params.sample_name_wbc}' -f {params.THRESHOLD_VarFreq} > {output}
        """

rule run_mutect2_somatic:
    input:
        cfDNA=DIR_bams + "/sorted/{wildcard}.bam",
        wbc=lambda wildcards: get_wbc_path(wildcards.wildcard, "/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam"),
    output:
        vcf=temp(DIR_results + "/variant_calling_somatic_raw/Mutect2/{wildcard}.vcf.gz"),
        stats=DIR_results + "/variant_calling_somatic_raw/Mutect2/{wildcard}.vcf.gz.stats",
    conda:
        "/groups/wyattgrp/users/amunzur/lu_chip/workflow/envs/chip_variantcalling.yaml"
    params:
        PATH_hg38=PATH_hg38,
        PATH_bed=PATH_bed,
        sample_name_wbc=lambda wildcards: get_wbc_name(wildcards.wildcard),
    threads: 12
    shell:
        "/home/amunzur/gatk-4.2.0.0/gatk Mutect2 \
        --reference {params.PATH_hg38} \
        --input {input.cfDNA} \
        --input {input.wbc} \
        --normal-sample {params.sample_name_wbc} \
        --output {output.vcf} \
        --force-active true \
        --initial-tumor-lod 0 \
        --tumor-lod-to-emit 0 \
        --intervals {params.PATH_bed}"

rule zip_vcf_files_vardict_somatic:
    input:
        DIR_results + "/variant_calling_somatic_raw/Vardict/{wildcard}.vcf",
    output:
        temp(DIR_results + "/variant_calling_somatic_raw/Vardict/{wildcard}.vcf.gz"),
    conda:
        "/groups/wyattgrp/users/amunzur/lu_chip/workflow/envs/bcftools.yaml"
    shell:
        "bgzip -c {input} > {output}"

rule index_vcf_files_somatic:
    input:
        DIR_results + "/variant_calling_somatic_raw/{variant_caller}/{wildcard}.vcf.gz",
    output:
        temp(DIR_results + "/variant_calling_somatic_raw/{variant_caller}/{wildcard}.vcf.gz.tbi"),
    conda:
        "/groups/wyattgrp/users/amunzur/lu_chip/workflow/envs/bcftools.yaml"
    shell:
        "tabix -p vcf {input}"


rule sort_vcf_somatic:
    input:
        vcf=DIR_results + "/variant_calling_somatic_raw/{variant_caller}/{wildcard}.vcf.gz",
        index=DIR_results + "/variant_calling_somatic_raw/{variant_caller}/{wildcard}.vcf.gz.tbi",
    output:
        temp(DIR_results + "/variant_calling_somatic_sorted/{variant_caller}/{wildcard}.vcf.gz"),
    conda:
        "/groups/wyattgrp/users/amunzur/lu_chip/workflow/envs/bcftools.yaml"
    shell:
        "bcftools sort {input.vcf} -Oz -o {output}"


rule index_sorted_vcf_files_somatic:
    input:
        DIR_results + "/variant_calling_somatic_sorted/{variant_caller}/{wildcard}.vcf.gz",
    output:
        temp(DIR_results + "/variant_calling_somatic_sorted/{variant_caller}/{wildcard}.vcf.gz.tbi"),
    conda:
        "/groups/wyattgrp/users/amunzur/lu_chip/workflow/envs/bcftools.yaml"
    shell:
        "tabix -p vcf {input}"

rule normalize_variants_somatic:
    input:
        vcf=DIR_results + "/variant_calling_somatic_sorted/{variant_caller}/{wildcard}.vcf.gz",
        index=DIR_results + "/variant_calling_somatic_sorted/{variant_caller}/{wildcard}.vcf.gz.tbi"
    output:
        temp(DIR_results + "/variant_calling_somatic_normalized/{variant_caller}/{wildcard}.vcf.gz"),
    params:
        PATH_hg38_dict="/groups/wyattgrp/reference/hg38/hg38.fa",
    conda:
        "/groups/wyattgrp/users/amunzur/lu_chip/workflow/envs/bcftools.yaml"
    shell:
        "bcftools norm \
        {input.vcf} \
        -m-both \
        -f {params} \
        -o {output}"

rule decompose_blocksubstitutions_somatic:
    input:
        DIR_results + "/variant_calling_somatic_normalized/{variant_caller}/{wildcard}.vcf.gz",
    output:
        DIR_results + "/variant_calling_somatic/{variant_caller}/{wildcard}.vcf.gz",
    conda:
        "/groups/wyattgrp/users/amunzur/lu_chip/workflow/envs/chip_variantcalling.yaml"
    shell:
        "vt decompose_blocksub {input} -o {output}"