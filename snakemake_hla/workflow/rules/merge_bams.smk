# Merges the HLA and non-HLA bams
rule merge_bams:
    input: 
        hla_bam=DIR_bams + "/sorted_hla_filtered/{wildcard}.bam",
        hla_bam_index=DIR_bams + "/sorted_hla_filtered/{wildcard}.bam.bai",
        hg38_bam=DIR_bams + "/sorted_hg38/{wildcard}.bam",
        hg38_bam_index=DIR_bams + "/sorted_hg38/{wildcard}.bam.bai",
        HLA_bed="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/hla_genes.bed"
    output:
        DIR_bams + "/final/{wildcard}.bam"
    shell:
        """
        samtools view -h -L {input.HLA_bed} -U /dev/null {input.hg38_bam} | \
        samtools merge - {input.hla_bam} - > {output}
        """

rule index_merged_bams:
    input:
        DIR_bams + "/final/{wildcard}.bam"
    output:
        DIR_bams + "/final/{wildcard}.bam.bai"
    run:
        shell(
            "samtools index {input}"
        )