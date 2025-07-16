rule run_ANNOVAR_somatic:
    input:
        DIR_results + "/variant_calling_somatic/{variant_caller}/{wildcard}.vcf.gz",
    output:
        DIR_results + "/data/annovar_outputs_somatic/{variant_caller}/{wildcard}.hg38_multianno.vcf",
        DIR_results + "/data/annovar_outputs_somatic/{variant_caller}/{wildcard}.hg38_multianno.txt",
    params:
        DIR_results + "/data/annovar_outputs_somatic/{variant_caller}/{wildcard}",
    shell:
        "perl /groups/wyattgrp/users/amunzur/software/annovar/table_annovar.pl {input} /groups/wyattgrp/users/amunzur/software/annovar/humandb \
        -vcfinput \
        -buildver hg38 \
        -out {params} \
        -remove \
        -protocol refGene,cosmic97_coding,avsnp150,clinvar_20221231,gnomad40_exome \
        -operation g,f,f,f,f \
        -nastring ."

rule vcfToTable_Vardict:
    input:
        DIR_results + "/data/annovar_outputs_somatic/Vardict/{wildcard}.hg38_multianno.vcf",
    output:
        DIR_results + "/data/variant_tables/somatic/Vardict/{wildcard}.tsv",
    shell:
        """
        /home/amunzur/gatk-4.2.0.0/gatk VariantsToTable \
        -V {input} \
        -F CHROM -F POS -F REF -F ALT -F TYPE -F FILTER -F STATUS -GF DP -GF VD -GF AF -GF ALD -GF RD -GF SBF -GF ODDRATIO \
        -F Func.refGene -F Gene.refGene -F ExonicFunc.refGene -F AAChange.refGene -F cosmic97_coding -F avsnp150 -F gnomad40_exome_AF -F CLNALLELEID -F CLNSIG \
        -O {output}
        """

rule vcfToTable_Mutect2_somatic:
    input:
        DIR_results + "/data/annovar_outputs_somatic/Mutect2/{wildcard}.hg38_multianno.vcf",
    output:
        DIR_results + "/data/variant_tables/somatic/Mutect2/{wildcard}.tsv",
    shell:
        """
        /home/amunzur/gatk-4.2.0.0/gatk VariantsToTable \
            -V {input} \
            -F CHROM -F POS -F REF -F ALT -F TYPE -F EVENTLENGTH -GF SB \
            -F Func.refGene -F Gene.refGene -F ExonicFunc.refGene -F AAChange.refGene -F AF -F cosmic97_coding -F avsnp150 -F gnomad40_exome_AF -F CLNALLELEID -F CLNSIG \
            -O {output}
        """