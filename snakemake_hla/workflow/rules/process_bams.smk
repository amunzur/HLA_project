rule mapBAM:
    input:
        R1=DIR_fastq + "/trimmed/{wildcard}_1.fq.gz",
        R2=DIR_fastq + "/trimmed/{wildcard}_2.fq.gz",
    params:
        PATH_hg38=PATH_hg38,
    output:
        temp(DIR_bams + "/raw/{wildcard}.bam")
    threads: 16
    run:
        shell("bwa mem {params.PATH_hg38} {input.R1} {input.R2} -Y -t {threads} | samtools view -hb - > {output}")

rule addRG:
    input:
        DIR_bams + "/raw/{wildcard}.bam",
    params:
        PATH_hg38=PATH_hg38,
        sample="{wildcard}",
    output:
        temp(DIR_bams + "/rg/{wildcard}.bam")
    conda: 
        "../envs/picard.yaml"
    shell:
        "picard AddOrReplaceReadGroups -Xmx32G I={input} O={output} RGID=A RGSM={params.sample} RGPL=illumina RGLB=lib1 RGPU=unit1"

rule fixmate:
    input:
        DIR_bams + "/rg/{wildcard}.bam",
    output:
        temp(DIR_bams + "/fixmate/{wildcard}.bam")
    conda: 
        "../envs/picard.yaml"
    shell:
        "picard -Xmx40g FixMateInformation I={input} O={output} SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT"

rule markDuplicates:
    input: 
        DIR_bams + "/fixmate/{wildcard}.bam"
    output:
        bam=temp(DIR_bams + "/markDuplicates/{wildcard}.bam"),
        metrics=DIR_metrics + "/markDuplicates/{wildcard}.txt"
    conda: 
        "../envs/picard.yaml"
    shell: 
        "picard -Xmx40g MarkDuplicates  I={input} O={output.bam} M={output.metrics}"

rule sort_markDuplicates:
    input:
        DIR_bams + "/markDuplicates/{wildcard}.bam"
    output:
        DIR_bams + "/sorted/{wildcard}.bam"
    shell:
        "samtools sort {input} -o {output} -m 12G -@ 12"

rule index_sorted:
    input:
        DIR_bams + "/sorted/{wildcard}.bam"
    output:
        DIR_bams + "/sorted/{wildcard}.bam.bai"
    shell:
        "samtools index {input}"