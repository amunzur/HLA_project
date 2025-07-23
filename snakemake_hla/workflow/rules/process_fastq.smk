rule run_fastqc_raw:
    input:
        DIR_fastq + "/raw/{wildcard}.fq.gz",
    output:
        output_zip=DIR_fastqc_raw + "/{wildcard}_fastqc.zip",
        output_html=DIR_fastqc_raw + "/{wildcard}_fastqc.html",
    threads: 5
    shell:
        "/home/amunzur/FastQC/fastqc {input} --outdir=`dirname {output.output_zip}`"

# mask low quality bases in fastq files
rule mask_fastq:
    input:
        DIR_fastq + "/raw/{wildcard}.fq.gz",
    output:
        temp(DIR_fastq + "/masked/{wildcard}.fq.gz"),
    params:
        min_base_quality=20,
    run:
        shell(
            "zcat {input} | /groups/wyattgrp/users/amunzur/gene_panel_pipeline/dependencies/fasta mask by quality - {params.min_base_quality} | gzip > {output}"
        )

rule trim_fastq:
    input:
        R1=DIR_fastq + "/masked/{wildcard}_1.fq.gz",
        R2=DIR_fastq + "/masked/{wildcard}_2.fq.gz",
    output:
        R1=DIR_fastq + "/trimmed/{wildcard}_1.fq.gz",
        R2=DIR_fastq + "/trimmed/{wildcard}_2.fq.gz",
        html_report="results/reports/fastp/{wildcard}.html",
        json_report="results/reports/fastp/{wildcard}.json",
    threads: 12
    params:
        minimum_read_length=50,
    shell:
        """
        /groups/wyattgrp/cargo/bin/fasta interleave <(zcat {input.R1} | sed 's-/1 - -g') <(zcat {input.R2} | sed 's-/2 - -g') |\
        fastp --stdin --interleaved_in -o {output.R1} -O {output.R2} \
        --length_required {params.minimum_read_length} \
        --disable_trim_poly_g \
        --disable_quality_filtering \
        --dont_overwrite \
        --html {output.html_report} \
        --json {output.json_report} \
        --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --umi --umi_loc per_read --umi_len 3 --umi_skip 2 \
        --thread {threads}
        """

# to run DASH later
rule combine_trimmed_fq:
    input:
        R1=DIR_fastq + "/trimmed/{wildcard}_1.fq.gz",
        R2=DIR_fastq + "/trimmed/{wildcard}_2.fq.gz"
    output:
        DIR_fastq + "/trimmed_combined/{wildcard}.fq.gz"
    shell:
        """
        zcat {input.R1} {input.R2} | gzip > {output}
        """
