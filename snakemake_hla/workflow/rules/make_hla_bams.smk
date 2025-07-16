rule align_to_patient_specific_hla_reference:
    input:
        R1=DIR_fastq + "/trimmed/{wildcard}_1.fq.gz",
        R2=DIR_fastq + "/trimmed/{wildcard}_2.fq.gz",
        patient_hla_fasta = lambda wildcards: return_hla_fasta_path(wildcards.wildcard, path_sample_list="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test.tsv"),
    output:
        DIR_bams + "/hla/RAW/{wildcard}.bam",
    threads: 4
    shell:
        """
        bwa index {input.patient_hla_fasta} && bwa mem {input.patient_hla_fasta} {input.R1} {input.R2} -Y -t {threads} | samtools view -b -h -o {output} - > {output}        
        """

# Add read groups to the mapped bam file
rule addRG_hla:
    input:
        DIR_bams + "/hla/RAW/{wildcard}.bam"
    params:
        sample="{wildcard}",
    output:
        temp(DIR_bams + "/hla/RG/{wildcard}.bam")
    threads: 12
    conda: 
        "../envs/picard.yaml"
    shell:
        "picard AddOrReplaceReadGroups -Xmx20G I={input} O={output} RGID=A RGSM={params.sample} RGPL=illumina RGLB=lib1 RGPU=unit1"

rule fixmate_hla:
    input:
        DIR_bams + "/hla/RG/{wildcard}.bam"
    output:
        temp(DIR_bams + "/hla/fixmate/{wildcard}.bam")
    threads: 12
    conda: 
        "../envs/picard.yaml"
    shell:
        "picard -Xmx40g FixMateInformation I={input} O={output} SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT"

rule markDuplicates_hla:
    input: 
        DIR_bams + "/hla/fixmate/{wildcard}.bam"
    output:
        bam=temp(DIR_bams + "/hla/markdup/{wildcard}.bam"),
        metrics=DIR_metrics + "/markDuplicates_hla_bams/{wildcard}.txt"
    conda: 
        "../envs/picard.yaml"
    shell: 
        "picard -Xmx40g MarkDuplicates  I={input} O={output.bam} M={output.metrics}"

# We save a copy of the hla bam with all the reads intact
rule sort_markDuplicates_hla_full:
    input:
        DIR_bams + "/hla/markdup/{wildcard}.bam"
    output:
        DIR_bams + "/sorted_hla_full/{wildcard}.bam"
    threads: 12
    shell:
        """
        samtools sort {input} -o {output}
        """

rule index_hla_bams:
    input:
        DIR_bams + "/sorted_hla_full/{wildcard}.bam"
    output:
        DIR_bams + "/sorted_hla_full/{wildcard}.bam.bai"
    run:
        shell(
            "samtools index {input}"
        )

# If a read is more than 20% soft clipped we exclude it. Returns a filtered bam. 
rule filter_soft_clipping:
    input:
        DIR_bams + "/hla/markdup/{wildcard}.bam"
    output:
        temp(DIR_bams + "/hla/filtered_softclipping/{wildcard}_filtered_softclipping.bam")  # Output filtered BAM file
    threads: 12
    run:
        shell(r"""
            samtools view -h {input} | \
            awk 'BEGIN {{OFS="\\t"}} 
                function calculate_softclip_fraction(cigar) {{
                    soft_len = total_len = 0;
                    while (match(cigar, /[0-9]+[MIDNSHP=X]/)) {{
                        len = substr(cigar, RSTART, RLENGTH-1);
                        type = substr(cigar, RLENGTH, 1);
                        cigar = substr(cigar, RSTART + RLENGTH);
                        if (type == "S") {{ soft_len += len }}
                        if (type ~ /[M=X]/) {{ total_len += len }}
                    }}
                    if (total_len > 0) {{
                        return soft_len / (soft_len + total_len);  # Return fraction of soft-clipped bases
                    }} else {{
                        return 0;  # Return 0 if no aligned portion
                    }}
                }}
                {{
                    if (/^@/) {{
                        print $0;  # Print header lines
                    }} else {{
                        cigar = $6;  # Get CIGAR string from column 6
                        softclip_fraction = calculate_softclip_fraction(cigar);  # Calculate soft-clip fraction
                        if (softclip_fraction <= 0.2) {{
                            print $0;  # Keep read if soft-clipping is <= 20%
                        }}
                    }}
                }}' | samtools view -Sb - > {output}
        """)

# Excludes reads with more than 2 mismatches. Returns filtered bam without those reads.
rule filter_mismatches:
    input:
        DIR_bams + "/hla/filtered_softclipping/{wildcard}_filtered_softclipping.bam"  # Input filtered BAM file
    output:
        temp(DIR_bams + "/hla/filtered_mismatch/{wildcard}_filtered_mismatch.bam")  # Output filtered BAM file
    threads: 12
    shell:
        r"""
        samtools view -h {input} | \
        awk 'BEGIN {{OFS="\\t"}} 
            {{
                if (/^@/) {{
                    print $0;  # Print header lines
                }} else {{
                    for (i = 12; i <= NF; i++) {{
                        if ($i ~ /^NM:i:/) {{
                            nm = substr($i, 6) + 0;  # Extract NM tag value and convert to number
                            if (nm <= 2) {{
                                print $0;  # Print read if it has 2 or fewer mismatches
                            }}
                            break;
                        }}
                    }}
                }}
            }}' | samtools view -Sb - > {output}
        """

# After filtering we sort the outputted bams
rule sort_markDuplicates_hla:
    input:
        DIR_bams + "/hla/filtered_mismatch/{wildcard}_filtered_mismatch.bam"
    output:
        DIR_bams + "/sorted_hla_filtered/{wildcard}.bam"
    threads: 12
    shell:
        """
        samtools sort {input} -o {output}
        """

rule index_hla_bams_filtered:
    input:
        DIR_bams + "/sorted_hla_filtered/{wildcard}.bam"
    output:
        DIR_bams + "/sorted_hla_filtered/{wildcard}.bam.bai"
    run:
        shell(
            "samtools index {input}"
        )