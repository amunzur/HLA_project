rule depth_exome:
    input:
        DIR_bams + "/sorted/{wildcard}.bam",
    output:
        DIR_metrics + "/depth_exome/{wildcard}.txt",
    params:
        PATH_bed="/groups/wyattgrp/users/amunzur/pipeline/resources/panel/hyperexome/KAPA_HyperExome_hg38_capture_targets.bed"
    run:
        shell(
            "samtools depth -b {params.PATH_bed} -a {input} > {output}"
            )

rule depth_exome_MEDIAN:
    input:
        DIR_metrics + "/depth_exome/{wildcard}.txt",
    output:
        DIR_metrics + "/depth_exome_MEDIAN/{wildcard}.txt",
    shell:
        "awk '{{a[NR]=$3}} END {{n=asort(a); printf(\"10th percentile: %s\\nMedian: %s\\n90th percentile: %s\\n\", a[int(n*0.1)], a[int(n/2)], a[int(n*0.9)])}}' {input} > {output}"

###########################################################
# Depth specifically in hla regions
rule depth_panel_hla:
    input:
        hla_bam=DIR_bams + "/sorted_hla_full/{wildcard}.bam",
        idxstats=DIR_metrics + "/idxstats/{wildcard}.txt"
    output:
        hla_a=DIR_metrics + "/depth_hla/{wildcard}_A.txt",
        hla_b=DIR_metrics + "/depth_hla/{wildcard}_B.txt",
        hla_c=DIR_metrics + "/depth_hla/{wildcard}_C.txt",
    shell:
        """
        # Extract HLA chromosome names from the idxstats file
        hla_chromosomes=$(grep -E "^hla_" {input.idxstats} | cut -f1)
        
        # Group chromosomes by type
        hla_a_chromosomes=$(echo "$hla_chromosomes" | grep "^hla_a")
        hla_b_chromosomes=$(echo "$hla_chromosomes" | grep "^hla_b")
        hla_c_chromosomes=$(echo "$hla_chromosomes" | grep "^hla_c")
        
        # Process each HLA group
        for group in hla_a hla_b hla_c; do
            # Set the appropriate output file based on the group
            if [ "$group" = "hla_a" ]; then
                output_file={output.hla_a}
                chromosomes="$hla_a_chromosomes"
            elif [ "$group" = "hla_b" ]; then
                output_file={output.hla_b}
                chromosomes="$hla_b_chromosomes"
            elif [ "$group" = "hla_c" ]; then
                output_file={output.hla_c}
                chromosomes="$hla_c_chromosomes"
            fi
            
            # Clear the output file if it exists
            > "$output_file"
            
            # Loop through the chromosomes and calculate depth
            for chr in $chromosomes; do
                samtools depth -r "$chr" {input.hla_bam} >> "$output_file"
            done
        done
        """

rule depth_panel_hla_MEDIAN:
    input:
        DIR_metrics + "/depth_hla/{wildcard}_{gene}.txt",
    output:
        DIR_metrics + "/depth_hla_MEDIAN/{wildcard}_{gene}.txt",
    shell:
        "awk '{{a[NR]=$3}} END {{n=asort(a); printf(\"10th percentile: %s\\nMedian: %s\\n90th percentile: %s\\n\", a[int(n*0.1)], a[int(n/2)], a[int(n*0.9)])}}' {input} > {output}"

###########################################################
# Calculates depth in all probes
rule run_depth_in_probes:
    input:
        DIR_bams + "/sorted/{wildcard}.bam",
    output:
        DIR_metrics + "/depth_probes/{wildcard}.txt",
    params:
        PATH_bed="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/hla_panel/probe_locations.bed"
    run:
        shell(
            "samtools depth -b {params.PATH_bed} -a {input} > {output}"
            )
###########################################################

# Calculates depth in nonHLA panel target regions, this is computed to calculate the median nonHLA depth so hg38 bam is used.
rule depth_panel_nonHLA:
    input:
        DIR_bams + "/sorted/{wildcard}.bam",
    output:
        DIR_metrics + "/depth_panel_nonHLA/{wildcard}.txt",
    params:
        PATH_bed="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/hla_panel/target_regions_nonHLA.bed"
    run:
        shell(
            "samtools depth -b {params.PATH_bed} -a {input} > {output}"
            )

rule depth_panel_nonHLA_MEDIAN:
    input:
        DIR_metrics + "/depth_panel_nonHLA/{wildcard}.txt",
    output:
        DIR_metrics + "/depth_panel_nonHLA_MEDIAN/{wildcard}.txt",
    shell:
        "awk '{{a[NR]=$3}} END {{n=asort(a); printf(\"10th percentile: %s\\nMedian: %s\\n90th percentile: %s\\n\", a[int(n*0.1)], a[int(n/2)], a[int(n*0.9)])}}' {input} > {output}"

################################################
rule trimmed_read_counts:
    input:
        DIR_fastq + "/trimmed/{wildcard}_1.fq.gz"
    output:
        DIR_metrics + "/trimmed_read_counts/{wildcard}.txt"
    shell:
        """
        n_raw_reads=$(wc -l {input} | awk '{{print $1}}')
        n_reads=$((n_raw_reads / 4))
        echo $n_reads > {output}
        """

rule raw_read_counts:
    input:
        DIR_fastq + "/raw/{wildcard}_1.fq.gz"
    output:
        DIR_metrics + "/raw_fastq_read_counts/{wildcard}.txt"
    shell:
        """
        n_raw_reads=$(zcat {input} | wc -l)
        n_reads=$((n_raw_reads / 4))
        echo "{wildcards.wildcard}\t$((n_reads * 2))" > {output}
        """

rule save_dup_perc:
    input:
        DIR_metrics + "/markDuplicates/{wildcard}.txt"
    output:
        DIR_metrics + "/duplicate_fraction/{wildcard}.txt"
    shell:
        """
        dup_perc=$(grep "lib1" {input} | cut -f9)
        echo "{wildcards.wildcard}\t$dup_perc" > {output}
        """

# to run dash later on using these read counts
rule trimmed_combined_read_counts:
    input:
        DIR_fastq + "/trimmed_combined/{wildcard}.fq.gz"
    output:
        DIR_metrics + "/trimmed_combined_read_counts/{wildcard}.txt"
    shell:
        """
        n_raw_reads=$(zcat {input} | wc -l | awk '{{print $1}}')
        n_reads=$((n_raw_reads / 4))
        echo $n_reads > {output}
        """

rule run_idxstats:
    input:
        bam=DIR_bams + "/sorted_hla_full/{wildcard}.bam",
        bam_index=DIR_bams + "/sorted_hla_full/{wildcard}.bam.bai",
    output:
        DIR_metrics + "/idxstats/{wildcard}.txt"
    shell:
        """
        samtools idxstats {input.bam} > {output}
        """

# rule run_hla_insert_size:
#     input:
#         DIR_bams + "/subsetted_sorted_hla/{wildcard}_HLA_{gene}.bam",
#     output:
#         hist = DIR_results + "/metrics/hla_fragment/HLA_{gene}_{wildcard}.png",
#         metrics = DIR_results + "/metrics/hla_fragment/HLA_{gene}_{wildcard}_metrics.txt",
#     threads: 12
#     conda:
#          "../envs/deeptools.yaml"
#     shell:
#         """
#         bamPEFragmentSize -b {input} \
#             --histogram {output.hist} \
#             --plotTitle "Fragment size" \
#             --maxFragmentLength 1000 \
#             --samplesLabel "{wildcards.wildcard}" > {output.metrics}
#         """

# rule run_WES_insert_size:
#     input:
#         DIR_bams + "/subsetted_sorted/{wildcard}.bam",
#     output:
#         hist = DIR_results + "/metrics/WES_fragment/{wildcard}.png",
#         metrics = DIR_results + "/metrics/WES_fragment/{wildcard}_metrics.txt",
#     threads: 12
#     conda:
#          "../envs/deeptools.yaml"
#     shell:
#         """
#         bamPEFragmentSize -b {input} \
#             --histogram {output.hist} \
#             --plotTitle "Fragment size" \
#             --maxFragmentLength 1000 \
#             --samplesLabel "{wildcards.wildcard}" > {output.metrics}
#         """