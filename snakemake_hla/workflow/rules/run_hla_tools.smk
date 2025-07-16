# wildcard must be the WBC sample
rule call_hla_type:
    input:
        DIR_bams + "/sorted/{wildcard}.bam",
    output:
        DIR_results + "/polysolver/hla_types/{wildcard}/winners.hla.txt", # the main output file needed
    params:
        DIR_output=DIR_results + "/polysolver/hla_types/{wildcard}", # output dir to feed into the command
    threads: 12
    conda: 
        "../envs/polysolver.yaml"
    shell:
            "/home/amunzur/anaconda3/envs/polysolver/bin/shell_call_hla_type \
            {input} \
            Unknown \
            1 \
            hg38 \
            STDFQ \
            0 \
            {params}"

# wildcard is tumor bam
rule call_hla_mutations_from_type:
    input:
        PATH_tumor=DIR_bams + "/sorted/{wildcard}.bam",
        PATH_hla_results=lambda wildcards: get_WBC_hla_results("{wildcard}".format(wildcard=wildcards.wildcard))
    output:
        DIR_results + "/polysolver/hla_mutations/{wildcard}/hla_mutations/hla.intervals",
        # directory(DIR_results + "/polysolver/hla_mutations/{wildcard}/hla_mutations") # dont remove, needed to run DASH.
    params:
        NAME_tumor="{wildcard}",
        DIR_output=DIR_results + "/polysolver/hla_mutations/{wildcard}/hla_mutations"
    conda: 
        "../envs/polysolver.yaml"
    threads: 12
    shell:
        """
        NAME_WBC=$(cat /groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test.tsv | cut -f2,3 | grep {params.NAME_tumor} | cut -f1).bam
        PATH_WBC=/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam/sorted/$NAME_WBC
        shell_call_hla_mutations_from_type ${{PATH_WBC}} {input.PATH_tumor} {input.PATH_hla_results} hg38 STDFQ {params.DIR_output}
        """

# wildcard is tumor bam
rule annotate_mutations:
    input:
        DIR_results + "/polysolver/hla_mutations/{wildcard}/hla_mutations/hla.intervals"
    output:
        DIR_results + "/polysolver/hla_mutations/{wildcard}/hla_mutations/{wildcard}.mutect.unfiltered.annotated",
        DIR_results + "/polysolver/hla_mutations/{wildcard}/hla_mutations/{wildcard}.strelka_indels.unfiltered.annotated",
    params:
        indiv="{wildcard}",
        DIR_output=DIR_results + "/polysolver/hla_mutations/{wildcard}/hla_mutations"
    conda: 
        "../envs/polysolver.yaml"
    threads: 12
    shell:
        """
        shell_annotate_hla_mutations {params.indiv} {params.DIR_output}
        """