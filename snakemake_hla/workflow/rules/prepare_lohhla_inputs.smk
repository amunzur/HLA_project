# wildcard is WBC bam
rule generate_hla_fasta:
    input:
        polysolver_hla="/groups/wyattgrp/users/amunzur/gillian_proj/hla-polysolver/data/abc_complete.fasta",
        hla_types=DIR_results + "/polysolver/hla_types/{wildcard}/winners.hla.txt"
    output:
        hla_fasta=DIR_resources + "/lohhla_fasta/{wildcard}/hla_fasta.fa",
        tmp_hla_fasta=DIR_temp + "/{wildcard}/winner.hla.txt",
    run:
        shell('tr "\t" "\n" < <(cat {input.hla_types}  | cut -f 2,3)> {output.tmp_hla_fasta}')
        shell("grep --no-group-separator -A 1 -f {output.tmp_hla_fasta} {input.polysolver_hla} > {output.hla_fasta}")

rule generate_hla_fasta_CDS:
    input:
        polysolver_hla="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/reference/abc_complete_CDS.fasta",
        hla_types=DIR_results + "/polysolver/hla_types/{wildcard}/winners.hla.txt"
    output:
        hla_fasta=DIR_resources + "/patient_hla_CDS/{wildcard}/hla_fasta.fa",
        tmp_hla_fasta=DIR_temp + "/{wildcard}/winners.hla.txt",
    run:
        shell('tr "\t" "\n" < <(cat {input.hla_types}  | cut -f 2,3)> {output.tmp_hla_fasta}')
        shell("grep --no-group-separator -A 1 -f {output.tmp_hla_fasta} {input.polysolver_hla} > {output.hla_fasta}")

rule index_hla_fasta:
    input:
        DIR_resources+"/lohhla_fasta/{wildcard}/hla_fasta.fa"
    output:
        DIR_resources+"/lohhla_fasta/{wildcard}/hla_fasta.fa.fai",
        DIR_resources+"/lohhla_fasta/{wildcard}/hla_fasta.fa.amb",
    shell:
        "samtools faidx {input} && bwa index {input}"

# wildcard is tumor name
rule make_symlinks: 
    input:
        PATH_tumorbam=DIR_bams + "/sorted/{wildcard}.bam",
        PATH_tumorbai=DIR_bams + "/sorted/{wildcard}.bam.bai",
    output: 
        PATH_tumorbam=DIR_results + "/data/bam/lohhla_bams/{wildcard}/{wildcard}.bam",
        PATH_tumorbai=DIR_results + "/data/bam/lohhla_bams/{wildcard}/{wildcard}.bam.bai",
    run:
        import pandas as pd
        import os
        def make_symlinks(PATH_tumor):
            tumorbam = os.path.basename(PATH_tumor)
            paired_samples = pd.read_csv("/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test.tsv", "\t", names = ["Patient_ID", "WBC_name", "Tumor_name"])
            paired_samples = paired_samples[["WBC_name", "Tumor_name"]]
            
            mask = paired_samples["Tumor_name"] == tumorbam.replace(".bam", "")
            wbcbam = paired_samples["WBC_name"][mask].tolist()[0]
            PATH_WBC = DIR_results + "/data/bam/sorted/" + wbcbam+".bam" # wbc path
            DIR_output = DIR_results + "/data/bam/lohhla_bams/" + tumorbam.replace(".bam", "") #make the dir where symlinks for wbc and tumor bams will placed 
            
            # bams
            cmd1 = " ".join(["mkdir -p", DIR_output])
            cmd2 = " ".join(["ln -s", PATH_tumor, os.path.join(DIR_output, tumorbam)])
            cmd3 = " ".join(["ln -s", PATH_WBC, os.path.join(DIR_output, wbcbam+".bam")])
            # bais
            cmd4 = " ".join(["ln -s", PATH_tumor + ".bai", os.path.join(DIR_output, tumorbam + ".bai")])
            cmd5 = " ".join(["ln -s", PATH_WBC + ".bai", os.path.join(DIR_output, wbcbam + ".bam.bai")])

            for cmd in [cmd1, cmd2, cmd3, cmd4, cmd5]:
                print(cmd)
                os.system(cmd)
        make_symlinks(input.PATH_tumorbam)

rule reformat_polysolver_winners:
    input:
        DIR_results + "/polysolver/hla_types/{wildcard}/winners.hla.txt"
    output:
        DIR_results + "/polysolver/hla_types/{wildcard}/winners.hla.LOHHLA.txt"
    shell:
        """
        awk '{{print $2; print $3}}' {input} > {output}
        """

# rule format_ploidy_and_purity:
#     input:
#         PATH_ploidy_purity="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/ct_fraction/purity_ploidy.tsv"
#     output:
#         DIR_results + "/data/ploidy_purity/{wildcard}/purity_ploidy.txt"
#     run:
#         import pandas as pd
#         purity_ploidy_df = pd.read_csv(input.PATH_ploidy_purity)
#         sample_ploidy = purity_ploidy_df[purity_ploidy_df["sample"] == '{wildcard}']["ploidy"]
#         sample_purity = purity_ploidy_df[purity_ploidy_df["sample"] == '{wildcard}']["ctdna"]
        
#         # Save to file
#         with open("example.txt", "w") as file:
#             file.write("Ploidy\ttumorPurity\ttumorPloidy\n")
#             file.write('{wildcard}')
#             file.write(sample_purity)
#             file.write(sample_ploidy)
#             file.write(sample_ploidy)

# rule modify_and_move_winners:
#     input:
#     output:
#     run:


#wildcard is wbc samples
# rule run_lohhla:
#     input:
#         tumorname="{wildcard}"
#         alternative_solutions=DIR_results + "/sequenza/{wildcard}/{wildcard}_alternative_solutions_LOHHLA.txt",
#         hla_fasta=DIR_resources + "/lohhla_fasta/{wildcard}/hla_fasta.fa",
#         lohhla_bam_dir=directory(DIR_results + "/data/bam/lohhla_bams/{wildcard}"),
#         hla_types= # winners from polysolvers
#     output:
#         directory(DIR_results + "/lohhla/{wildcard}"),
#     conda:
#         "../envs/lohhla.yaml"
#     threads: 12
#     run:
#         tumorbam = os.path.basename(PATH_tumor)
#         paired_samples = pd.read_csv("/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_table.tsv", "\t", names = ["patient_id", "normalname", "tumorname", "wbcbam", "tumorbam"])
#         paired_samples = paired_samples[["wbcbam", "tumorbam"]]

#         mask = paired_samples["tumorbam"] == tumorbam
#         wbcbam = paired_samples["wbcbam"][mask].tolist()[0]

#         cmd = 'Rscript /groups/wyattgrp/users/amunzur/gillian_proj/lohhla/LOHHLAscript.R --patientId input.tumorname --outputDir output --normalBAMfile input.normal --BAMDir input.lohhla_bam_dir --hlaPath input.hla_types --HLAfastaLoc input.hla_fasta --CopyNumLoc input.alternative_solutions --mappingStep TRUE --minCoverageFilter 10 --fishingStep TRUE --cleanUp FALSE --gatkDir /home/amunzur/anaconda3/envs/polysolver/jar --novoDir /home/amunzur/anaconda3/envs/lohhla/bin'