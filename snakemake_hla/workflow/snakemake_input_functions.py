# input functions 
def get_WBC_path(wildcard, DIR_bams):
    paired_samples = pd.read_csv("/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test.tsv", sep="\t")
    paired_samples = paired_samples[["WBC_name", "Tumor_name"]]
    
    mask = paired_samples["Tumor_name"] == wildcard
    wbcname = paired_samples["WBC_name"][mask].tolist()[0]
    wbcpath=DIR_bams + "/sorted/" + wbcname + ".bam"
    
    return(wbcpath)

def get_WBC_chrom(wildcard):
    '''
    Given a tumor name as a wildcard without the ".bam" extension and a chromosome name, return the path to the matching wbc chromosome bam path.
    '''
    print(wildcard)
    chrom="1"
    DIR_bams="/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam"
    wildcard = re.escape(wildcard)
    # wildcard=wildcard + ".bam"
    paired_samples = pd.read_csv("/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test.tsv", sep="\t")
    paired_samples = paired_samples[["WBC_name", "Tumor_name"]]
    mask = paired_samples["Tumor_name"] == wildcard
    wbcname = paired_samples["WBC_name"][mask].tolist()[0].split(".")[0]
    wbc_path = os.path.join(DIR_bams, "split", wbcname, wbcname + ".REF_chr" + chrom +".bam")
    outlist.append(wbc_path)
    print(wildcard)
    print(wbc_path)
    print(" ")
    return(wbc_path)

def return_split_chroms_tumor(tumor_name): 
    '''
    Given a tumor name (without the .bam extension), return a list of path to the split chromosomes produced in the split_bam rule.
    '''
    DIR_bams="/groups/wyattgrp/users/amunzur/hla_pipeline/results/data/bam"
    DIR_bams + "/split/${tumor_name}/{tumor_name}.REF_chr1.bam",
    paths_list = []

    for x in ['REF_chr{}'.format(x) for x in range(1, 23)]:
        sample = ".".join([tumor_name, x, "bam"])
        sample_index = ".".join([tumor_name, x, "bam", "bai"])
        sample_path = os.path.join(DIR_bams, "split", tumor_name, sample)
        sample_index_path = os.path.join(DIR_bams, "split", tumor_name, sample_index)
        paths_list.append(sample_path)
        paths_list.append(sample_index_path)

    return(paths_list)

def return_seqz_per_chrom(tumor_name, DIR_results):
    '''
    Given a tumor name, return the paths to the seqz files generated for each chromosome.
    '''
    DIR_seqz = os.path.join(DIR_results, tumor_name)
    paths_list = []

    for x in ['REF_chr{}'.format(x) for x in range(1, 23)]:
        sample = ".".join([tumor_name, x, "seqz"])
        sample_path = os.path.join(DIR_results, "sequenza", tumor_name, sample)
        paths_list.append(sample_path)

    return(paths_list)

def get_wbc_name(tumor_name): 
    paired_samples = pd.read_csv("/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test.tsv", sep="\t")
    paired_samples = paired_samples[["WBC_name", "Tumor_name"]]

    mask = paired_samples["Tumor_name"] == tumor_name
    wbc_row = paired_samples[mask]
    if wbc_row.empty:
        raise ValueError("Could not find WBC sample for tumor '{}'".format(tumor_name))
    wbcname = wbc_row["WBC_name"].iloc[0]

    return(wbcname)

def get_WBC_hla_results(tumor_name):
    wbcname = get_wbc_name(tumor_name)
    DIR_results = "/groups/wyattgrp/users/amunzur/hla_pipeline/results"
    PATH_hla_types = os.path.join(DIR_results, "polysolver/hla_types", wbcname, "winners.hla.txt")
    return(PATH_hla_types)

# Given the name of a cfDNA sample, return the path to the combine fq file for the WBC sample.
def wbc_trimmed_combined_fq(tumor_name): 
    wbcname = get_wbc_name(tumor_name)
    DIR_results = "/groups/wyattgrp/users/amunzur/hla_pipeline/results"    
    PATH_wbc_fq = os.path.join(DIR_results, "data/fq/trimmed_combined", wbcname + ".fq")
    return(PATH_wbc_fq)

def get_read_counts(tumor_name, sample_type):
    wbcname = get_wbc_name(tumor_name)
    DIR_results = "/groups/wyattgrp/users/amunzur/hla_pipeline/results"
    DIR_readcounts = os.path.join(DIR_results, "metrics", "trimmed_combined_read_counts")
    if sample_type == "cfDNA":
        path = os.path.join(DIR_readcounts, tumor_name + ".txt")
    else: 
         path = os.path.join(DIR_readcounts, wbcname + ".txt")
    return(path)

def get_all_seqz_outputs(tumor_name, DIR_results):
    '''
    Given a sample name, return all of the separate chr.seqz.gz file paths.
    '''
    files = [os.path.join(DIR_results, "sequenza_results", tumor_name, f) for f in os.listdir(os.path.join(DIR_results, "sequenza_results", tumor_name)) if f.startswith('chr') and f.endswith('.gz')]
    return(files)

def get_input_combine_seqz(wildcard):
    seqz_dir = os.path.join(DIR_results, "sequenza_results", str(wildcard))
    return [f for f in os.listdir(seqz_dir) if f.startswith('chr') and f.endswith('.gz')]

def return_hla_fasta_path(sample_name, path_sample_list="/groups/wyattgrp/users/amunzur/hla_pipeline/resources/sample_list_panel_test.tsv"): 
    """
    Given a sample name (could be both cfDNA or WBC, doesn't matter.) return the path to the HLA fasta file. 
    Used in generating the patient HLA bams.
    """
    dir_hla_fasta = "/groups/wyattgrp/users/amunzur/hla_pipeline/resources/lohhla_fasta"  # where all fastas are
    print(sample_name)
    if "WBC" in sample_name:
        path_fasta = os.path.join(dir_hla_fasta, sample_name, "hla_fasta.fa")
    else:
        df = pd.read_csv(path_sample_list, sep="\t")
        wbc_name = df[df["Tumor_name"] == sample_name]["WBC_name"].iloc[0]
        path_fasta = os.path.join(dir_hla_fasta, wbc_name, "hla_fasta.fa")
    # path_fasta_index = path_fasta + ".fai"
    return path_fasta

def return_polysolver_output(sample_name, path_paired_samples):
    """
    """
    print(sample_name)
    if "WBC" in sample_name:
        sample_wbc=sample_name
    else:
        paired_samples=pd.read_csv(path_paired_samples, sep="\t")
        sample_wbc=paired_samples[paired_samples["Tumor_name"]==sample_name]["WBC_name"].values[0]
    
    path_sample_hla_calls_polysolver=os.path.join("/groups/wyattgrp/users/amunzur/hla_pipeline/results/polysolver/hla_types", sample_wbc, "winners.hla.txt")
    sample_hla_calls=pd.read_csv(path_sample_hla_calls_polysolver, sep="\t", header=None, names=["Gene", "Allele1", "Allele2"])
    alleles=pd.concat([sample_hla_calls["Allele1"], sample_hla_calls["Allele2"]], ignore_index=True)
    return(alleles)