import pandas as pd
import os

# create wildcard variables for workflow
metadata = pd.read_csv("inputs/metadata.tsv", sep = "\t")   # read in metadata as a pandas dataframe
TX = metadata['tx_run_accession'].unique().tolist()         # make run accession in tx col into list
GX = metadata['gx_accession'].unique().tolist()             # make gx accession into list
KSIZES = [21, 31, 51]                                       # create a list of k-mer sizes for the workflow
gx_tmp = metadata['gx_accession'].tolist()
TX_MINUS_GX = [x + '-minus-' + y for x, y in zip(TX, gx_tmp)]   # create list that binds tx to gx

rule all:
    input: 
        expand("outputs/subtract_sourmash_sketch_describe/{tx_minus_gx}-k{ksize}.csv", tx_minus_gx = TX_MINUS_GX, ksize = KSIZES),
        expand("outputs/subtract_sourmash_sketch_filtered_describe/{tx_minus_gx}-k{ksize}.csv", tx_minus_gx = TX_MINUS_GX, ksize = KSIZES),
        expand("outputs/gx_sourmash_sketch_describe/{gx_accession}.csv", gx_accession = GX),
        expand("outputs/tx_sourmash_sketch_describe/{tx_run_accession}.csv", tx_run_accession = TX),

rule download_doryteuthis_genome:
    output: "inputs/genomes/GCA_023376005.1_genomic.fna.gz"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/023/376/005/GCA_023376005.1_UCB_Dpea_1/GCA_023376005.1_UCB_Dpea_1_genomic.fna.gz
    '''

rule download_mus_genome:
    output: "inputs/genomes/GCA_000001635.9_genomic.fna.gz"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCA_000001635.9_GRCm39/GCA_000001635.9_GRCm39_genomic.fna.gz
    '''

rule download_octopus_genome:
    output: "inputs/genomes/GCF_001194135.1_genomic.fna.gz"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/194/135/GCF_001194135.1_Octopus_bimaculoides_v2_0/GCF_001194135.1_Octopus_bimaculoides_v2_0_genomic.fna.gz
    '''
rule sourmash_sketch_gx:
    input: "inputs/genomes/{gx_accession}_genomic.fna.gz"
    output: "outputs/gx_sourmash_sketch/{gx_accession}.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,scaled=100000,abund --name {wildcards.gx_accession} -o {output} {input}
    '''

rule sourmash_sketch_tx:
    """
    Use fastq-dump to download sequencing data for a given accession.
    Pipe the sequencing data directly to a sketch without writing to disk.
    Using fastq-dump like this dumps both single end and paired end reads into std out, which can then be piped into sourmash sketch
    """
    output: "outputs/tx_sourmash_sketch/{tx_run_accession}.sig"
    conda: 'envs/sourmash.yml'
    shell:'''
    fastq-dump --disable-multithreading --fasta 0 --skip-technical --readids --read-filter pass --dumpbase --split-spot --clip -Z {wildcards.tx_run_accession} | 
        sourmash sketch dna -p k=21,k=31,k=51,scaled=100000,abund --name {wildcards.tx_run_accession} -o {output} -
    '''

rule sourmash_sig_filter_tx:
    """
    remove hashes of abundance 1 as these are highly likely to be errors
    """
    input: "outputs/tx_sourmash_sketch/{tx_run_accession}.sig"
    output: "outputs/tx_sourmash_sketch_filtered/{tx_run_accession}.sig"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash sig filter -m 2 -o {output} {input}
    '''

rule calculate_tx_not_in_gx:
    """
    Use the sourmash CLI to subtract a genome sketch from a transcriptome from the same species, retaining transcriptome abundances.
    Snakemake will auto-parse the wildcard mtx_minus_mgx on the "-minus-" string to back-propagate the correct wildcards to *x_run_accession.
    Providing the accessions in this way limits the number of ways they can combine;
    If TX solved for tx_run_accession and GX solved for gx_run_accession, then snakemake would execute this rule for all combinations of TX and GX.
    Instead, by binding pairs together in TX_MINUS_GX, only pairs of transcriptomes and genomes are subtracted.
    """
    input:
        tx_sig = 'outputs/tx_sourmash_sketch_filtered/{tx_run_accession}.sig', 
        gx_sig = 'outputs/gx_sourmash_sketch/{gx_accession}.sig', 
    output: "outputs/subtract_sourmash_sketch/{tx_run_accession}-minus-{gx_accession}-k{ksize}.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig subtract -k {wildcards.ksize} -o {output} -A {input.tx_sig} {input.tx_sig} {input.gx_sig}
    '''

rule sourmash_sig_describe_sketches_tx:
    """
    Use the sourmash CLI to report detailed information about all sketches, including number of hashes.
    Output the information as a csv file. 
    """
    input: "outputs/tx_sourmash_sketch/{tx_run_accession}.sig"
    output: "outputs/tx_sourmash_sketch_describe/{tx_run_accession}.csv"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash sig describe --csv {output} {input}
    '''

rule sourmash_sig_describe_sketches_gx:
    """
    Use the sourmash CLI to report detailed information about all sketches, including number of hashes.
    Output the information as a csv file. 
    """
    input: "outputs/gx_sourmash_sketch/{gx_accession}.sig"
    output: "outputs/gx_sourmash_sketch_describe/{gx_accession}.csv"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash sig describe --csv {output} {input}
    '''

rule sourmash_sig_describe_subtracted_sketches:
    """
    Use the sourmash CLI to report detailed information about all sketches, including number of hashes.
    Output the information as a csv file. 
    """
    input: "outputs/subtract_sourmash_sketch/{tx_minus_gx}-k{ksize}.sig"
    output: "outputs/subtract_sourmash_sketch_describe/{tx_minus_gx}-k{ksize}.csv"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash sig describe --csv {output} {input}
    '''
